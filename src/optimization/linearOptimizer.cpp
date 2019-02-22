//
//  linearOptimizer.cpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include <iomanip>

#include "armadillo"
#include "linearOptimizer.hpp"
#include "../coronagraph/coronagraph.hpp"

double linearOptimizer::eval_mu_contrast(const std::vector<double> &x, std::vector<double> &grad, void *parentPointer)
{
    linearOptimizer* thisOptimizer = reinterpret_cast<linearOptimizer*>(parentPointer);

    arma::vec optVec0 = thisOptimizer->optVec;
    arma::vec curOptVec;
    arma::vec curDeltaOptVec;

    int nJacRows = thisOptimizer->nJacRows;
    int n_elem = thisOptimizer->optVec.n_elem;

    thisOptimizer->objMat.submat(nJacRows, 0, nJacRows + n_elem-1, n_elem-1) = x[0]*arma::eye<arma::mat>(n_elem, n_elem);
    arma::solve(curDeltaOptVec, thisOptimizer->objMat, thisOptimizer->objVec);
    curOptVec = optVec0 + curDeltaOptVec;
    globalCoronagraph->set_optimization_data(thisOptimizer->componentName, thisOptimizer->dataName, curOptVec);
    thisOptimizer->compute_contrast();
    std::cout << x[0] << " " << thisOptimizer->contrast << std::endl;
    
    return thisOptimizer->contrast;
}

linearOptimizer::linearOptimizer() {
}

linearOptimizer::linearOptimizer(initCommandSet*& cmdBlocks) {
    init(cmdBlocks);
}

void linearOptimizer::init(initCommandSet*& cmdBlock) {
    
    std::vector<initCommandSet*> subBlocks = cmdBlock->find_command_blocks();
    
    for (int j=0; j<subBlocks.size(); j++) {
        if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "linearOptimizer")) {
            init_block(subBlocks[j]);
        }
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "regionContrast")) {
            region = new regionContrast(subBlocks[j]);
            region->print("region");
        }
    }
    print("initialized nonlinear optimizer");
}

void linearOptimizer::init_block(initCommandSet*& cmdBlock) {
    std::vector<initCommandSet*> subBlocks = cmdBlock->find_command_blocks();
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
}

void linearOptimizer::set(std::string fieldName, const char *arg) {
//    std::cout << "fieldname: " << fieldName << ", arg: " << arg << std::endl;
    if (fieldName == "linearOptimizer")
        ;
    else if (fieldName == "dataToOptimize") {
        // arg is two strings separated by a comma,
        // the first string is the name of the component to be optimized
        // the second string is the name of the data to be optimized
        const char *cPtr = strchr(arg, ',');
        int componentNameLen = cPtr - arg;
        int dataNameLen = strlen(arg) - componentNameLen + 1;
        componentName = new char[componentNameLen + 1];
        dataName = new char[dataNameLen + 1];
        strncpy(componentName, arg, componentNameLen);
        componentName[componentNameLen] = '\0';
        strncpy(dataName, cPtr + 1, dataNameLen);
        dataName[dataNameLen] = '\0';
        std::cout << "optimize component " << componentName << ", data " << dataName << std::endl;
        //        draw("initial FPM");
        
        globalCoronagraph->get_optimization_data(componentName, dataName, startVec);
        optVecSize = startVec.n_elem;
        optVec.set_size(optVecSize);
    }
    else if (fieldName == "initOptVec") {
        startVec.set_size(optVecSize);
        if (!strcmp(arg, "random")) {
            startVec = globalUpperBound*(arma::randu(optVecSize) - 0.5);
            std::cout << "mean(startVec) = " << mean(startVec) << " min(startVec) = " << min(startVec) << " max(startVec) = " << max(startVec) << std::endl;
        } else {
            float v1;
            // arg is a comma separated pair of doubles
            sscanf(arg, "%f", &v1);
            startVec = v1*arma::ones<arma::vec>(optVecSize);
        }
    }
    else if (fieldName == "contrastRelTolerance") {
        // arg is a single float
        contrastRelTolerance = atof(arg);
    }
    else if (fieldName == "jacobianDx") {
        // arg is a single float
        jacobianDx = atof(arg);
    }
    else if (fieldName == "globalBounds") {
        float v1, v2;
        // arg is a comma separated pair of doubles
        sscanf(arg, "%f, %f", &v1, &v2);
        globalLowerBound = v1;
        globalUpperBound = v2;
    }
    // set stopping criteria
    else if (fieldName == "maxIterations") {
        // arg is a single int
        maxIterations = atoi(arg);
    }
    else if (fieldName == "stopContrast") {
        // arg is a single float
        stopContrastValue = atof(arg);
    }
    else if (fieldName == "deltaOptVecTolerance") {
        // arg is a single float
        deltaOptVecTolerance = atof(arg);
    }
    else if (fieldName == "regularizationMu") {
        // arg is a single float
        regularizationMu = atof(arg);
    }
    else if (fieldName == "stepSizeScale") {
        // arg is a single float
        stepSizeScale = atof(arg);
    }
    // set output parameters
    else if (fieldName == "outputInterval") {
        // arg is a single int
        outputInterval = atoi(arg);
    }
    else if (fieldName == "outputDirectory") {
        // arg is a string
        outputDirectory = new char[strlen(arg)+1];
        strcpy(outputDirectory, arg);
    }
    else if (fieldName == "saveState") {
        // arg is a single int
        saveState = atoi(arg);
    }
    else
        std::cout << "!!! linearOptimizer bad set field name: " << fieldName << std::endl;
    
}

void linearOptimizer::reset_Efield(efield *& E) {
    delete E;
    E = new efield(*initialEfield);
    E->set("name", "FPM Efield");
}

void linearOptimizer::optimize(void) {
    arma::wall_clock timer;
    
    calibEfield = new efield(*initialEfield);
    calibEfield->set("name", "calibration Efield");
    fullEfield = new efield(*initialEfield);
    fullEfield->set("name", "FPM Efield");
    
    arma::mat calibIntensity;
    calibMaxIntensity = region->compute_intensity(calibEfield, calibIntensity, true);
    save_mat("calibration_PSF.fits", calibIntensity);
    std::cout << "calibMaxIntensity = " << calibMaxIntensity << std::endl;
    
    region->get_region_pixels(calibEfield, regionPixelIndex);
    // set the objective matrix size = 2*# of polarizations*# of sources*# of wavelengths by # of points in the region of interest
    // + the # of optimized data points for the diagonal regularization matrix
    int nJacRows = 2*calibEfield->E.size()*calibEfield->E[0].size()*calibEfield->E[0][0]->n_slices*regionPixelIndex.n_elem;
    objMat.set_size(nJacRows + optVec.n_elem, optVec.n_elem);
//    std::cout << "size(objMat) = " << size(objMat)  << ", optVec.n_elem = " << optVec.n_elem << std::endl;
//    std::cout << "nJacRows = " << nJacRows << ", nJacRows + optVec.n_elem-1 = " << nJacRows + optVec.n_elem-1 << std::endl;
    objMat.submat(nJacRows, 0, nJacRows + optVec.n_elem-1, optVec.n_elem-1) = regularizationMu*arma::eye<arma::mat>(optVec.n_elem, optVec.n_elem);

    objVec.set_size(nJacRows + optVec.n_elem);
//    std::cout << "size(objVec) = " << size(objVec)  << ", optVec.n_elem = " << optVec.n_elem << std::endl;
//    std::cout << "nJacRows = " << nJacRows << ", nJacRows + optVec.n_elem-1 = " << nJacRows + optVec.n_elem-1 << std::endl;
    objVec.subvec(nJacRows, nJacRows+optVec.n_elem-1) = arma::zeros<arma::vec>(optVec.n_elem);
    
    std::cout << "============= starting optimization =============" << std::endl;
    std::string fname = (std::string)outputDirectory + "/optBestVal.txt";
    bestValFid = fopen(fname.c_str(), "w");
    
    fname = (std::string)outputDirectory + "/optValHistory.txt";
    histValFid = fopen(fname.c_str(), "w");
    
    optVec = startVec;
    deltaOptVec = arma::ones(size(optVec));
    lastDeltaOptVec = optVec;
    optimizationTimer.tic();
    // perform the optimization
    timer.tic();
    while (!stopping_criteria()) {
//        arma::wall_clock setOpttimer;
//        setOpttimer.tic();
        globalCoronagraph->set_optimization_data(componentName, dataName, optVec);
//        std::cout << "time to set optimization data: " << setOpttimer.toc() << std::endl;

        compute_contrast();
        
        // output data from this step
        
        fprintf(get_histValFid(), "%g %g\n", optimizationTimer.toc(), contrast);
        fflush(get_histValFid());
        if (contrast < get_bestVal()) {
            fprintf(get_bestValFid(), "%g %g\n", optimizationTimer.toc(), contrast);
            fflush(get_bestValFid());
            set_bestVal(contrast);
            globalCoronagraph->save_optimization_data(componentName, dataName, get_ouput_directory());
        }
        if (fabs(lastContrast - contrast) < contrastRelTolerance)
            contrastNotChanging = true;
        if (contrast < stopContrastValue)
            contrastReachedStopValue = true;
        lastContrast = contrast;

        if (norm(abs(deltaOptVec)) < deltaOptVecTolerance)
            optVecNotChanging = true;
        

        std::cout << "evaluation " << nIterations << ", " << timer.toc()/60 << " minutes, contrast = " << contrast << std::endl;
        timer.tic();
        if (!stopping_criteria()) {
            // solve for the change in optimization parameters
            bool iterationStatus = false;
//            // test contrast by taking a single step with the old contrast
//            arma::vec testOptVec = optVec + lastDeltaOptVec;
//            globalCoronagraph->set_optimization_data(componentName, dataName, &testOptVec);
//
//            reset_Efield(calibEfield);
//            calibMaxIntensity = region->compute_intensity(calibEfield, calibIntensity, true);
//            std::cout << "calibMaxIntensity = " << calibMaxIntensity << std::endl;
//            reset_Efield(fullEfield);
//            fullMaxIntensity = region->compute_intensity(fullEfield, fullIntensity, false)/calibMaxIntensity;
//            fullIntensity /= calibMaxIntensity;
//            double testContrast = arma::mean(fullIntensity(regionPixelIndex));
//            // if contrast decreases by more than 10%, or contrast increases, recompute the Jacobian
//            std::cout << "0.1*contrast = " << 0.1*contrast << ", lastContrast - testContrast = " << lastContrast - testContrast << std::endl;
//            if (lastContrast - testContrast > 0.1*contrast | lastContrast - testContrast < 0)  {
//                iterationStatus = compute_step();
//                lastDeltaOptVec = deltaOptVec;
//            } else {
//                iterationStatus = true;
//                deltaOptVec = lastDeltaOptVec;
//            }
            iterationStatus = compute_step();
            if (iterationStatus) {
//                std::cout << "optVec: " << optVec << std::endl;
//                std::cout << "deltaOptVec: " << deltaOptVec << std::endl;
                
                // this step size thing didn't work
//                double maxDelta = max(abs(deltaOptVec));
//                if (maxDelta > stepSizeScale*globalUpperBound)
//                    optVec += stepSizeScale*globalUpperBound*deltaOptVec/maxDelta;
//                else
//                    optVec += deltaOptVec;
                
                char outName[200];
                sprintf(outName, "%s/deltaOptVec_%d.fits", get_ouput_directory(), nIterations);
                save_vec(outName, deltaOptVec);
                
                optVec += deltaOptVec;
                // enforce bounds
                for (int i=0; i<optVec.n_elem; i++) {
                    if (optVec[i] > globalUpperBound) {
                        std::cout << "enforcing upper bound" << std::endl;
                        optVec[i] = globalUpperBound;
                    } else if (optVec[i] < globalLowerBound) {
                        std::cout << "enforcing lower bound" << std::endl;
                        optVec[i] = globalLowerBound;
                    }
                }
            } else
                std::cout << "solver failed" << std::endl;
            
            std::cout << "norm deltaOptVec = " << norm(abs(deltaOptVec)) << ", max(abs(deltaOptVec)) = " << max(abs(deltaOptVec))  << ", mean(abs(deltaOptVec)) = " << mean(abs(deltaOptVec)) << std::endl;
        }
        nIterations++;
    }
    fclose(bestValFid);
    fclose(histValFid);
    
    std::cout << "nIterations: " << nIterations << ", time " << timer.toc() << " seconds" << std::endl;
}

bool linearOptimizer::stopping_criteria(void) {
    if (contrastNotChanging)
        std::cout << "stopping because of contrastNotChanging" << std::endl;
    if (contrastReachedStopValue)
        std::cout << "stopping because of contrastReachedStopValue" << std::endl;
    if (optVecNotChanging)
        std::cout << "stopping because of optVecNotChanging" << std::endl;
    if (nIterations > maxIterations)
        std::cout << "stopping because of nIterations > maxIterations" << std::endl;
    
    return (contrastNotChanging | contrastReachedStopValue | optVecNotChanging | nIterations > maxIterations);
}

bool linearOptimizer::compute_step(void)
{
    arma::wall_clock timer;
    bool testMode = false;
    
    // compute the region intensity via the regionContrast object, which calls the coronagraph
    //    std::cout << "calibMaxIntensity = " << calibMaxIntensity << std::endl;
    timer.tic();
    
    //
    //
    // Compte the jacobian
    //
    
    // if the error is too big or if the contrast is no longer decreasing, compute the gradient
    if (true) {
        compute_jacobian();
//        load_mat("objMat_part3.fits", objMat);
    } else {
        std::cout << "using previous objMat" << std::endl;
        objMat = lastObjMat;
    }
    
    compute_objective_vector();
    save_vec("objVec.fits", objVec);
//    load_vec("objVec_part3.fits", objVec);
/*
    nlopt::opt opt(nlopt::LN_BOBYQA, 1);
    opt.set_min_objective(linearOptimizer::eval_mu_contrast, this);
    opt.set_upper_bounds(1e9);
    opt.set_lower_bounds(1e7);
    opt.set_ftol_rel(1e-10);
    opt.set_maxeval(10);
    std::vector<double>x(1,regularizationMu);
//    opt.set_initial_step(1e7);
    double minf;
    try{
        nlopt::result result = opt.optimize(x, minf);
        std::cout << "found minimum contrast " << std::setprecision(10) << minf << std::endl;
        std::cout << "result " << result << std::endl;
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    regularizationMu = x[0];
    
    int nJacRows = nJacRows;
    int n_elem = optVec.n_elem;
    objMat.submat(nJacRows, 0, nJacRows + n_elem-1, n_elem-1) = x[0]*arma::eye<arma::mat>(n_elem, n_elem);
*/

    int nJacRows = 2*calibEfield->E.size()*calibEfield->E[0].size()*calibEfield->E[0][0]->n_slices*regionPixelIndex.n_elem;
    arma::vec optVec0 = optVec;
    arma::vec curOptVec;
    double lastContrast = 1e6;
    double lastRm = 0;
    for (double rm = 1e7; rm <= 3e8; rm += 1e7) {
        objMat.submat(nJacRows, 0, nJacRows + optVec.n_elem-1, optVec.n_elem-1) = rm*arma::eye<arma::mat>(optVec.n_elem, optVec.n_elem);
        arma::solve(deltaOptVec, objMat, objVec);
        curOptVec = optVec0 + deltaOptVec;
        globalCoronagraph->set_optimization_data(componentName, dataName, curOptVec);
        compute_contrast();
        std::cout << rm << " " << contrast << std::endl;
        if (contrast < lastContrast) {
            lastContrast = contrast;
            lastRm = rm;
        }
        else {
            contrast = lastContrast;
            break;
        }
    }
    regularizationMu = lastRm;
    objMat.submat(nJacRows, 0, nJacRows + optVec.n_elem-1, optVec.n_elem-1) = regularizationMu*arma::eye<arma::mat>(optVec.n_elem, optVec.n_elem);
    std::cout << "selected contrast = " << contrast << " at mu = " << regularizationMu << std::endl;
    objMat.save("objMat.txt", arma::arma_ascii);
    save_mat("objMat.fits", objMat);
    
    lastObjMat = objMat;
    lastObjMatContrast = contrast;

    // solve for the new optimization parameters
    bool solveStatus = arma::solve(deltaOptVec, objMat, objVec);
    
    return solveStatus;
}

void linearOptimizer::compute_contrast(void)
{
    arma::mat calibIntensity;

    reset_Efield(calibEfield);
    calibMaxIntensity = region->compute_intensity(calibEfield, calibIntensity, true);
//    std::cout << "calibMaxIntensity = " << calibMaxIntensity << std::endl;

    arma::mat fullIntensity;
    reset_Efield(fullEfield);
    double fullMaxIntensity = region->compute_intensity(fullEfield, fullIntensity, false)/calibMaxIntensity;
    // normalize
    fullIntensity /= calibMaxIntensity;
    //    draw_mat(log10(fullIntensity), "fullIntensity", "matlab");
    // compute the contrast
    contrast = arma::mean(fullIntensity(regionPixelIndex));
}

void linearOptimizer::compute_jacobian(void)
{
    
    //
    //
    // Compte the jacobian
    //
    
    arma::cx_mat tmpE(fullEfield->E[0][0]->n_rows, fullEfield->E[0][0]->n_cols);
    std::cout << "computing jacobian: ";
    arma::vec dx;
    
    //        save_mat("fullEfield.fits", fullEfield->E[0][0]->slice(0));
    //        char outName[200];
    for (int i=0; i<optVec.n_elem; i++) {
        dx = optVec;
        dx(i) = dx(i) + jacobianDx;
        if (!(i%10) | optVec.n_elem <= 20) {
            std::cout << i << " ";
            fflush(stdout);
        }
        globalCoronagraph->set_optimization_data(componentName, dataName, dx);
        reset_Efield(jacEfield);
        
        globalCoronagraph->execute(jacEfield, 0);
        //            sprintf(outName, "jacEfield_%d.fits", i);
        //            save_mat(outName, jacEfield->E[0][0]->slice(0));
        // distribute the result
        // each i is a column of the jac matrix
        // loop over polarizations, sources and wavelengths
        // each column i contains:
        // Re(dE/dx_i) for source 0, polarization 0, wavelength 0
        // Im(dE/dx_i) for source 0, polarization 0, wavelength 0
        // ...
        // Re(dE/dx_i) for source 0, polarization 0, wavelength calibEfield->E[0][0]->n_slices - 1
        // Im(dE/dx_i) for source 0, polarization 0, wavelength calibEfield->E[0][0]->n_slices - 1
        // ... ...
        // Re(dE/dx_i) for source 0, polarization E->E[0].size(), wavelength 0
        // Im(dE/dx_i) for source 0, polarization E->E[0].size(), wavelength 0
        // ...
        // Re(dE/dx_i) for source 0, polarization E->E[0].size(), wavelength calibEfield->E[0][0]->n_slices - 1
        // Im(dE/dx_i) for source 0, polarization E->E[0].size(), wavelength calibEfield->E[0][0]->n_slices - 1
        // ... ... ...
        // Re(dE/dx_i) for source E->E.size(), polarization 0, wavelength calibEfield->E[0][0]->n_slices - 1
        // Im(dE/dx_i) for source E->E.size(), polarization 0, wavelength calibEfield->E[0][0]->n_slices - 1
        // ... ...
        // Re(dE/dx_i) source E->E.size(), polarization E->E[0].size(), wavelength 0
        // Im(dE/dx_i) source E->E.size(), polarization E->E[0].size(), wavelength 0
        // ...
        // Re(dE/dx_i) source E->E.size(), polarization E->E[0].size(), wavelength calibEfield->E[0][0]->n_slices - 1
        // Im(dE/dx_i) source E->E.size(), polarization E->E[0].size(), wavelength calibEfield->E[0][0]->n_slices - 1
        for (int s=0; s<jacEfield->E.size(); s++) {
            for (int p=0; p<jacEfield->E[s].size(); p++) {
                for (int k=0; k<jacEfield->E[s][p]->n_slices; k++) {
                    int n = regionPixelIndex.n_elem;
                    int nWavelengths = jacEfield->E[s][p]->n_slices;
                    int nPolarizations = jacEfield->E[s].size();
                    int firstRow = 2*k*n + p*2*n*nWavelengths + s*2*n*nWavelengths*nPolarizations;
                    int lastRow = firstRow + regionPixelIndex.n_elem - 1;
                    tmpE = (jacEfield->E[s][p]->slice(k) - fullEfield->E[s][p]->slice(k))/jacobianDx;
                    //                        sprintf(outName, "tmpE_%d.fits", i);
                    //                        save_mat(outName, tmpE);
                    //                        std::cout << std::endl;
                    //                        std::cout << "size(objMat) = " << size(objMat) << ", size(regionPixelIndex) = " << size(regionPixelIndex) << std::endl;
                    //                        std::cout << "k=" << k << ", i=" << i << ": firstRow = " << firstRow << ", lastRow = " << lastRow << std::endl;
                    //                        std::cout << "i=" << i << ": firstRow + optVec.n_elem = " << firstRow + regionPixelIndex.n_elem << ", lastRow + optVec.n_elem = " << lastRow + regionPixelIndex.n_elem << std::endl;
                    //                        fflush(stdout);
                    objMat.submat(firstRow, i, lastRow, i) = arma::real(tmpE(regionPixelIndex));
                    //                        std::cout << "set real part" << std::endl;
                    objMat.submat(firstRow + regionPixelIndex.n_elem, i, lastRow + regionPixelIndex.n_elem, i) = arma::imag(tmpE(regionPixelIndex));
                }
            }
        }
    }
    std::cout << std::endl;
}

void linearOptimizer::compute_objective_vector(void)
{

    arma::cx_mat tmpE(fullEfield->E[0][0]->n_rows, fullEfield->E[0][0]->n_cols);
    // set up objective vector
    for (int s=0; s<fullEfield->E.size(); s++) {
        for (int p=0; p<fullEfield->E[s].size(); p++) {
            for (int k=0; k<fullEfield->E[s][p]->n_slices; k++) {
                int n = regionPixelIndex.n_elem;
                int nWavelengths = jacEfield->E[s][p]->n_slices;
                int nPolarizations = jacEfield->E[s].size();
                int firstRow = 2*k*n + p*2*n*nWavelengths + s*2*n*nWavelengths*nPolarizations;
                int lastRow = firstRow + regionPixelIndex.n_elem - 1;
                tmpE = fullEfield->E[s][p]->slice(k);
                //                std::cout << "size(tmpE(regionPixelIndex)) = " << size(tmpE(regionPixelIndex)) << std::endl;
                //                fflush(stdout);
                objVec.subvec(firstRow, lastRow) = -arma::real(tmpE(regionPixelIndex));
                objVec.subvec(firstRow + regionPixelIndex.n_elem, lastRow + regionPixelIndex.n_elem) = -arma::imag(tmpE(regionPixelIndex));
            }
        }
    }
}

void linearOptimizer::print(const char *hdr) {
    if (hdr != NULL)
        std::cout << hdr << ":" << std::endl;
//    std::cout << "initX: " << initX << std::endl;
    std::cout << "deltaOptVecTolerance: " << deltaOptVecTolerance << std::endl;
//    std::cout << "global bounds: " << arma::min(globalLowerBound) << " to " << arma::max(globalUpperBound) << std::endl;
//    std::cout << "init bounds: " << arma::min(initLowerBound) << " to " << arma::max(initUpperBound) << std::endl;
//    std::cout << "maxIterations: " << maxIterations << std::endl;
//    std::cout << "stopContrastValue: " << stopContrastValue << std::endl;
//    std::cout << "outputInterval: " << outputInterval << std::endl;
//    std::cout << "outputDirectory: " << outputDirectory << std::endl;
//    std::cout << "saveState: " << saveState << std::endl;
    fflush(stdout);
}
