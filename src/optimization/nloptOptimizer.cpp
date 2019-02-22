//
//  nloptOptimizer.cpp
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
#include "nloptOptimizer.hpp"
#include "../coronagraph/coronagraph.hpp"

double nloptOptimizer::eval_contrast(const std::vector<double> &x, std::vector<double> &grad, void *parentPointer)
{
    arma::wall_clock timer;
    bool testMode = false;

    nloptOptimizer* thisOptimizer = reinterpret_cast<nloptOptimizer*>(parentPointer);
    
    //
    //
    // run the coronagraph to get the contrast for the current x
    //
    
    // send the x data to the coronagraph as an arma::vec
    arma::vec optVec;
    std_vec_to_arma_vec(x, optVec);
    std::cout << "optVec: " << optVec << std::endl;

    // compute the region intensity via the regionContrast object, which calls the coronagraph
//    std::cout << "calibMaxIntensity = " << thisOptimizer->calibMaxIntensity << std::endl;
    timer.tic();
    // reinitialize the E field
    
    arma::wall_clock setOpttimer;
    setOpttimer.tic();
//    globalCoronagraph->set_optimization_data(thisOptimizer->componentName, thisOptimizer->dataName, optVec);
    globalCoronagraph->set_optimization_data(optVec);
    std::cout << "time to set optimization data: " << setOpttimer.toc() << std::endl;
    double contrast;
    arma::vec tx(optVec.n_elem);
    arma::vec ty(optVec.n_elem);
    arma::vec testGoal(optVec.n_elem);
    if (testMode) {
        for (double i=0; i<optVec.n_elem; i++) {
            double gridSize = sqrt(optVec.n_elem);
            arma::uvec xy = ind2sub(arma::size(gridSize, gridSize), i);
            tx(i) = ((double) xy(0))/(32 - 1);
            ty(i) = ((double) xy(1))/(32 - 1);
        }
        testGoal = (0.5e-7*(tx%tx + ty%ty) + 1e-8);
        contrast = norm(optVec - testGoal);
    } else {
        arma::mat fullIntensity;
        arma::mat calibIntensity;
        thisOptimizer->reset_calibEfield();
        double calibMaxIntensity = thisOptimizer->region->compute_intensity(thisOptimizer->calibEfield, calibIntensity, true);
        thisOptimizer->reset_fullEfield();
        double fullMaxIntensity = thisOptimizer->region->compute_intensity(thisOptimizer->fullEfield, fullIntensity, false)/calibMaxIntensity;
        // normalize
        fullIntensity /= calibMaxIntensity;
        
    //    draw_mat(log10(fullIntensity), "fullIntensity", "matlab");
        // compute the contrast
        contrast = arma::mean(fullIntensity(thisOptimizer->regionPixelIndex));
    }
    
    // output data from this step
    thisOptimizer->increment_nIterations();
    std::cout << "evaluation " << thisOptimizer->get_nIterations() << ", " << timer.toc() << " seconds, contrast = " << contrast << std::endl;
    
    fprintf(thisOptimizer->get_histValFid(), "%g %g\n", thisOptimizer->optimizationTimer.toc(), contrast);
    fflush(thisOptimizer->get_histValFid());
    if (contrast < thisOptimizer->get_bestVal()) {
        fprintf(thisOptimizer->get_bestValFid(), "%g %g\n", thisOptimizer->optimizationTimer.toc(), contrast);
        fflush(thisOptimizer->get_bestValFid());
        thisOptimizer->set_bestVal(contrast);
        globalCoronagraph->save_optimization_data(thisOptimizer->componentName, thisOptimizer->dataName, thisOptimizer->get_ouput_directory());
    }

    //
    //
    // if the optimization method needs a gradient, compute it.
    //
    if (!grad.empty()) {
        // compute the error beween the current contrast and the prediction using the previous gradient
        double gradientError = 1e6;
        if (thisOptimizer->previousGradientX.n_elem > 0) {
            arma::vec deltaX = thisOptimizer->previousGradientX - optVec;
            double contrastEstimate = thisOptimizer->previousGradientContrast + 0.1*dot(thisOptimizer->previousGradient, deltaX);
            gradientError = (contrastEstimate - contrast)/contrast;
            std::cout << "relative gradient error = " << gradientError << std::endl;
        }
        
        // if the error is too big or if the contrast is no longer decreasing, compute the gradient
        std::cout << "contrast = " << contrast << ", lastContrast = " << thisOptimizer->lastContrast << std::endl;
        std::cout << "fabs(gradientError) = " << fabs(gradientError) << ", gradientErrorThreshold = " << thisOptimizer->gradientErrorThreshold << ", norm(previousGradient) = " << norm(thisOptimizer->previousGradient) << std::endl;
        if (fabs(gradientError) > thisOptimizer->gradientErrorThreshold | contrast > thisOptimizer->lastContrast) {
            std::cout << "computing gradient, gradientDx = " << thisOptimizer->gradientDx << ": ";
            arma::vec dx;
            arma::vec gradient(optVec.n_elem);
            arma::mat gradientIntensity;
            for (int i=0; i<optVec.n_elem; i++) {
                dx = optVec;
                dx(i) = dx(i) + thisOptimizer->gradientDx;
                double c;
                if (testMode) {
                    c = norm(dx - testGoal);
                } else {
                    if (!(i%10) | optVec.n_elem <= 20) {
                        std::cout << i << " ";
                        fflush(stdout);
                    }
//                    globalCoronagraph->set_optimization_data(thisOptimizer->componentName, thisOptimizer->dataName, dx);
                    globalCoronagraph->set_optimization_data(dx);
                    thisOptimizer->reset_fullEfield();
                    double gradientMaxIntensity = thisOptimizer->region->compute_intensity(thisOptimizer->fullEfield, gradientIntensity, false)/thisOptimizer->calibMaxIntensity;
                    gradientIntensity /= thisOptimizer->calibMaxIntensity;
                    c = arma::mean(gradientIntensity(thisOptimizer->regionPixelIndex));
                }
                gradient(i) = -(contrast - c)/thisOptimizer->gradientDx;
            }
            std::cout << std::endl;
            
//            std::cout << "gradient: " << gradient << std::endl;
            std::cout << "norm of gradient = " << norm(gradient) << std::endl;

            arma_vec_to_std_vec(gradient, grad);
            thisOptimizer->previousGradient = gradient;
            thisOptimizer->previousGradientX = optVec;
            thisOptimizer->previousGradientContrast = contrast;
        } else {
            arma::vec tg = 0.1*thisOptimizer->previousGradient;
            arma_vec_to_std_vec(tg, grad);
        }
    }
    thisOptimizer->lastContrast = contrast;

    return contrast;
}

nloptOptimizer::nloptOptimizer() {
}

nloptOptimizer::nloptOptimizer(initCommandSet*& cmdBlocks) {
    init(cmdBlocks);
}

void nloptOptimizer::init(initCommandSet*& cmdBlock) {
    
    std::vector<initCommandSet*> subBlocks = cmdBlock->find_command_blocks();
    
    for (int j=0; j<subBlocks.size(); j++) {
        if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "nloptOptimizer")) {
            init_block(subBlocks[j]);
        }
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "regionContrast")) {
            region = new regionContrast(subBlocks[j]);
            region->print("region");
        }
    }
    
    globalCoronagraph->get_optimization_data(startVec);
    optVecSize = startVec.n_elem;
    initX = startVec;
    
    print("initialized nlopt optimizer");
}

void nloptOptimizer::init_block(initCommandSet*& cmdBlock) {
    std::vector<initCommandSet*> subBlocks = cmdBlock->find_command_blocks();
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
}

void nloptOptimizer::set(std::string fieldName, const char *arg) {
//    std::cout << "fieldname: " << fieldName << ", arg: " << arg << std::endl;
    if (fieldName == "nloptOptimizer")
        ;
    else if (fieldName == "dataToOptimize") {
        // arg is two strings separated by a comma,
        // the first string is the name of the component to be optimized
        // the second string is the name of the data to be optimized
/*
 const char *cPtr = strchr(arg, ',');
        int componentNameLen = cPtr - arg;
        int dataNameLen = strlen(arg) - componentNameLen + 1;
        componentName = new char[componentNameLen + 1];
        dataName = new char[dataNameLen + 1];
        strncpy(componentName, arg, componentNameLen);
        componentName[componentNameLen] = '\0';
        strncpy(dataName, cPtr + 1, dataNameLen);
        dataName[dataNameLen] = '\0';
 */
        componentName = new char[100];
        dataName = new char[100];
        float lb;
        float ub;
        sscanf(arg, "%s %s %f %f", componentName, dataName, &lb, &ub);
        std::cout << "optimize component " << componentName << ", data " << dataName << ", lower bound " << lb << ", upper bound " << ub << std::endl;
        //        draw("initial FPM");
        
        globalCoronagraph->add_optimization_data(componentName, dataName, lb, ub);
    }
    else if (fieldName == "optMethod") {
        // arg is a string
        if (!strcmp(arg, "nelderMead")) {
            optMethod = nlopt::LN_NELDERMEAD;
        }
        else if (!strcmp(arg, "subPlex")) {
            optMethod = nlopt::LN_SBPLX;
        }
        else if (!strcmp(arg, "cobyla")) {
            optMethod = nlopt::LN_COBYLA;
        }
        else if (!strcmp(arg, "bobyqa")) {
            optMethod = nlopt::LN_BOBYQA;
        }
        else if (!strcmp(arg, "newuoa")) {
            optMethod = nlopt::LN_NEWUOA_BOUND;
        }
        else if (!strcmp(arg, "mma")) {
            optMethod = nlopt::LD_MMA;
        }
        else if (!strcmp(arg, "bfgs")) {
            optMethod = nlopt::LD_LBFGS;
        }
        else if (!strcmp(arg, "tNewtonPrecondRestart")) {
            optMethod = nlopt::LD_TNEWTON_PRECOND_RESTART;
        }
        else if (!strcmp(arg, "tNewtonPrecond")) {
            optMethod = nlopt::LD_TNEWTON_PRECOND;
        }
        else if (!strcmp(arg, "tNewtonRestart")) {
            optMethod = nlopt::LD_TNEWTON_RESTART;
        }
        else if (!strcmp(arg, "tNewton")) {
            optMethod = nlopt::LD_TNEWTON;
        }
        else {
            std::cout << "bad optimization method: " << arg << std::endl;
            assert(NULL);
        }
    }
    else if (fieldName == "initX") {
//        assert(optVecSize == 2);
//        initX.set_size(optVecSize);
//        float v1, v2;
//        // arg is a comma separated pair of doubles
//        sscanf(arg, "%f, %f", &v1, &v2);
//        initX[0] = v1;
//        initX[1] = v2;
//        assert(optVecSize == 1);
        initX.set_size(optVecSize);
        if (!strcmp(arg, "random")) {
            initX = globalUpperBound*(arma::randu(optVecSize) - 0.5);
            std::cout << "mean(initX) = " << mean(initX) << " min(initX) = " << min(initX) << " max(initX) = " << max(initX) << std::endl;
        } else {
            float v1;
            // arg is a comma separated pair of doubles
            sscanf(arg, "%f", &v1);
            initX = v1*arma::ones<arma::vec>(optVecSize);
        }
    }
    else if (fieldName == "xRelTolerance") {
        // arg is a single float
        xRelTolerance = atof(arg);
    }
    else if (fieldName == "fRelTolerance") {
        // arg is a single float
        fRelTolerance = atof(arg);
    }
    else if (fieldName == "gradientDx") {
        // arg is a single float
        gradientDx = atof(arg);
    }
    else if (fieldName == "gradientErrorThreshold") {
        // arg is a single float
        gradientErrorThreshold = atof(arg);
    }
    else if (fieldName == "globalBounds") {
        float v1, v2;
        // arg is a comma separated pair of doubles
        sscanf(arg, "%f, %f", &v1, &v2);
        globalLowerBound = v1;
        globalUpperBound = v2;
    }
    else if (fieldName == "initBounds") {
        float v1, v2;
        // arg is a comma separated pair of doubles
        sscanf(arg, "%f, %f", &v1, &v2);
        initLowerBound = v1;
        initUpperBound = v2;
    }
    // set stopping criteria
    else if (fieldName == "maxIterations") {
        // arg is a single int
        maxIterations = atoi(arg);
    }
    else if (fieldName == "optVecSize") {
        // arg is a single int
        optVecSize = atoi(arg);
        std::cout << "optVecSize: " << optVecSize << std::endl;
    }
    else if (fieldName == "stopObjectiveValue") {
        // arg is a single float
        stopObjectiveValue = atof(arg);
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
        std::cout << "!!! nloptOptimizer bad set field name: " << fieldName << std::endl;
    
}

void nloptOptimizer::reset_fullEfield(void){
    if (fullEfield != NULL)
        delete fullEfield;
    fullEfield = new efield(*initialEfield);
    fullEfield->set("name", "FPM Efield");
}

void nloptOptimizer::reset_calibEfield(void){
    if (calibEfield != NULL)
        delete calibEfield;
    calibEfield = new efield(*initialEfield);
    calibEfield->set("name", "calibration Efield");
}

void nloptOptimizer::optimize(void) {
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


    nlopt::opt opt(optMethod, optVecSize);
    opt.set_min_objective(nloptOptimizer::eval_contrast, this);
    opt.set_stopval(stopObjectiveValue);
    opt.set_maxeval(maxIterations);
    opt.set_ftol_rel(fRelTolerance);
    opt.set_xtol_rel(xRelTolerance);
    
    globalCoronagraph->get_optimization_bounds(lowerBounds, upperBounds);
    std::vector<double> lbv;
    std::vector<double> ubv;
    arma_vec_to_std_vec(lowerBounds, lbv);
    arma_vec_to_std_vec(upperBounds, ubv);
    opt.set_lower_bounds(lbv);
    opt.set_upper_bounds(ubv);
    double minf;
    std::vector<double> x;
    arma_vec_to_std_vec(initX, x);
    std::cout << "x: " << x[0] << std::endl;
    std::vector<double> stepSizeV;
    stepSizeV.resize(optVecSize, 0.0);
    opt.get_initial_step(x, stepSizeV);
    arma::vec stepSize;
    std_vec_to_arma_vec(stepSizeV, stepSize);
//    opt.set_initial_step(1e-8);
    std::cout << "mean step size = " << mean(stepSize) << std::endl;
    std::cout << "max step size = " << max(stepSize) << std::endl;
    std::cout << "min step size = " << min(stepSize) << std::endl;

    std::cout << "============= starting optimization =============" << std::endl;
    std::string fname = (std::string)outputDirectory + "/optBestVal.txt";
    bestValFid = fopen(fname.c_str(), "w");
    
    fname = (std::string)outputDirectory + "/optValHistory.txt";
    histValFid = fopen(fname.c_str(), "w");
    
    timer.tic();
    optimizationTimer.tic();
    try{
        nlopt::result result = opt.optimize(x, minf);
        std::cout << "found minimum contrast " << std::setprecision(10) << minf << std::endl;
        std::cout << "result " << result << std::endl;
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    fclose(bestValFid);
    fclose(histValFid);
    std::cout << "nIterations: " << nIterations << ", time " << timer.toc() << " seconds" << std::endl;
//    std::cout << "original opt values: " << startVec << std::endl;
    globalCoronagraph->get_optimization_data(componentName, dataName, finalVec);
//    std::cout << "final opt values: " << finalVec << std::endl;
    std::cout << "max difference between original and final: " << max(abs(finalVec - startVec)) << std::endl;
    globalCoronagraph->save_optimization_data(componentName, dataName);
}

void nloptOptimizer::print(const char *hdr) {
    if (hdr != NULL)
        std::cout << hdr << ":" << std::endl;
//    std::cout << "initX: " << initX << std::endl;
    std::cout << "xRelTolerance: " << xRelTolerance << std::endl;
//    std::cout << "global bounds: " << arma::min(globalLowerBound) << " to " << arma::max(globalUpperBound) << std::endl;
//    std::cout << "init bounds: " << arma::min(initLowerBound) << " to " << arma::max(initUpperBound) << std::endl;
//    std::cout << "maxIterations: " << maxIterations << std::endl;
//    std::cout << "stopObjectiveValue: " << stopObjectiveValue << std::endl;
//    std::cout << "outputInterval: " << outputInterval << std::endl;
//    std::cout << "outputDirectory: " << outputDirectory << std::endl;
//    std::cout << "saveState: " << saveState << std::endl;
    fflush(stdout);
}
