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

    nloptOptimizer* thisOptimizer = reinterpret_cast<nloptOptimizer*>(parentPointer);
    
    arma::vec optVec;
    std_vec_to_arma_vec(x, optVec);
//    std::cout << "optVec: " << optVec << std::endl;
    
    globalCoronagraph->set_optimization_data(thisOptimizer->componentName, thisOptimizer->dataName, &optVec);

//    std::cout << "calibMaxIntensity = " << thisOptimizer->calibMaxIntensity << std::endl;
    timer.tic();
    arma::mat fullIntensity;
    thisOptimizer->reset_fullEfield();
    double fullMaxIntensity = thisOptimizer->region->compute_intensity(thisOptimizer->fullEfield, fullIntensity, false)/thisOptimizer->calibMaxIntensity;
    fullIntensity /= thisOptimizer->calibMaxIntensity;
    
//    draw_mat(log10(fullIntensity), "fullIntensity", "matlab");
    double contrast = arma::mean(fullIntensity(thisOptimizer->regionPixelIndex));
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
        const char *cPtr = strchr(arg, ',');
        int componentNameLen = cPtr - arg;
        int dataNameLen = strlen(arg) - componentNameLen + 1;
        std::cout << "componentNameLen " << componentNameLen << ", dataNameLen " << dataNameLen << std::endl;
        componentName = new char[componentNameLen + 1];
        dataName = new char[dataNameLen + 1];
        strncpy(componentName, arg, componentNameLen);
        componentName[componentNameLen] = '\0';
        strncpy(dataName, cPtr + 1, dataNameLen);
        dataName[dataNameLen] = '\0';
        std::cout << "optimize component " << componentName << ", data " << dataName << std::endl;
        //        draw("initial FPM");
        
        globalCoronagraph->get_optimization_data(componentName, dataName, &startVec);
        optVecSize = startVec.n_elem;
        initX = startVec;
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
        float v1;
        // arg is a comma separated pair of doubles
        sscanf(arg, "%f", &v1);
        initX = v1*arma::ones<arma::vec>(optVecSize);
    }
    else if (fieldName == "xRelTolerance") {
        // arg is a single float
        xRelTolerance = atof(arg);
    }
    else if (fieldName == "fRelTolerance") {
        // arg is a single float
        fRelTolerance = atof(arg);
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
    delete fullEfield;
    fullEfield = new efield(*initialEfield);
    fullEfield->set("name", "FPM Efield");
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
    opt.set_upper_bounds(globalUpperBound);
    opt.set_lower_bounds(globalLowerBound);
    double minf;
    std::vector<double> x;
    arma_vec_to_std_vec(initX, x);
    std::cout << "x: " << x[0] << std::endl;
    
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
    globalCoronagraph->get_optimization_data(componentName, dataName, &finalVec);
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
