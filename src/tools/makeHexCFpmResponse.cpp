//
//  makeHexCFpmResponse.cpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>

#include "armadillo"
#include "makeHexCFpmResponse.hpp"
#include "../coronagraph/coronagraph.hpp"
#include "../telescope/telescope.hpp"
#include "../lib/csim_lib.hpp"
#include "../data/responseData.hpp"

makeHexCFpmResponse::makeHexCFpmResponse() {
}

makeHexCFpmResponse::makeHexCFpmResponse(initCommandSet*& cmdBlocks) {
    init(cmdBlocks);
}

void makeHexCFpmResponse::init(initCommandSet*& cmdBlock) {
    
    std::vector<initCommandSet*> subBlocks = cmdBlock->find_command_blocks();
    
    for (int j=0; j<subBlocks.size(); j++) {
        if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "makeHexCFpmResponse")) {
            init_block(subBlocks[j]);
        }
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "regionContrast")) {
            region = new regionContrast(subBlocks[j]);
            region->print("region");
        }
    }
}

void makeHexCFpmResponse::init_block(initCommandSet*& cmdBlock) {
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
}

void makeHexCFpmResponse::set(std::string fieldName, const char *arg) {
    if (fieldName == "makeHexCFpmResponse")
        ;
    else if (fieldName == "DataToOptimize") {
        // arg is two strings separated by a comma,
        // the first string is the name of the component to be optimized
        // the second string is the name of the data to be optimized
        char *cPtr = strchr(arg, ',');
        int componentNameLen = cPtr - arg;
        int dataNameLen = strlen(arg) - componentNameLen + 1;
        std::cout << "componentNameLen " << componentNameLen << ", dataNameLen " << dataNameLen << std::endl;
        componentName = new char[componentNameLen + 1];
        dataName = new char[dataNameLen + 1];
        strncpy(componentName, arg, componentNameLen);
        componentName[componentNameLen] = '\0';
        strncpy(dataName, cPtr + 1, dataNameLen);
        dataName[dataNameLen] = '\0';
        std::cout << "component " << componentName << ", data " << dataName << std::endl;
        //        draw("initial FPM");
    }
    else if (fieldName == "draw") {
        // arg is a string
        if (!strcmp(arg, "on")) {
            draw = true;
        }
        else if (!strcmp(arg, "off")) {
            draw = false;
        }
        else
            std::cout << "!!! makeHexCFpmResponse draw bad set field name: " << fieldName << std::endl;
    }
    else
        std::cout << "!!! makeHexCFpmResponse bad set field name: " << fieldName << std::endl;
    
}

void makeHexCFpmResponse::compute_response(void) {
    arma::wall_clock timer;
    
    // compute the calibration normalization
    efield *calibEfield = new efield(*initialEfield);
    calibEfield->set("name", "calibration Efield");
    
    globalCoronagraph->set_calibration_state(true);
    timer.tic();
    globalCoronagraph->execute(calibEfield, 0);
    std::cout << "contrast curve calibration csim execution time: " << timer.toc() << " seconds" << std::endl;
    
    arma::cube calibIntensity;
    calibIntensity.zeros(size(*(calibEfield->E[0][0])));
    for (int s=0; s<calibEfield->E.size(); s++) {
        for (int p=0; p<calibEfield->E[s].size(); p++) {
            calibIntensity += real(*(calibEfield->E[s][p]) % arma::conj(*(calibEfield->E[s][p])));
        }
    }
    arma::mat calibIntensitySum = sum(calibIntensity, 2);
    double calibMaxIntensity = max(max(calibIntensitySum));
    std::cout << "calibMaxIntensity = " << calibMaxIntensity << std::endl;
    
    // get the data type to be computed (sags in this case)
    // in order to know how many sags we have
    arma::vec origHexSags;
    globalCoronagraph->get_optimization_data(componentName, dataName, &origHexSags);
    
    int nSags = origHexSags.n_elem;
//    nSags = 1;
    std::cout << "nSags = " << nSags << std::endl;
   
    // set all sags to zero
    arma::vec curSag = arma::zeros(size(origHexSags));
    globalCoronagraph->set_optimization_data(componentName, dataName, &curSag);
    
    region->print("response region:");
    arma::uvec regionPixelIndex;
    // be careful to use the result of an execution so we're using the output CCD geometry
    region->get_region_pixels(calibEfield, regionPixelIndex);
//    std::cout << "there are " << regionPixelIndex.n_elem << " pixels: " << std::endl;
//    std::cout << regionPixelIndex << std::endl;
    
//    responseData *response = new responseData(regionPixelIndex.n_elem, nSags, calibEfield->E[0][0]->n_slices, calibEfield->E[0].size(), calibEfield->E.size());
    responseData *response = new responseData(regionPixelIndex.n_elem, nSags+1, calibEfield->E[0][0]->n_slices, calibEfield->E[0].size(), calibEfield->E.size());
    response->set((std::string)"outputDirectory", "testResponse");
    response->print("response data");
    
    globalCoronagraph->set_calibration_state(false);
    for (int sag=0; sag<=nSags; sag++){
//    for (int sag=0; sag<nSags; sag++){
        globalCoronagraph->set_optimization_data(componentName, "useOnlyThisHex", &sag);
        
        bool doBabinet = false;
        if (sag < nSags)
            doBabinet = false;
        else if (sag == nSags)
            doBabinet = true;
            
        globalCoronagraph->set_optimization_data(componentName, "setBabinet", &doBabinet);
        
        // run the coronagraph
        efield *fullEfield = new efield(*initialEfield);
        fullEfield->set("name", "FPM Efield");

        timer.tic();
        globalCoronagraph->execute(fullEfield, 0);
        std::cout << "make hex response function csim execution time: " << timer.toc() << " seconds" << std::endl;
        
        for (int s=0; s<fullEfield->E.size(); s++) {
            for (int p=0; p<fullEfield->E[s].size(); p++) {
                for (int sl=0; sl<response->M[s][p]->n_slices; sl++) {
                    arma::cx_mat smat = fullEfield->E[s][p]->slice(sl);
                    arma::cx_mat tmat = response->M[s][p]->slice(sl);
                    tmat(arma::span::all, sag) = smat(regionPixelIndex);
                    response->M[s][p]->slice(sl) = tmat;
                }
            }
        }
        

        delete fullEfield;
       
    }
//    draw_mat(real(response->M[0][0]->slice(0)), "real responese");
//    draw_mat(imag(response->M[0][0]->slice(0)), "imag responese");
//    save_mat("response0.fits", response->M[0][0]->slice(0));
    response->save();
    
}

