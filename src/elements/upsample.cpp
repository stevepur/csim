//
//  upsample.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "upsample.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include "../lib/csim_lib.hpp"

upsample::upsample() {
}

upsample::upsample(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

efield* upsample::execute(efield* E, celem* prev, celem* next, double time) {
//    std::cout << "executed a upsample " << name << std::endl;
    std::complex<double> i1(0, 1);

    pre_execute(E, prev, next, time);
    
//    draw_mat(abs(E->E[0][0]->slice(0)), "E before downsampling");

    int nRows = E->E[0][0]->n_rows;
    int nCols = E->E[0][0]->n_cols;
    
    arrayGeom newGeom;
    int newNRows = (int) ceil(nRows*sampleFactor);
    int newNCols = (int) ceil(nCols*sampleFactor);
    double newPixelScaleX = E->arrayGeometry.pixelSizeX/sampleFactor;
    double newPixelScaleY = E->arrayGeometry.pixelSizeY/sampleFactor;
    newGeom.set_geometry(newNRows, newNCols, newPixelScaleX, newPixelScaleY, E->arrayGeometry.originX, E->arrayGeometry.originY);
//    E->arrayGeometry.print("original geometry");
//    newGeom.print("new geometry");
    

    for (int s=0; s<E->E.size(); s++) {
        for (int p=0; p<E->E[s].size(); p++) {
            arma::cx_cube *outE = new arma::cx_cube;
            outE->zeros(newNRows, newNCols, E->E[0][0]->n_slices);

            for (int k=0; k<E->E[s][p]->n_slices; k++) {
                arma::mat tmpMat(newNRows, newNCols);
                arma::mat tmpEMat;
                arma::cx_mat intE(size(tmpMat));
                
                tmpEMat = real(E->E[s][p]->slice(k));
                interp2(tmpEMat, tmpMat, E->arrayGeometry.pixelX, E->arrayGeometry.pixelY, newGeom.pixelX, newGeom.pixelY);
                intE.set_real(tmpMat);
                tmpEMat = imag(E->E[s][p]->slice(k));
                interp2(tmpEMat, tmpMat, E->arrayGeometry.pixelX, E->arrayGeometry.pixelY, newGeom.pixelX, newGeom.pixelY);
                intE.set_imag(tmpMat);
                
                outE->slice(k) = intE;
            }
            E->E[s].at(p) = outE;

        }
    }
    E->arrayGeometry = newGeom;
    
//    draw_mat(abs(E->E[0][0]->slice(0)), "E after downsampling");

    post_execute(E, prev, next, time);
    
    return E;
}

void upsample::init(initCommandSet*& cmdBlock) {
//    std::cout << "initing a upsample" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    
    post_init();
}

void upsample::set(std::string fieldName, const char *arg) {

    bool found = celem::set(fieldName, arg);
    if (fieldName == "upsample")
        ;
    else if (fieldName == "sampleFactor") {
        // arg is a single int
        sampleFactor = atoi(arg);
        std::cout << "sampleFactor = " << sampleFactor << std::endl;
        if (sampleFactor < 1.0) {
            std::cout << "!!! sampleFactor < 1 not supported!: " << sampleFactor << std::endl;
            assert(NULL);
        }
    }
    else if (!found)
        std::cout << "!!! upsample bad set field name: " << fieldName << std::endl;
}

