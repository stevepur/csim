//
//  downsample.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "downsample.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include "../lib/csim_lib.hpp"

downsample::downsample() {
}

downsample::downsample(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

efield* downsample::execute(efield* E, celem* prev, celem* next, double time) {
//    std::cout << "executed a downsample " << name << std::endl;
    std::complex<double> i1(0, 1);

    pre_execute(E, prev, next, time);
    
//    draw_mat(abs(E->E[0][0]->slice(0)), "E before downsampling");

    int nRows = E->E[0][0]->n_rows;
    int nCols = E->E[0][0]->n_cols;
    arma::mat kernel = arma::ones<arma::mat>(sampleFactor, sampleFactor);
    
    arrayGeom newGeom;
    int newNRows = (int) ceil(nRows/sampleFactor);
    int newNCols = (int) ceil(nCols/sampleFactor);
    double newPixelScaleX = sampleFactor*E->arrayGeometry.pixelSizeX;
    double newPixelScaleY = sampleFactor*E->arrayGeometry.pixelSizeY;
    newGeom.set_geometry(newNRows, newNCols, newPixelScaleX, newPixelScaleY, E->arrayGeometry.originX, E->arrayGeometry.originY);
//    E->arrayGeometry.print("original geometry");
//    newGeom.print("new geometry");
    

    for (int s=0; s<E->E.size(); s++) {
        for (int p=0; p<E->E[s].size(); p++) {
            arma::cx_cube *outE = new arma::cx_cube;
            outE->zeros(newNRows, newNCols, E->E[0][0]->n_slices);

            for (int k=0; k<E->E[s][p]->n_slices; k++) {
#if 0
                arma::cx_mat paddedE = arma::zeros<arma::cx_mat>(nRows + 2*sampleFactor, nCols + 2*sampleFactor);
//                std::cout << "size paddedE: " << size(paddedE) << std::endl;
                arma::cx_mat Eslice = E->E[s][p]->slice(k);
//                std::cout << "filling paddedE" << std::endl;
//                std::cout << "size submat paddedE: " << size(paddedE.submat(sampleFactor, sampleFactor, paddedE.n_rows-sampleFactor-1, paddedE.n_cols-sampleFactor-1)) << std::endl;
                paddedE.submat(sampleFactor, sampleFactor, paddedE.n_rows-sampleFactor-1, paddedE.n_cols-sampleFactor-1) = Eslice;
//                std::cout << "setting low rows" << std::endl;
                paddedE.submat(0, sampleFactor, sampleFactor-1, paddedE.n_cols-sampleFactor-1) = repmat(Eslice.row(0), sampleFactor, 1);
//                std::cout << "setting high rows" << std::endl;
                paddedE.submat(paddedE.n_rows-sampleFactor, sampleFactor, paddedE.n_rows-1, paddedE.n_cols-sampleFactor-1) = repmat(Eslice.row(nRows-1), sampleFactor, 1);
//                std::cout << "setting low cols" << std::endl;
                paddedE.submat(sampleFactor, 0, paddedE.n_rows-sampleFactor-1, sampleFactor-1) = repmat(Eslice.col(0), 1, sampleFactor);
//                std::cout << "setting high cols" << std::endl;
                paddedE.submat(sampleFactor, paddedE.n_cols-sampleFactor, paddedE.n_rows-sampleFactor-1, paddedE.n_cols-1) = repmat(Eslice.col(nCols-1), 1, sampleFactor);
                
                arma::mat tmpC = arma::conv2(real(paddedE), kernel, "same");
                arma::mat tmpE = tmpC.submat(sampleFactor, sampleFactor, paddedE.n_rows-sampleFactor-1, paddedE.n_cols-sampleFactor-1)/sum(sum(kernel));
#endif
                arma::mat tmpMat(newNRows, newNCols);
                arma::mat tmpEMat;
                arma::cx_mat intE(size(tmpMat));
                
//                tmpEMat = real(tmpE);
                tmpEMat = real(E->E[s][p]->slice(k));
                interp2(tmpEMat, tmpMat, E->arrayGeometry.pixelX, E->arrayGeometry.pixelY, newGeom.pixelX, newGeom.pixelY);
                intE.set_real(tmpMat);
//                tmpEMat = imag(tmpE);
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

void downsample::init(initCommandSet*& cmdBlock) {
//    std::cout << "initing a downsample" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    
    post_init();
}

void downsample::set(std::string fieldName, const char *arg) {

    bool found = celem::set(fieldName, arg);
    if (fieldName == "downsample")
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
        std::cout << "!!! downsample bad set field name: " << fieldName << std::endl;
}

