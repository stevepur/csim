//
//  zoomFftPropagator.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "zoomFftPropagator.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include "../coronagraph/coronagraph.hpp"

zoomFftPropagator::zoomFftPropagator() {
}

zoomFftPropagator::zoomFftPropagator(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

efield* zoomFftPropagator::execute(efield* E, celem* prev, celem* next, double time) {
    double z;
    
    
    std::cout << "calibrating: " << globalCoronagraph->get_calibration_state() << ", disableForCalibration: " << disableForCalibration << std::endl;
    if (!(disableForCalibration && globalCoronagraph->get_calibration_state())) {
        int nRowsE = E->E[0][0]->n_rows;
        int nColsE = E->E[0][0]->n_cols;
        arrayGeom agOut;
        
        pre_execute(E, prev, next, time);

        std::cout << "executed zoomFftPropagator " << name << std::endl;
        print("execute parameters:");
        E->arrayGeometry.print("zoomFftPropagator input E geometry: ");
        double lambdaFocalLength;
        if (dir == 1) {
            double dxprime = 2*(E->beamRadiusPhysical/(E->pixelScale*nRowsE*zoomFactor))*E->lambdaData[0].lambda*80;
            std::cout << "dxprime = " << dxprime << std::endl;
            
            agOut.pixelSizeX = dxprime;
            agOut.pixelSizeY = dxprime;
            int nRows = padFactor*nRowsE;
            int nCols = padFactor*nColsE;
            
            agOut.pixelX = arma::zeros<arma::vec>(nCols);
            for (int i=0; i<nCols; i++)
                agOut.pixelX[i] = ((double) i - ((double)nCols)/2.)*agOut.pixelSizeX;
            agOut.pixelY = arma::zeros<arma::vec>(nRows);
            for (int i=0; i<nRows; i++)
                agOut.pixelY[i] = ((double) i - ((double)nRows)/2.)*agOut.pixelSizeY;
            agOut.physicalSize = agOut.pixelX.n_elem*agOut.pixelSizeX;
            agOut.physicalSizeX = agOut.pixelX.n_elem*agOut.pixelSizeX;
            agOut.physicalSizeY = agOut.pixelY.n_elem*agOut.pixelSizeY;
            lambdaFocalLength = E->lambdaData[0].lambda*3.52;
        } else if (dir == -1) {
            agOut.pixelSizeX = E->pixelScale;
            agOut.pixelSizeY = E->pixelScale;
            int nRows = nRowsE/padFactor;
            int nCols = nColsE/padFactor;
            
            agOut.pixelX = arma::zeros<arma::vec>(nCols);
            for (int i=0; i<nCols; i++)
                agOut.pixelX[i] = ((double) i - ((double)nCols)/2.)*agOut.pixelSizeX;
            agOut.pixelY = arma::zeros<arma::vec>(nRows);
            for (int i=0; i<nRows; i++)
                agOut.pixelY[i] = ((double) i - ((double)nRows)/2.)*agOut.pixelSizeY;
            agOut.physicalSize = agOut.pixelX.n_elem*agOut.pixelSizeX;
            agOut.physicalSizeX = agOut.pixelX.n_elem*agOut.pixelSizeX;
            agOut.physicalSizeY = agOut.pixelY.n_elem*agOut.pixelSizeY;
            lambdaFocalLength = -E->lambdaData[0].lambda*3.52;
        }
        agOut.print("agOut: ");
        propZoomFft.init(E->arrayGeometry, agOut, &lambdaFocalLength, E->E[0][0]->n_slices);
//        propZoomFft.init(E->arrayGeometry.pixelX.n_elem, E->arrayGeometry.pixelY.n_elem, zoomFactor, dir, E->E[0][0]->n_slices);
        
        for (int s=0; s<E->E.size(); s++) {
            for (int p=0; p<E->E[s].size(); p++) {
                *(E->E[s][p]) = propZoomFft.execute(*(E->E[s][p]));
            }
        }
        E->arrayGeometry = agOut;
        E->print("post zoomFftPropagator E: ");
        
        post_execute(E, prev, next, time);
    }
    
    return E;
}

void zoomFftPropagator::init(initCommandSet*& cmdBlock) {
    std::cout << "initing a zoomFftPropagator" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    post_init();
}

void zoomFftPropagator::set(std::string fieldName, const char *arg) {
    bool found = celem::set(fieldName, arg);
    if (fieldName == "zoomFftPropagator")
    ;
    else if (fieldName == "dir") {
        dir = atoi(arg);
    }
    else if (fieldName == "zoomFactor") {
        zoomFactor = atof(arg);
    }
    else if (fieldName == "padFactor") {
        padFactor = atof(arg);
    }
    else if (!found)
        std::cout << "!!! zoomFftPropagator bad set field name: " << fieldName << std::endl;
}

void zoomFftPropagator::print(const char *hdr) {
    std::cout << "zoomFftPropagator " << hdr << std::endl;
    std::cout << "dir = " << dir << std::endl;
    std::cout << "zoomFactor = " << zoomFactor << std::endl;
}

