//
//  doubleFraunhoferInterpFPM.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "doubleFraunhoferInterpFPM.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include <stdio.h>
#include "../coronagraph/coronagraph.hpp"

doubleFraunhoferInterpFPM::doubleFraunhoferInterpFPM() {
}

doubleFraunhoferInterpFPM::doubleFraunhoferInterpFPM(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

efield* doubleFraunhoferInterpFPM::execute(efield* E, celem* prev, celem* next, double time) {
    std::complex<double> i1(0, 1);
    
    
    std::cout << "calibrating: " << globalCoronagraph->get_calibration_state() << ", disableForCalibration: " << disableForCalibration << std::endl;
    if (!(disableForCalibration && globalCoronagraph->get_calibration_state())) {
        
        std::cout << "executing doubleFraunhoferInterpFPM " << name << std::endl;
        E->arrayGeometry.print("doubleFraunhoferInterpFPM input E geometry: ");
        E->print("at top of doubleFraunhoferInterpFPM execute:");
        
        pre_execute(E, prev, next, time);

        for (int s=0; s<E->E.size(); s++) {
            for (int p=0; p<E->E[s].size(); p++) {
                for (int sl=0; sl<E->E[s][p]->n_slices; sl++) {
                
                    // the actual execution
                    arma::cx_mat tcMat;

                    std::cout << "executing doubleFraunhoferInterpFPM slice " << sl << std::endl;
                    std::cout << "lambda =  " << E->lambdaData[sl].lambda << std::endl;

                    // Fraunhofer propagation 1
                    tcMat = propZoomFftIn.execute(E->E[s][p]->slice(sl), sl);
//                    save_mat("tcMat_fp1", tcMat, "reIm");
                    
                    // apply the mask
                    hexFPM->execute(E, tcMat, sl, time);

                    // Fraunhofer propagation 2
                    E->E[s][p]->slice(sl) = propZoomFftOut.execute(tcMat, sl);
                }
            }
        }
        E->print("post doubleFraunhoferInterpFPM E: ");
        
        post_execute(E, prev, next, time);
    }
    
    return E;
}

void doubleFraunhoferInterpFPM::init(initCommandSet*& cmdBlock) {
    std::cout << "initing a doubleFraunhoferInterpFPM" << std::endl;
    
    std::vector<initCommandSet*> subBlocks = cmdBlock->find_command_blocks();
    for (int i=0; i<subBlocks.size(); i++) {
        std::cout << "processing " << subBlocks[i]->commandList[0]->getCmdStr() << std::endl;
        // define the FPM
        if (!strcmp(subBlocks[i]->commandList[0]->getCmdStr(), "complexHexMaskFPM")) {
            hexFPM = new complexHexMaskFPM(subBlocks[i], zoomFactor);
        }
        else if (!strcmp(subBlocks[i]->commandList[0]->getCmdStr(), "doubleFraunhoferInterpFPM")) {
            for (int c=0; c<subBlocks[i]->commandList.size(); c++) {
                set(subBlocks[i]->commandList[c]->getCmdStr(),
                    subBlocks[i]->commandList[c]->getArgStr());
            }
        }
    }
    init_geom();
    post_init();
}

void doubleFraunhoferInterpFPM::init_geom(void) {
    int nRowsE = initialEfield->E[0][0]->n_rows;
    int nColsE = initialEfield->E[0][0]->n_cols;
    arrayGeom EGeom = initialEfield->arrayGeometry;
    arrayGeom zoomGeom;
    fpmScale = hexFPM->get_fpmScale();
    
    assert(fpmScale != -1);  // make sure init_mask was run first

    double *lambdaFocalLengthIn = new double[initialEfield->E[0][0]->n_slices];
    double *lambdaFocalLengthOut = new double[initialEfield->E[0][0]->n_slices];
    for (int sl=0; sl<initialEfield->E[0][0]->n_slices; sl++) {
        lambdaFocalLengthIn[sl] = initialEfield->lambdaData[sl].lambda*3.52;
        lambdaFocalLengthOut[sl] = -lambdaFocalLengthIn[sl];
    }
    
    // initialize the zoomed grid geometry
    
    zoomGeom.pixelSizeX = fpmScale;
    zoomGeom.pixelSizeY = fpmScale;
    int nRows = padFactor*nRowsE;
    int nCols = padFactor*nColsE;
    
    zoomGeom.pixelX = arma::zeros<arma::vec>(nCols);
    for (int i=0; i<nCols; i++)
        zoomGeom.pixelX[i] = ((double) i - ((double)nCols)/2.)*zoomGeom.pixelSizeX;
    zoomGeom.pixelY = arma::zeros<arma::vec>(nRows);
    for (int i=0; i<nRows; i++)
        zoomGeom.pixelY[i] = ((double) i - ((double)nRows)/2.)*zoomGeom.pixelSizeY;
    zoomGeom.physicalSize = zoomGeom.pixelX.n_elem*zoomGeom.pixelSizeX;
    zoomGeom.physicalSizeX = zoomGeom.pixelX.n_elem*zoomGeom.pixelSizeX;
    zoomGeom.physicalSizeY = zoomGeom.pixelY.n_elem*zoomGeom.pixelSizeY;

    zoomGeom.print("zoomGeom: ");
    
    propZoomFftIn.init(EGeom, zoomGeom, lambdaFocalLengthIn, initialEfield->E[0][0]->n_slices);
    propZoomFftOut.init(zoomGeom, EGeom, lambdaFocalLengthOut, initialEfield->E[0][0]->n_slices);
}

void doubleFraunhoferInterpFPM::set(std::string fieldName, const char *arg) {
    bool found = celem::set(fieldName, arg);
    if (fieldName == "doubleFraunhoferInterpFPM")
    ;
    else if (fieldName == "zoomFactor") {
        zoomFactor = atof(arg);
    } else if (fieldName == "padFactor") {
        padFactor = atof(arg);
    } else if (!found)
        std::cout << "!!! doubleFraunhoferInterpFPM bad set field name: " << fieldName << std::endl;
}

void doubleFraunhoferInterpFPM::print(const char *hdr) {
    std::cout << "doubleFraunhoferInterpFPM " << hdr << std::endl;
    std::cout << "zoomFactor = " << zoomFactor << std::endl;
    std::cout << "padFactor = " << padFactor << std::endl;
    hexFPM->print("FPM:");
}

void doubleFraunhoferInterpFPM::draw(const char *title) {
    char str[200];
    sprintf(str, "%s complexMask %s", title, name);
}


