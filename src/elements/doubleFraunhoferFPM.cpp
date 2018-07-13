//
//  doubleFraunhoferFPM.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "doubleFraunhoferFPM.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include "../coronagraph/coronagraph.hpp"

doubleFraunhoferFPM::doubleFraunhoferFPM() {
}

doubleFraunhoferFPM::doubleFraunhoferFPM(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

efield* doubleFraunhoferFPM::execute(efield* E, celem* prev, celem* next, double time) {
    std::complex<double> i1(0, 1);
    
    
    std::cout << "calibrating: " << globalCoronagraph->get_calibration_state() << ", disableForCalibration: " << disableForCalibration << std::endl;
    if (!(disableForCalibration && globalCoronagraph->get_calibration_state())) {
        int nRowsE = E->E[0][0]->n_rows;
        int nColsE = E->E[0][0]->n_cols;
        arrayGeom EGeom = E->arrayGeometry;
        arrayGeom zoomGeom;
        
        std::cout << "executing doubleFraunhoferFPM " << name << std::endl;
        E->arrayGeometry.print("doubleFraunhoferFPM input E geometry: ");
        E->print("at top of doubleFraunhoferFPM execute:");
        
        pre_execute(E, prev, next, time);

        for (int s=0; s<E->E.size(); s++) {
            for (int p=0; p<E->E[s].size(); p++) {
                double *lambdaFocalLengthIn = new double[E->E[0][0]->n_slices];
                double *lambdaFocalLengthOut = new double[E->E[0][0]->n_slices];
                for (int sl=0; sl<E->E[s][p]->n_slices; sl++) {
                    lambdaFocalLengthIn[sl] = E->lambdaData[sl].lambda*3.52;
                    lambdaFocalLengthOut[sl] = -lambdaFocalLengthIn[sl];
                }

                for (int sl=0; sl<E->E[s][p]->n_slices; sl++) {
                    // initialize the zoomed grid geometry
                    double dxprime = 2*(E->beamRadiusPhysical/(E->pixelScale*nRowsE*zoomFactor))*E->lambdaData[sl].lambda*80;
                    std::cout << "lambda = " << E->lambdaData[sl].lambda << ", dxprime = " << dxprime << std::endl;
                    
                    zoomGeom.pixelSizeX = dxprime;
                    zoomGeom.pixelSizeY = dxprime;
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
                    if (sl == 0)
                        zoomGeom.print("zoomGeom: ");
                    
                    propZoomFftIn.init(EGeom, zoomGeom, lambdaFocalLengthIn, E->E[0][0]->n_slices);
                    propZoomFftOut.init(zoomGeom, EGeom, lambdaFocalLengthOut, E->E[0][0]->n_slices);
                
                
                    // the actual execution
                    arma::cx_mat tcMat;

                    std::cout << "executing doubleFraunhoferFPM slice " << sl << std::endl;

                    // Fraunhofer propagation 1
                    tcMat = propZoomFftIn.execute(E->E[s][p]->slice(sl), sl);
                    
                    // apply the complex mask
                    int row0 = tcMat.n_rows/2 - complexMaskCube.n_rows/2;
                    int col0 = tcMat.n_cols/2 - complexMaskCube.n_cols/2;
                    if (maskIndex == -1)
                        maskIndex = sl;
                    std::cout << "maskIndex = " << maskIndex << std::endl;

#pragma omp parallel
                    {
#pragma omp for
                        for (int c=0; c<complexMaskCube.n_cols; ++c) {
                            for (int r=0; r<complexMaskCube.n_rows; ++r)
                                tcMat(row0 + r, col0 + c) *= complexMaskCube(r,c,maskIndex);
                        }
                    }
                    
                    // Fraunhofer propagation 2
                    E->E[s][p]->slice(sl) = propZoomFftOut.execute(tcMat, sl);
                }
            }
        }
        E->print("post doubleFraunhoferFPM E: ");
        
        post_execute(E, prev, next, time);
    }
    
    return E;
}

void doubleFraunhoferFPM::init(initCommandSet*& cmdBlock) {
    std::cout << "initing a doubleFraunhoferFPM" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    post_init();
}

void doubleFraunhoferFPM::initMask(const char *filenameAmp, const char *filenamePh) {
    std::complex<double> i1(0, 1);
    
    load_cube(filenameAmp, complexMaskMatAmp);
    load_cube(filenamePh, complexMaskMatPh);
    complexMaskCube.set_size(size(complexMaskMatAmp));
    for (int sl=0; sl<complexMaskMatAmp.n_slices; ++sl) {
        complexMaskCube.slice(sl) = complexMaskMatAmp.slice(sl) % exp(-i1*complexMaskMatPh.slice(sl));
    }
}

void doubleFraunhoferFPM::set(std::string fieldName, const char *arg) {
    bool found = celem::set(fieldName, arg);
    if (fieldName == "doubleFraunhoferFPM")
    ;
    else if (fieldName == "zoomFactor") {
        zoomFactor = atof(arg);
    } else if (fieldName == "padFactor") {
        padFactor = atof(arg);
    } else if (fieldName == "maskIndex") {
        // arg is a single integer
        maskIndex = atoi(arg);
    } else if (fieldName == "maskFilename") {
        // arg is two filenames separated by a comma, each giving
        // the .fits amplitude and phase filename that contains the complexMask definition
        const char *cPtr = strchr(arg, ',');
        int strAmpLen = cPtr - arg;
        int strPhLen = strlen(arg) - strAmpLen + 1;
        std::cout << "strAmpLen " << strAmpLen << ", strPhLen " << strPhLen << std::endl;
        char *strAmp = new char[strAmpLen + 1];
        char *strPh = new char[strPhLen + 1];
        strncpy(strAmp, arg, strAmpLen);
        strAmp[strAmpLen] = '\0';
        strncpy(strPh, cPtr + 1, strPhLen);
        strPh[strPhLen] = '\0';
        std::cout << "loading " << strAmp << " and " << strPh << std::endl;
        initMask(strAmp, strPh);
        //        draw("initial FPM");
    } else if (!found)
        std::cout << "!!! doubleFraunhoferFPM bad set field name: " << fieldName << std::endl;
}

void doubleFraunhoferFPM::print(const char *hdr) {
    std::cout << "doubleFraunhoferFPM " << hdr << std::endl;
    std::cout << "zoomFactor = " << zoomFactor << std::endl;
    std::cout << "padFactor = " << padFactor << std::endl;
}

void doubleFraunhoferFPM::draw(const char *title) {
    char str[200];
    sprintf(str, "%s complexMask %s", title, name);
    draw_mat(complexMaskMatAmp.slice(maskIndex), str, "rainbow");
    draw_mat(complexMaskMatPh.slice(maskIndex), str, "rainbow");
}


