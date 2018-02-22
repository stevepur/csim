//
//  celem.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "complexMask.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include "../lib/csim_lib.hpp"
#include "../coronagraph/coronagraph.hpp"

complexMask::complexMask() {
}

complexMask::complexMask(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

efield* complexMask::execute(efield* E, celem* prev, celem* next, double time) {
    std::complex<double> i1(0, 1);
    pre_execute(E, prev, next, time);
    
    std::cout << "calibrating: " << globalCoronagraph->get_calibration_state() << ", disableForCalibration: " << disableForCalibration << std::endl;
    if (!(disableForCalibration && globalCoronagraph->get_calibration_state())) {
        std::cout << "executed complexMask " << name << std::endl;
        fflush(stdout);
        for (int s=0; s<E->E.size(); s++) {
            for (int p=0; p<E->E[s].size(); p++) {
                for (int k=0; k<E->E[s][p]->n_slices; k++) {
                    arma::cx_mat *tcMat = &E->E[s][p]->slice(k);
                    int row0 = tcMat->n_rows/2 - complexMaskCube.n_rows/2;
                    int col0 = tcMat->n_cols/2 - complexMaskCube.n_cols/2;
#pragma omp parallel
                    {
#pragma omp for
                        for (int c=0; c<complexMaskCube.n_cols; ++c) {
                            for (int r=0; r<complexMaskCube.n_rows; ++r)
                                (*tcMat)(row0 + r, col0 + c) *= complexMaskCube(r,c,maskIndex);
                        }
                    }
//                    E->E[s][p]->slice(k) = tcMat;
                }
            }
        }
    }
    
    post_execute(E, prev, next, time);
    
    return E;
}

void complexMask::init(const char *filenameAmp, const char *filenamePh) {
    std::complex<double> i1(0, 1);

    load_cube(filenameAmp, complexMaskMatAmp);
    load_cube(filenamePh, complexMaskMatPh);
    complexMaskCube.set_size(size(complexMaskMatAmp));
    for (int sl=0; sl<complexMaskMatAmp.n_slices; ++sl) {
        complexMaskCube.slice(sl) = complexMaskMatAmp.slice(sl) % exp(-i1*complexMaskMatPh.slice(sl));
    }
}

void complexMask::init(initCommandSet*& cmdBlock) {
    std::cout << "initing a complexMask" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    post_init();
}

void complexMask::set(std::string fieldName, const char *arg) {
    std::cout << "complexMask::set: " << fieldName << ":" << arg << std::endl;
    bool found = celem::set(fieldName, arg);
    if (fieldName == "complexMask")
        ;
    else if (fieldName == "maskIndex") {
        // arg is a single integer
        maskIndex = atoi(arg);
    } else if (fieldName == "filename") {
        // arg is two filenames separated by a comma, each giving
        // the .fits amplitude and phase filename that contains the complexMask definition
        char *cPtr = strchr(arg, ',');
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
        init(strAmp, strPh);
//        draw("initial FPM");
    } else if (!found)
        std::cout << "!!! complexMask bad set field name: " << fieldName << std::endl;
}

void complexMask::draw(const char *title) {
    char str[200];
    sprintf(str, "%s complexMask %s", title, name);
    draw_mat(complexMaskMatAmp.slice(maskIndex), str, "rainbow");
    draw_mat(complexMaskMatPh.slice(maskIndex), str, "rainbow");
}

