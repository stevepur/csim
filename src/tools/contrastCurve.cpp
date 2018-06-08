//
//  contrastCurve.cpp
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
#include "contrastCurve.hpp"
#include "../coronagraph/coronagraph.hpp"
#include "../telescope/telescope.hpp"
#include "../lib/csim_lib.hpp"

contrastCurve::contrastCurve() {
}

contrastCurve::contrastCurve(initCommandSet*& cmdBlocks) {
    init(cmdBlocks);
}

void contrastCurve::init(initCommandSet*& cmdBlock) {
//    std::cout << "initing a contrastCurve" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
}

void contrastCurve::set(std::string fieldName, const char *arg) {
    if (fieldName == "makeContrastCurve")
        ;
    else if (fieldName == "filename") {
        // arg is a string
        filename = new char[strlen(arg)+1];
        strcpy(filename, arg);
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
            std::cout << "!!! contrastCurve draw bad set field name: " << fieldName << std::endl;
    }
    else if (fieldName == "drawTo") {
        // arg is a string
        if (!strcmp(arg, "full")) {
            drawTo = -1;
        }
        else {
            drawTo = atof(arg);
        }
    }
    else
        std::cout << "!!! contrastCurve bad set field name: " << fieldName << std::endl;
}

void contrastCurve::make_contrast_curve(void) {
    arma::wall_clock timer;
    
    efield *calibEfield = new efield(*initialEfield);
    calibEfield->set("name", "calibration Efield");
    efield *fullEfield = new efield(*initialEfield);
    fullEfield->set("name", "FPM Efield");
    
//    initialEfield->print("init efield");
//    calibEfield->print("calibration Efield");
//    fullEfield->print("FPM Efield");
    
    globalCoronagraph->set_calibration_state(true);
    timer.tic();
    globalCoronagraph->execute(calibEfield, 0);
    std::cout << "contrast curve calibration csim execution time: " << timer.toc() << " seconds" << std::endl;
    
    globalCoronagraph->set_calibration_state(false);
    timer.tic();
    globalCoronagraph->execute(fullEfield, 0);
    std::cout << "contrast curve csim execution time: " << timer.toc() << " seconds" << std::endl;
    
    double lambda = calibEfield->lambdaData[(int) round(calibEfield->E[0][0]->n_slices/2)].lambda;
    double loD = globalTelescope->compute_loD(lambda);
    std::cout << "lambda = " << lambda << ", loD = " << loD << std::endl;
    
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
//    if (draw) {
//        arma::vec calPx = calibEfield->arrayGeometry.pixelX/loD;
//        arma::vec calPy = calibEfield->arrayGeometry.pixelY/loD;
//        draw_mat(log10(calibIntensitySum), calPx[0], calPx[calPx.n_elem-1], calPy[0], calPy[calPy.n_elem-1], "calibration PSF", "matlab");
//    }
    save_mat("calibration_PSF.fits", calibIntensitySum);
    
    arma::cube fullIntensity;
    fullIntensity.zeros(size(*(calibEfield->E[0][0])));
    for (int s=0; s<fullEfield->E.size(); s++) {
        for (int p=0; p<fullEfield->E[s].size(); p++) {
            fullIntensity += real(*(fullEfield->E[s][p]) % arma::conj(*(fullEfield->E[s][p])));
        }
    }
    arma::vec fullPx = fullEfield->arrayGeometry.pixelX/loD;
    arma::vec fullPy = fullEfield->arrayGeometry.pixelY/loD;
    arma::mat fullIntensitySum = sum(fullIntensity, 2)/calibMaxIntensity;
    std::cout << "fullMaxIntensity = " << max(max(fullIntensitySum)) << std::endl;
    if (draw) {
        draw_mat(log10(fullIntensitySum), fullPx[0], fullPx[fullPx.n_elem-1], fullPy[0], fullPy[fullPy.n_elem-1], "normalized full coronagraph PSF", "matlab");
    }
    
    save_mat("normalized_PSF.fits", fullIntensitySum);
    arma::vec ccSamplePoints = arma::linspace<arma::vec>(0., arma::max(fullPx), (int) ((fullPx.n_elem-2)/pixelSampling));
    arma::vec ccAvg = arma::zeros<arma::vec>(ccSamplePoints.n_elem-2);
    arma::vec ccPos = arma::zeros<arma::vec>(ccSamplePoints.n_elem-2);
    arma::mat RR = fullEfield->arrayGeometry.pixelRR/loD;
    
    std::cout << "computing contrast curve" << std::endl;
    timer.tic();
    #pragma omp parallel for
    for (int i=1; i<ccSamplePoints.n_elem-1; i++) {
        arma::umat inSample = (RR > ccSamplePoints[i-1]) % (RR < ccSamplePoints[i+1]);
        ccAvg[i-1] = arma::mean(fullIntensitySum(find(inSample)));
        ccPos[i-1] = ccSamplePoints[i];
    }
    std::cout << "execution time: " << timer.toc() << " seconds" << std::endl;
    std::cout << "finished computing contrast curve" << std::endl;
    
    assert(filename != NULL);
    FILE *fid = fopen(filename, "w");
    for (int i=1; i<ccAvg.n_elem; i++)
        fprintf(fid, "%f %e\n", ccPos[i], ccAvg[i]);
    fclose(fid);

    if (draw) {
        if (drawTo < 0) {
            plot_vec(ccPos, ccAvg, "contrast curve", "semilogy");
        }
        else {
            plot_vec(ccPos, ccAvg, 0, drawTo, -1, -1, "contrast curve", "semilogy");
        }
    }
    
}


