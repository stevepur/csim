//
//  regionContrast.cpp
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
#include "regionContrast.hpp"
#include "../coronagraph/coronagraph.hpp"
#include "../telescope/telescope.hpp"
#include "../lib/csim_lib.hpp"

regionContrast::regionContrast() {
}

regionContrast::regionContrast(initCommandSet*& cmdBlocks) {
    init(cmdBlocks);
}

void regionContrast::init(initCommandSet*& cmdBlock) {
//    std::cout << "initing a regionContrast" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
}

void regionContrast::set(std::string fieldName, const char *arg) {
    if (fieldName == "regionContrast")
        ;
    else if (fieldName == "filename") {
        // arg is a string
        filename = new char[strlen(arg)+1];
        strcpy(filename, arg);
    }
    else if (fieldName == "referenceLambda") {
        // arg is a single float
        referenceLambda = atof(arg);
        set_loD();
    }
    else if (fieldName == "anglesInDegrees") {
        float v1, v2;
        // arg is a comma separated pair of doubles
        sscanf(arg, "%f, %f", &v1, &v2);
        angle1 = (double) v1;
        angle2 = (double) v2;
        std::cout << "region angles: " << angle1 << " to " << angle2 << std::endl;
    }
    else if (fieldName == "radii") {
        float v1, v2;
        // arg is a comma separated pair of doubles
        sscanf(arg, "%f, %f", &v1, &v2);
        radius1 = (double) v1;
        radius2 = (double) v2;
        std::cout << "region radii: " << radius1 << " to " << radius2 << std::endl;
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
            std::cout << "!!! regionContrast draw bad set field name: " << fieldName << std::endl;
    }
    else
        std::cout << "!!! regionContrast bad set field name: " << fieldName << std::endl;
    
}

void regionContrast::get_region_pixels(efield *E, arma::uvec& pixelIndex, arma::vec& pixelX, arma::vec& pixelY) {
    
    arma::umat inSample = (E->arrayGeometry.pixelRR/loD > radius1)
    % (E->arrayGeometry.pixelRR/loD < radius2)
    % (E->arrayGeometry.pixelTT > angle1*M_PI/180)
    % (E->arrayGeometry.pixelTT < angle2*M_PI/180);
    
    pixelIndex = find(inSample);
    pixelX = E->arrayGeometry.pixelXX(pixelIndex);
    pixelY = E->arrayGeometry.pixelYY(pixelIndex);
}

void regionContrast::set_loD(void) {
    
    // if referenceLambda is not defined, set it to the average of lambdas in E
    if (referenceLambda == -1.0) {
        int nLambdas = initialEfield->E[0][0]->n_slices;
        double lambdaSum = 0.0;
        for (int i=0; i<nLambdas; i++)
            lambdaSum += initialEfield->lambdaData[i].get_wavelength();
        referenceLambda = lambdaSum/nLambdas;
    }
    loD = globalTelescope->get("primaryfLength")*referenceLambda/globalTelescope->get("primaryDiameter");
    std::cout << "referenceLambda = " << referenceLambda << ", loD = " << loD << std::endl;
}

void regionContrast::compute_contrast(void) {
    arma::wall_clock timer;
    
    efield *calibEfield = new efield(*initialEfield);
    calibEfield->set("name", "calibration Efield");
    efield *fullEfield = new efield(*initialEfield);
    fullEfield->set("name", "FPM Efield");
    
    globalCoronagraph->set_calibration_state(true);
    timer.tic();
    globalCoronagraph->execute(calibEfield, 0);
    std::cout << "contrast curve calibration csim execution time: " << timer.toc() << " seconds" << std::endl;
    
    globalCoronagraph->set_calibration_state(false);
    timer.tic();
    globalCoronagraph->execute(fullEfield, 0);
    std::cout << "contrast curve csim execution time: " << timer.toc() << " seconds" << std::endl;
    
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
    
    arma::uvec pixelIndex;
    arma::vec pixelX;
    arma::vec pixelY;
    get_region_pixels(fullEfield, pixelIndex, pixelX, pixelY);

    if (draw) {
        arma::umat matIndex = arma::ind2sub(size(fullIntensitySum), pixelIndex);
        int minRow = min(matIndex(0,arma::span::all));
        int maxRow = max(matIndex(0,arma::span::all));
        int minCol = min(matIndex(1,arma::span::all));
        int maxCol = max(matIndex(1,arma::span::all));
        
        std::cout << "bounding indices: " << minRow << ", " << maxRow << ", " << minCol << ", " << maxCol << std::endl;
        
        arma::mat drawContrastMat = zeros(arma::size(fullIntensitySum));
        for (int i=0; i<pixelIndex.n_elem; i++)
            drawContrastMat(pixelIndex(i)) = 1.0;
        drawContrastMat = drawContrastMat % fullIntensitySum;
        arma::mat contrastRegionExtract = drawContrastMat(arma::span(minRow,maxRow),arma::span(minCol,maxCol));
        draw_mat(log10(contrastRegionExtract), arma::min(pixelX/loD), arma::max(pixelX/loD), arma::min(pixelY/loD), arma::max(pixelY/loD), "normalized full coronagraph PSF", "matlab");
    }
    
    
    save_mat("normalized_PSF.fits", fullIntensitySum);
    
    double contrast = arma::mean(fullIntensitySum(pixelIndex));
    std::cout << "region contrast: " << contrast << std::endl;

    assert(filename != NULL);
    FILE *fid = fopen(filename, "w");
    fprintf(fid, "%f, %f\n", angle1, angle2);
    fprintf(fid, "%f, %f\n", radius1, radius2);
    fprintf(fid, "%e\n", contrast);
    fclose(fid);
}

void regionContrast::print(const char *hdr) {
    if (hdr != NULL)
        std::cout << hdr << ":" << std::endl;
    std::cout << "referenceLambda: " << referenceLambda << std::endl;
    std::cout << "region angles: " << angle1 << " to " << angle2 << std::endl;
    std::cout << "region radii: " << radius1 << " to " << radius2 << std::endl;
}
