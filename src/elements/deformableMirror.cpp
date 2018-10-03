//
//  deformableMirror.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "deformableMirror.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>

deformableMirror::deformableMirror() {
}

deformableMirror::deformableMirror(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

efield* deformableMirror::execute(efield* E, celem* prev, celem* next, double time) {
//    std::cout << "executed a deformableMirror " << name << std::endl;
    std::complex<double> i1(0, 1);
    if (mirrorMat.n_rows != E->E[0][0]->n_rows) {
        //        draw_mat(mirrorMat, "mask before downsampling");
        downsample_via_convolution(mirrorMat, mirrorMat, E->initArrayGeometry, E->arrayGeometry);
        //        draw_mat(mirrorMat, "mask after downsampling");
    }

    pre_execute(E, prev, next, time);
    
    for (int s=0; s<E->E.size(); s++) {
        for (int p=0; p<E->E[s].size(); p++) {
            for (int k=0; k<E->E[s][p]->n_slices; k++) {
                double lambda = E->lambdaData[k].lambda;
                E->E[s][p]->slice(k) = E->E[s][p]->slice(k) % exp((-mirrorSign*2*2*M_PI/lambda)*i1*mirrorMat); // element-wise multiplication
            }
        }
    }
    
    post_execute(E, prev, next, time);
    
    return E;
}

void deformableMirror::read_mirror_file(const char *filename) {
    load_mat(filename, mirrorMat);
}

void deformableMirror::read_actuator_file(const char *filename) {
    load_mat(filename, actuatorMat);
    nActuatorRows = actuatorMat.n_rows;
    nActuatorCols = actuatorMat.n_cols;
    std::cout << "read_actuator_file: nActuatorRows = " << nActuatorRows << ", nActuatorCols = " << nActuatorCols << std::endl;
}

void deformableMirror::init(initCommandSet*& cmdBlock) {
//    std::cout << "initing a deformableMirror" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    
    actuatorUpsampleFft = new fft;
    influenceFunctionFft = new fft;
    surfIfft = new ifft;
    
    compute_influence_function();
    compute_deformable_mirror_surface();
//    draw("init");
    
    post_init();
}

void deformableMirror::set(std::string fieldName, const char *arg) {

    bool found = celem::set(fieldName, arg);
    if (fieldName == "deformableMirror")
        ;
    else if (fieldName == "filename") {
        // arg is the .fits filename that contains the deformableMirror definition
        read_mirror_file(arg);
    }
    else if (fieldName == "actuatorFilename") {
        // arg is the .fits filename that contains the deformableMirror actuator data
        read_actuator_file(arg);
    }
    else if (fieldName == "sign") {
        double mSign = atof(arg);
        if (mSign == -1 | mSign == 1)
            mirrorSign = mSign;
        else
            std::cout << "!!!! illegal propagationSign, ignoring" << std::endl;
    }
    else if (fieldName == "sigma") {
        sigma = atof(arg);
    }
    else if (fieldName == "nActuatorRows") {
        nActuatorRows = atoi(arg);
    }
    else if (fieldName == "nActuatorCols") {
        nActuatorCols = atoi(arg);
    }
    else if (!found)
        std::cout << "!!! deformableMirror bad set field name: " << fieldName << std::endl;
}

void deformableMirror::get_optimization_data(const char *dataName, void *data) {
    if (!strcmp(dataName, "actuatorValues")) {
        arma::vec dataToReturn;
        dataToReturn.set_size(optNRows*optNCols);
        int dataCount = 0;
        for (int i=0; i<optNRows; i++)
            for (int j=0; j<optNCols; j++) {
                dataToReturn[dataCount] = actuatorMat(optr0 + i*optStride, optc0 + j*optStride);
                dataCount++;
            }
        *(arma::vec *)data = dataToReturn;
//        *(arma::vec *)data = vectorize(actuatorMat);
    }
}

void deformableMirror::set_optimization_data(const char *dataName, void *data) {
    if (!strcmp(dataName, "actuatorValues")) {
        arma::vec returnedData = *(arma::vec *)data;
        int dataCount = 0;
        for (int i=0; i<optNRows; i++)
            for (int j=0; j<optNCols; j++) {
                actuatorMat(optr0 + i*optStride, optc0 + j*optStride) = returnedData[dataCount];
                dataCount++;
            }
        //        fpmSags = *(arma::vec *)data;
        compute_deformable_mirror_surface();
        //        draw_mat(mirrorMat, "mirrorMat", "matlab");
    }
}

void deformableMirror::save_optimization_data(const char *dataName, char *outputDirectory) {
    if (!strcmp(dataName, "actuatorValues")) {
        
        std::string fname = (std::string)outputDirectory + "/optimalActuators.fits";
        save_mat(fname.c_str(), actuatorMat);
    }
}

void deformableMirror::compute_influence_function(void) {
    assert(initialEfield);
    assert(nActuatorRows > 0);
    
    // assumes the actuator array is square
    double influenceSigma = sigma*initialEfield->arrayGeometry.physicalSize/nActuatorRows;
    
    influenceFunction = exp(-4*log(2)*(square(initialEfield->arrayGeometry.pixelXX) + square(initialEfield->arrayGeometry.pixelYY))/(influenceSigma*influenceSigma));
//    std::cout << "influenceFunction size: " << size(influenceFunction) << std::endl;
}

void deformableMirror::compute_deformable_mirror_surface(void) {
    std::complex<double> i1(0, 1);

    assert(nActuatorRows > 0);
    assert(nActuatorCols > 0);

    // upsample the actuator data by evenly placing the data as individual points in an array of the same size as the influence function
    int upSampleRowFactor = ceil((influenceFunction.n_rows + 1)/nActuatorRows);
    int upSampleColFactor = ceil((influenceFunction.n_cols + 1)/nActuatorCols);

//    std::cout << "compute_deformable_mirror_surface: upSampleRowFactor = " << upSampleRowFactor << ", upSampleColFactor = " << upSampleColFactor << std::endl;
    
    actuatorUpsample = arma::zeros<arma::mat>(upSampleRowFactor*nActuatorRows, upSampleColFactor*nActuatorCols);
//    std::cout << "actuatorUpsample size: " << size(actuatorUpsample) << std::endl;
//    draw_mat(actuatorMat, "actuatorMat");
    
    for (int r=0; r<actuatorMat.n_rows; r++) {
        for (int c=0; c<actuatorMat.n_cols; c++) {
            actuatorUpsample(r*upSampleRowFactor,c*upSampleColFactor) = actuatorMat(r,c);
        }
    }
    FTactuatorUpsample = arma::zeros<arma::cx_mat>(size(actuatorUpsample));
    FTactuatorUpsample.set_real(actuatorUpsample);
    actuatorUpsampleFft->execute(FTactuatorUpsample);
//    save_mat("FTactuatorUpsample_fft.fits", FTactuatorUpsample);

    // account for the fact that the actuators operate in the center of cells by shifting half a cell in space = multiplying by exponential phase factor
//    std::cout << "size rowvecs: " << size(arma::linspace<arma::rowvec>(0., (double) ceil(FTactuatorUpsample.n_rows/2)-1, ceil(FTactuatorUpsample.n_rows/2))) << ", " << size(arma::linspace<arma::rowvec>(-(double) ceil(FTactuatorUpsample.n_rows/2)-1, -1.,  ceil(FTactuatorUpsample.n_rows/2))) << std::endl;
    arma::rowvec fx = arma::join_rows(arma::linspace<arma::rowvec>(0., (double) ceil(FTactuatorUpsample.n_rows/2)-1, ceil(FTactuatorUpsample.n_rows/2)), arma::linspace<arma::rowvec>(-(double) ceil(FTactuatorUpsample.n_rows/2), -1.,  ceil(FTactuatorUpsample.n_rows/2)))/influenceFunction.n_rows;
//    std::cout << "size colvecs: " << size(arma::linspace<arma::colvec>(0., (double) ceil(FTactuatorUpsample.n_cols/2)-1, ceil(FTactuatorUpsample.n_cols/2))) << ", " << size(arma::linspace<arma::colvec>(-(double) ceil(FTactuatorUpsample.n_cols/2)-1, -1., ceil(FTactuatorUpsample.n_cols/2))) << std::endl;
    arma::colvec fy = arma::join_cols(arma::linspace<arma::colvec>(0., (double) ceil(FTactuatorUpsample.n_cols/2)-1, ceil(FTactuatorUpsample.n_cols/2)), arma::linspace<arma::colvec>(-(double) ceil(FTactuatorUpsample.n_cols/2), -1., ceil(FTactuatorUpsample.n_cols/2)))/influenceFunction.n_cols;
    
    arma::mat fxs = arma::ones(FTactuatorUpsample.n_rows, 1)*fx;
    arma::mat fys = fy*arma::ones(1, FTactuatorUpsample.n_cols);
    
    arma::mat phaseMat = fxs*(influenceFunction.n_rows*(1-1/nActuatorRows)/2) + fys*(influenceFunction.n_cols*(1-1/nActuatorCols)/2);
//    save_mat("phaseMat.fits", phaseMat);
    FTactuatorUpsample = FTactuatorUpsample % exp(2*M_PI*i1*phaseMat);
//    save_mat("FTactuatorUpsample_phaseShifted.fits", FTactuatorUpsample);

//    std::cout << "truncating, size(FTactuatorUpsample): " << size(FTactuatorUpsample) << std::endl;
    arma::cx_mat FTactuatorTruncated = fft_truncate(FTactuatorUpsample, influenceFunction.n_rows, influenceFunction.n_cols);
//    std::cout << "done truncating, size(FTactuatorTruncated): " << size(FTactuatorTruncated) << std::endl;

    FTinfluenceFunction = arma::zeros<arma::cx_mat>(size(influenceFunction));
    FTinfluenceFunction.set_real(influenceFunction);
    influenceFunctionFft->execute(FTinfluenceFunction);
//    save_mat("FTinfluenceFunction.fits", FTinfluenceFunction);
    arma::cx_mat FTMirrorMat = FTactuatorTruncated % FTinfluenceFunction;
//    save_mat("truncatedFTMirrorMat.fits", FTinfluenceFunction);

    surfIfft->execute(FTMirrorMat);
    mirrorMat = arma::real(FTMirrorMat);
//    save_mat("mirrorMat.fits", mirrorMat);
}

void deformableMirror::draw(const char *title) {
    char str[200];
    sprintf(str, "%s deformableMirror %s", title, name);
    draw_mat(mirrorMat, str, "gray");
}

