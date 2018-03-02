//
//  fpmPupToLyot.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "fpmPupToLyot.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include "../telescope/telescope.hpp"
#include "../coronagraph/coronagraph.hpp"

fpmPupToLyot::fpmPupToLyot() {
}

fpmPupToLyot::fpmPupToLyot(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

efield* fpmPupToLyot::execute(efield* E, celem* prev, celem* next, double time) {
//    std::cout << "executed a fpmPupToLyot " << name << std::endl;
    std::complex<double> i1(0, 1);
    std::complex<double> r1(1, 0);
    int nRowsE = E->E[0][0]->n_rows;
    int nColsE = E->E[0][0]->n_cols;

    pre_execute(E, prev, next, time);
    
    assert(E->beamRadiusPhysical > 0);
    
    // create space for a zero-padded E array, with a 2x padding
    paddedE = arma::zeros<arma::cx_mat>(2*nRowsE, 2*nColsE);
    paddedGeom.set_xy(paddedE.n_rows, paddedE.n_cols, E->arrayGeometry.physicalSize); // last arg is radius, so we divide by two but we're doubling
    fftMHatTimesfftPaddedE = arma::zeros<arma::cx_mat>(size(paddedE));
    double *lambda = new double[E->E[0][0]->n_slices];
    double *lambdaFocalLength = new double[E->E[0][0]->n_slices];
    mask->set_geometry(this, E, lambda, lambdaFocalLength);
    propZoomFft.init(maskGeom, paddedGeom, lambdaFocalLength, E->E[0][0]->n_slices);
//    std::cout << "lambda = " << lambda[0] << ", lambdaFocalLength = " << lambdaFocalLength[0] << std::endl;
//    maskGeom.print("maskGeom");
    for (int sl=0; sl<E->E[0][0]->n_slices; sl++) {
//        std::cout << "fpmPupToLyot::execute slice " << sl << ", lambda = " << lambda[sl] << std::endl;
        
//        std::cout << "lambda = " << lambda[sl] << ", lambdaFocalLength = " << lambdaFocalLength[sl] << std::endl;
//        std::cout << "calibrating: " << globalCoronagraph->get_calibration_state() << ", disableForCalibration: " << disableForCalibration << std::endl;
        mask->set_fpmMatAmp(this, lambda[sl], sl);
        // matlab: mask.M_hat = zoomFFT_realunits(mask.x, mask.y, mask.M - 1, pupil_ext.x, pupil_ext.y, pupil.f, lambda);
//        if (sl == E->E[0][0]->n_slices-1)
//            save_mat("fpmMatAmp.fits", fpmMatAmp, "reIm");
        maskHat = zoomFftSign*propZoomFft.execute(fpmMatAmp, sl);
//        if (sl == E->E[0][0]->n_slices-1)
//            save_mat("maskHat.fits", maskHat, "reIm");
        
        // matlab: FFT_M_hat = fft2((fftshift(mask.M_hat)));
        fftMaskHat = fft_shift(maskHat);
        maskFft.execute(fftMaskHat);
//        if (sl == E->E[0][0]->n_slices-1)
//            save_mat("fftMaskHat.fits", fftMaskHat, "reIm");
        for (int s=0; s<E->E.size(); s++) {
            for (int p=0; p<E->E[s].size(); p++) {
                // embed E in the zero-padded cube
                paddedE(arma::span(nRowsE/2, 3*nRowsE/2-1), arma::span(nColsE/2, 3*nColsE/2-1)) = E->E[s][p]->slice(sl);
                
                // matlab: FFT_pupil = fft2(pupil_ext.E);
                fftPaddedE = paddedE;
                paddedEFft.execute(fftPaddedE);
//                if (sl == E->E[0][0]->n_slices-1)
//                    save_mat("fftPaddedE.fits", fftPaddedE, "reIm");

                // matlab: E = ifft2(FFT_M_hat.*FFT_pupil)/(-1i*lambda*pupil.f)*(pupil.x(2) - pupil.x(1))^2 + pupil_ext.E;
                fftMHatTimesfftPaddedE = fftMaskHat%fftPaddedE;
                myIfft.execute(fftMHatTimesfftPaddedE);
//                if (sl == E->E[0][0]->n_slices-1)
//                    save_mat("fftMHatTimesfftPaddedE.fits", fftMHatTimesfftPaddedE, "reIm");
                fftMHatTimesfftPaddedE = fftMHatTimesfftPaddedE*pow(paddedGeom.pixelSizeX, 2)/(-i1*fRatioSign*lambdaFocalLength[sl]);
//                if (sl == E->E[0][0]->n_slices-1)
//                    save_mat("fftMHatTimesfftPaddedE_mult.fits", fftMHatTimesfftPaddedE, "reIm");
                mask->apply_babinet(this);
//                if (sl == E->E[0][0]->n_slices-1)
//                    save_mat("fftMHatTimesfftPaddedE_bab.fits", fftMHatTimesfftPaddedE, "reIm");
                
                // extract the central result for the final answer
                E->E[s][p]->slice(sl) = fftMHatTimesfftPaddedE(arma::span(nRowsE/2, 3*nRowsE/2-1), arma::span(nColsE/2, 3*nColsE/2-1));
            }
        }
    }
    post_execute(E, prev, next, time);
//    E->print("exiting fpmPupToLyot::execute: ");

    return E;
}

void fpmPupToLyot::init(initCommandSet*& cmdBlock) {
    //    std::cout << "initing a fpmPupToLyot" << std::endl;
//    cmdBlock->print("----- fpmPupToLyot command block");
    std::vector<initCommandSet*> subBlocks = cmdBlock->find_command_blocks();
    for (int j=0; j<subBlocks.size(); j++) {
//        char str[200];
//        sprintf(str, ">>> block %d", j);
//        subBlocks[j]->print(str);
        
        if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "fpmPupToLyot")) {
            init_block(subBlocks[j]);
        }
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "fpmCMCForPupToLyot")) {
            mask = new fpmCMCForPupToLyot(this, subBlocks[j]);
        }
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "fpmIntHexCMCForPupToLyot")) {
            mask = new fpmIntHexCMCForPupToLyot(this, subBlocks[j]);
        }
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "fpmBinaryForPupToLyot")) {
            mask = new fpmBinaryForPupToLyot(this, subBlocks[j]);
        }
    }
    
//    print("fpmPupToLyot");
    post_init();
}

void fpmPupToLyot::init_block(initCommandSet*& cmdBlock) {
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
}

void fpmPupToLyot::set(std::string fieldName, const char *arg) {

    bool found = celem::set(fieldName, arg);
    if (fieldName == "fpmPupToLyot") {
        ;
    } else if (fieldName == "focalRatio") {
        // arg is a single double
        focalRatio = atof(arg);
        std::cout << "focalRatio = " << focalRatio << std::endl;
    } else if (fieldName == "fRatioSign") {
        // arg is a single double
        fRatioSign = atof(arg);
        std::cout << "fRatioSign = " << fRatioSign << std::endl;
    } else if (fieldName == "zoomFftSign") {
        // arg is a single double
        zoomFftSign = atof(arg);
        std::cout << "zoomFftSign = " << zoomFftSign << std::endl;
    } else if (!found)
        std::cout << "!!! fpmPupToLyot bad set field name: " << fieldName << std::endl;
}

void fpmPupToLyot::get_optimization_data(const char *dataName, void *data) {
    mask->get_optimization_data(dataName, data);
}

void fpmPupToLyot::set_optimization_data(const char *dataName, void *data){
    mask->set_optimization_data(dataName, data);
}

void fpmPupToLyot::draw(const char *title) {
    mask->draw(title);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

fpmCMCForPupToLyot::fpmCMCForPupToLyot(fpmPupToLyot *p2l, initCommandSet*& cmdBlock) {
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(p2l, cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
}

void fpmCMCForPupToLyot::set(fpmPupToLyot *p2l, std::string fieldName, const char *arg) {
    
    if (fieldName == "fpmCMCForPupToLyot")
        ;
    else if (fieldName == "maskIndex") {
        // arg is a single integer
        maskIndex = atoi(arg);
    } else if (fieldName == "maskFilename") {
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
        initMask(strAmp, strPh);
        //        draw("initial FPM");
    } else if (fieldName == "fpmMatScale") {
        // arg is a single double
        fpmMatScale = atof(arg);
        std::cout << "fpmMatScale = " << fpmMatScale << std::endl;
    } else
        std::cout << "!!! fpmCMCForPupToLyot bad set field name: " << fieldName << std::endl;
}


void fpmCMCForPupToLyot::initMask(const char *filenameAmp, const char *filenamePh) {
    std::complex<double> i1(0, 1);
    
    load_cube(filenameAmp, complexMaskMatAmp);
    load_cube(filenamePh, complexMaskMatPh);
    complexMaskCube.set_size(size(complexMaskMatAmp));
    for (int sl=0; sl<complexMaskMatAmp.n_slices; ++sl) {
        complexMaskCube.slice(sl) = complexMaskMatAmp.slice(sl) % exp(-i1*complexMaskMatPh.slice(sl));
    }
}

void fpmCMCForPupToLyot::set_geometry(fpmPupToLyot *p2l, efield* E, double *lambda, double *lambdaFocalLength) {
    // matlab:
    // mask.dx = mask.x(2) - mask.x(1);
    // mask.dxprime = 1.098e-06; (fpmMatScale)
    // mask.resample = mask.dx / mask.dxprime;
    // mask.x = mask.x / mask.resample;
    p2l->maskGeom.set_xy(complexMaskCube.n_rows, complexMaskCube.n_cols, (E->arrayGeometry.physicalSize/2.)/(E->arrayGeometry.pixelSizeX/fpmMatScale));
    
    //    E->print("in set_geometry");
    double cFRatio = p2l->focalRatio*2*E->beamRadiusPhysical;
//    std::cout << "cFRatio set to: " << cFRatio << std::endl;
    for (int i=0; i<E->E[0][0]->n_slices; i++) {
        // lambdaFocalLength[i] = E->lambdaData[i].lambda*globalTelescope->get("primaryfRatio");
        //        std::cout << "lambda from E: " << E->lambdaData[i].lambda << std::endl;
        lambda[i] = E->lambdaData[i].lambda;
        //        std::cout << "lambda set to: " << lambda[i] << " for i = " << i << std::endl;
        
        lambdaFocalLength[i] = p2l->fRatioSign*lambda[i]*cFRatio;
    }
}

void fpmCMCForPupToLyot::set_fpmMatAmp(fpmPupToLyot *p2l, double lambda, int sl) {
    std::complex<double> i1(0, 1);
    
    if (maskIndex == -1)
        maskIndex = sl;
    if (globalCoronagraph->get_calibration_state() & p2l->disableForCalibration)
        p2l->fpmMatAmp = exp(i1*arma::zeros<arma::mat>(size(complexMaskMatAmp.slice(0))));
    else
        p2l->fpmMatAmp = complexMaskCube.slice(maskIndex);
    
    // for the CMC we pass mask - 1 to zoomFFT
    p2l->fpmMatAmp -= 1;
}

void fpmCMCForPupToLyot::apply_babinet(fpmPupToLyot *p2l){
    p2l->fftMHatTimesfftPaddedE += p2l->paddedE;
}

void fpmCMCForPupToLyot::draw(const char *title) {
    char str[200];
    //    sprintf(str, "%s fpmMat %s", title, "CMC mask");
    //    draw_mat(fpmMat, str, "gray");
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

fpmIntHexCMCForPupToLyot::fpmIntHexCMCForPupToLyot(fpmPupToLyot *p2l, initCommandSet*& cmdBlock) {
    
    std::vector<initCommandSet*> subBlocks = cmdBlock->find_command_blocks();
    for (int i=0; i<subBlocks.size(); i++) {
        std::cout << "processing " << subBlocks[i]->commandList[0]->getCmdStr() << std::endl;
        // define the FPM
        if (!strcmp(subBlocks[i]->commandList[0]->getCmdStr(), "complexHexMaskFPM")) {
            assert(zoomFactor > 0.0);
            hexFPM = new complexHexMaskFPM(subBlocks[i], zoomFactor);
        }
        else if (!strcmp(subBlocks[i]->commandList[0]->getCmdStr(), "fpmIntHexCMCForPupToLyot")) {
            for (int c=0; c<subBlocks[i]->commandList.size(); c++) {
                set(p2l, subBlocks[i]->commandList[c]->getCmdStr(),
                    subBlocks[i]->commandList[c]->getArgStr());
            }
        }
    }
}

void fpmIntHexCMCForPupToLyot::set(fpmPupToLyot *p2l, std::string fieldName, const char *arg) {
    
    if (fieldName == "fpmIntHexCMCForPupToLyot")
        ;
    else if (fieldName == "zoomFactor") {
        // arg is a single double
        zoomFactor = atof(arg);
        std::cout << "zoomFactor = " << zoomFactor << std::endl;
    } else
        std::cout << "!!! fpmIntHexCMCForPupToLyot bad set field name: " << fieldName << std::endl;
}

void fpmIntHexCMCForPupToLyot::set_geometry(fpmPupToLyot *p2l, efield* E, double *lambda, double *lambdaFocalLength) {
    // matlab:
    // mask.dx = mask.x(2) - mask.x(1);
    // mask.dxprime = 1.098e-06; (fpmMatScale)
    // mask.resample = mask.dx / mask.dxprime;
    // mask.x = mask.x / mask.resample;
    p2l->maskGeom.set_xy(hexFPM->interpNSubPix.n_rows, hexFPM->interpNSubPix.n_cols, (E->arrayGeometry.physicalSize/2.)/(E->arrayGeometry.pixelSizeX/hexFPM->fpmScale));
    
    //    E->print("in set_geometry");
    double cFRatio = p2l->focalRatio*2*E->beamRadiusPhysical;
//    std::cout << "cFRatio set to: " << cFRatio << std::endl;
    for (int i=0; i<E->E[0][0]->n_slices; i++) {
        // lambdaFocalLength[i] = E->lambdaData[i].lambda*globalTelescope->get("primaryfRatio");
        //        std::cout << "lambda from E: " << E->lambdaData[i].lambda << std::endl;
        lambda[i] = E->lambdaData[i].lambda;
        //        std::cout << "lambda set to: " << lambda[i] << " for i = " << i << std::endl;
        
        lambdaFocalLength[i] = p2l->fRatioSign*lambda[i]*cFRatio;
    }
}

void fpmIntHexCMCForPupToLyot::set_fpmMatAmp(fpmPupToLyot *p2l, double lambda, int sl) {
    std::complex<double> i1(0, 1);
    
    if (globalCoronagraph->get_calibration_state() & p2l->disableForCalibration)
        p2l->fpmMatAmp = exp(i1*arma::zeros<arma::mat>(size(hexFPM->interpNSubPix.slice(0))));
    else
        p2l->fpmMatAmp = hexFPM->make_complex_intpolated_mask(lambda, 0.0);
    
    // for the CMC we pass mask - 1 to zoomFFT
    p2l->fpmMatAmp -= 1;
}

void fpmIntHexCMCForPupToLyot::get_optimization_data(const char *dataName, void *data) {
    hexFPM->get_optimization_data(dataName, data);
}

void fpmIntHexCMCForPupToLyot::set_optimization_data(const char *dataName, void *data) {
    hexFPM->set_optimization_data(dataName, data);
}

void fpmIntHexCMCForPupToLyot::apply_babinet(fpmPupToLyot *p2l){
    p2l->fftMHatTimesfftPaddedE += p2l->paddedE;
}

void fpmIntHexCMCForPupToLyot::draw(const char *title) {
    char str[200];
    //    sprintf(str, "%s fpmMat %s", title, "CMC mask");
    //    draw_mat(fpmMat, str, "gray");
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

fpmBinaryForPupToLyot::fpmBinaryForPupToLyot(fpmPupToLyot *p2l, initCommandSet*& cmdBlock) {
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(p2l, cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
}

void fpmBinaryForPupToLyot::set(fpmPupToLyot *p2l, std::string fieldName, const char *arg) {
    
    if (fieldName == "fpmBinaryForPupToLyot")
        ;
    else if (fieldName == "referenceLambda") {
        // arg is a single double
        referenceLambda = atof(arg);
    }
    else if (fieldName == "innerRadiusFld") {
        // arg is a single double
        innerRadiusFld = atof(arg);
    }
    else if (fieldName == "outerRadiusFld") {
        // arg is a single double
        outerRadiusFld = atof(arg);
    }
    else if (fieldName == "maskFilename") {
        // arg is the .fits filename that contains the fpm mask data
        load_mat(arg, fpmMat);
    }
    else if (fieldName == "calibFilename") {
        // arg is the .fits filename that contains the calibration mask data
        load_mat(arg, calibMat);
    }
    else
        std::cout << "!!! fpmBinaryForPupToLyot bad set field name: " << fieldName << std::endl;
}


void fpmBinaryForPupToLyot::set_geometry(fpmPupToLyot *p2l, efield* E, double *lambda, double *lambdaFocalLength) {
    // matlab:
    // mask.flD = system.params.fpm.flD;%system.params.beam.f*lambda/system.optics.elem(i-1).Dx;
    // mask.gridsize = system.params.fpm.FOVflD*mask.flD;
    // mask.dx = mask.gridsize / system.optics.elem(i).N;
    // mask.dy = mask.dx;
    // mask.x = system.optics.elem(i).x;
    // mask.y = system.optics.elem(i).y;
    // mask.xx = system.optics.elem(i).xx;
    // mask.yy = system.optics.elem(i).yy;
    
    // params.fpm.flD = params.fpm.f*params.fpm.lambdaRef/params.primary.D; % define focal plane mask units of flD
    double flD = globalTelescope->get("primaryfRatio")*referenceLambda;
    double fovFld = 2*outerRadiusFld;
    p2l->maskGeom.set_xy_m1(fpmMat.n_rows, fpmMat.n_cols, outerRadiusFld*flD);
    p2l->maskGeom.print("puplToLyot mask geometry");
    
    //    E->print("in set_geometry");
    for (int i=0; i<E->E[0][0]->n_slices; i++) {
        lambda[i] = E->lambdaData[i].lambda;
        lambdaFocalLength[i] = p2l->fRatioSign*lambda[i]*globalTelescope->get("primaryfLength");
    }
}

void fpmBinaryForPupToLyot::set_fpmMatAmp(fpmPupToLyot *p2l, double lambda, int sl) {
    std::complex<double> r1(1, 0);
    
    if (globalCoronagraph->get_calibration_state() & p2l->disableForCalibration)
        p2l->fpmMatAmp = r1*calibMat;
    else
        p2l->fpmMatAmp = r1*fpmMat;
    
    // matlab: mask.M_hat = zoomFFT_realunits(mask.x, mask.y, mask.M, pupil_ext.x, pupil_ext.y, pupil.f, lambda);
}

void fpmBinaryForPupToLyot::draw(const char *title) {
    char str[200];
    sprintf(str, "%s fpmMat %s", title, "opaque mask");
    draw_mat(fpmMat, str, "gray");
}

