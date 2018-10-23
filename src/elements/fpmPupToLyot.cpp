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
    static int p2lCount = 0;
    
    assert(E->beamRadiusPhysical > 0);

    std::complex<double> i1(0, 1);
    std::complex<double> r1(1, 0);
    int nRowsE = E->E[0][0]->n_rows;
    int nColsE = E->E[0][0]->n_cols;

    pre_execute(E, prev, next, time);
    
//    E->arrayGeometry.print("E entering pupil2Lyot");
    
    omp_set_nested(1);
    if (!maskIsInited) {
        // this needs to be computed whenever the mask changes
//        std::cout << "------------- initing fpmPupToLyot ---------------" << std::endl;
        fftMHatTimesfftPaddedEScale.set_size(E->E[0][0]->n_slices);

        double *lambda = new double[E->E[0][0]->n_slices];
        if (lambdaFocalLength == NULL)
            lambdaFocalLength = new double[E->E[0][0]->n_slices];
        fflush(stdout);
        // create space for a zero-padded E array, with a 2x padding
        paddedE = arma::zeros<arma::cx_cube>(2*nRowsE, 2*nColsE, E->E[0][0]->n_slices);
        for (int sl=0; sl<E->E[0][0]->n_slices; sl++) {
            paddedGeom.set_xy(paddedE.slice(sl), E->arrayGeometry.pixelSizeX);
        }
        mask->set_geometry(this, E, lambda, lambdaFocalLength);
        fpmMatAmpCalib = arma::zeros<arma::cx_cube>(maskGeom.pixelX.n_elem, maskGeom.pixelY.n_elem, E->E[0][0]->n_slices);
        fpmMatAmp = arma::zeros<arma::cx_cube>(maskGeom.pixelX.n_elem, maskGeom.pixelY.n_elem, E->E[0][0]->n_slices);
        fftMaskHatCalib = arma::zeros<arma::cx_cube>(2*nRowsE, 2*nColsE, E->E[0][0]->n_slices);
        fftMaskHat = arma::zeros<arma::cx_cube>(2*nRowsE, 2*nColsE, E->E[0][0]->n_slices);

        if (propZoomFft == NULL)
            propZoomFft = new zoomFft[omp_get_max_threads()];
        fflush(stdout);
        if (paddedEFft == NULL)
            paddedEFft = new fft[omp_get_max_threads()];
        fflush(stdout);
        if (myIfft == NULL)
            myIfft = new ifft[omp_get_max_threads()];
        fflush(stdout);
        fftw_plan_with_nthreads(omp_get_max_threads());
        for (int t=0; t<omp_get_max_threads(); t++) {
            paddedEFft[t].init(paddedE.slice(0));
            myIfft[t].init(paddedE.slice(0));
        }
        fftw_plan_with_nthreads(1);
//        maskGeom.print("pupil2Lyot mask:");
//        paddedGeom.print("pupil2Lyot padded:");
        #pragma omp parallel for
        for (int sl=0; sl<E->E[0][0]->n_slices; sl++) {
            arma::cx_mat maskHat;
            fft maskFft;
            
            propZoomFft[omp_get_thread_num()].init(maskGeom, paddedGeom, lambdaFocalLength, E->E[0][0]->n_slices);
//            std::cout << "---- set_fpmMatAmp" << std::endl;
            mask->set_fpmMatAmp(this, lambda[sl], sl);
            // matlab: mask.M_hat = zoomFFT_realunits(mask.x, mask.y, mask.M - 1, pupil_ext.x, pupil_ext.y, pupil.f, lambda);
            // matlab: FFT_M_hat = fft2((fftshift(mask.M_hat)));
            if (!fftMaskHatCalibComputed) {
                // don't need to execute whenever the mask changes
                maskHat = zoomFftSign*propZoomFft[omp_get_thread_num()].execute(fpmMatAmpCalib.slice(sl), sl);
                fftMaskHatCalib.slice(sl) = fft_shift(maskHat);
                maskFft.execute(fftMaskHatCalib.slice(sl));
                fftMaskHatCalibComputed = true;
            }
        
            maskHat = zoomFftSign*propZoomFft[omp_get_thread_num()].execute(fpmMatAmp.slice(sl), sl);
            fftMaskHat.slice(sl) = fft_shift(maskHat);
            maskFft.execute(fftMaskHat.slice(sl));
            
            fftMHatTimesfftPaddedEScale[sl] = pow(paddedGeom.pixelSizeX, 2)/(-i1*fRatioSign*lambdaFocalLength[sl]);
        }
//        save_cube("maskHat", fftMaskHat);
//        save_cube("maskHatCalib", fftMaskHatCalib);
//        save_cube("mask", fpmMatAmp);
//        save_cube("maskCalib", fpmMatAmpCalib);
        
//        load_cube("donutPiaaLuvoir/matlabMaskHat", fftMaskHat);
//        fftPaddedE.zeros(2*nRowsE, 2*nColsE);
//        fftMHatTimesfftPaddedE.zeros(size(paddedE.slice(0)));

        maskIsInited = true;
    }

    for (int sl=0; sl<E->E[0][0]->n_slices; sl++) {
//        arma::wall_clock timer1;
//        timer1.tic();
//        fftPaddedE.zeros();
//        fftMHatTimesfftPaddedE.zeros();
//        std::cout << "zeroing padded E:" << timer1.toc() << " s" << std::endl;

//        timer1.tic();
        for (int s=0; s<E->E.size(); s++) {
            for (int p=0; p<E->E[s].size(); p++) {
//                arma::wall_clock timer2;
                // embed E in the zero-padded cube
//                timer2.tic();
                paddedE(arma::span(nRowsE/2, 3*nRowsE/2-1), arma::span(nColsE/2, 3*nColsE/2-1), arma::span(sl, sl)) = E->E[s][p]->slice(sl);
//                std::cout << "setting paddedE:" << timer2.toc() << " s" << std::endl;

//                timer2.tic();
                // matlab: FFT_pupil = fft2(pupil_ext.E);
                fftPaddedE = paddedE.slice(sl);
                paddedEFft[omp_get_thread_num()].execute(fftPaddedE);
//                std::cout << "computing paddedEFft:" << timer2.toc() << " s" << std::endl;

                // matlab: E = ifft2(FFT_M_hat.*FFT_pupil)/(-1i*lambda*pupil.f)*(pupil.x(2) - pupil.x(1))^2 + pupil_ext.E;
//                timer2.tic();
                if (globalCoronagraph->get_calibration_state() & disableForCalibration)
                    fftMHatTimesfftPaddedE = fftMaskHatCalib.slice(sl)%fftPaddedE;
                else
                    fftMHatTimesfftPaddedE = fftMaskHat.slice(sl)%fftPaddedE;
//                std::cout << "setting up fftMHatTimesfftPaddedE:" << timer2.toc() << " s" << std::endl;
//                timer2.tic();
                myIfft[omp_get_thread_num()].execute(fftMHatTimesfftPaddedE);
//                std::cout << "computing fftMHatTimesfftPaddedE:" << timer2.toc() << " s" << std::endl;

//                timer2.tic();
                fftMHatTimesfftPaddedE = fftMHatTimesfftPaddedE*fftMHatTimesfftPaddedEScale[sl];
//                std::cout << "fftMHatTimesfftPaddedEScale:" << fftMHatTimesfftPaddedEScale[sl] << ", formula:" << pow(paddedGeom.pixelSizeX, 2)/(-i1*fRatioSign*lambdaFocalLength[sl]) << std::endl;
//                fftMHatTimesfftPaddedE = fftMHatTimesfftPaddedE*pow(paddedGeom.pixelSizeX, 2)/(-i1*fRatioSign*lambdaFocalLength[sl]);
//                std::cout << "fixing fftMHatTimesfftPaddedE:" << timer2.toc() << " s" << std::endl;

                if (mask->doBabinet)
                    fftMHatTimesfftPaddedE += paddedE.slice(sl);
                
                // extract the central result for the final answer
//                timer2.tic();
                E->E[s][p]->slice(sl) = fftMHatTimesfftPaddedE(arma::span(nRowsE/2, 3*nRowsE/2-1), arma::span(nColsE/2, 3*nColsE/2-1));
//                std::cout << "extracting E:" << timer2.toc() << " s" << std::endl;
            }
        }
//        std::cout << "pupil2lyot execute:" << timer1.toc() << " s" << std::endl;
    }
    post_execute(E, prev, next, time);
//    E->print("exiting fpmPupToLyot::execute: ");
    p2lCount++;
    
    return E;
}

void fpmPupToLyot::init(initCommandSet*& cmdBlock) {
//        std::cout << "initing a fpmPupToLyot" << std::endl;
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
//    std::cout << "initing a fpmPupToLyot command block" << std::endl;

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
    if (!strcmp(dataName, "reinit"))
        maskIsInited = false;
    else
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
    set_babinet(true); // turn on babinet for this mask
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
    } else if (fieldName == "maskFilenameReIm") {
        // arg is two filenames separated by a comma, each giving
        // the .fits amplitude and phase filename that contains the complexMask definition
        const char *cPtr = strchr(arg, ',');
        int strReLen = cPtr - arg;
        int strImLen = strlen(arg) - strReLen + 1;
        std::cout << "strReLen " << strReLen << ", strImLen " << strImLen << std::endl;
        char *strRe = new char[strReLen + 1];
        char *strIm = new char[strImLen + 1];
        strncpy(strRe, arg, strReLen);
        strRe[strReLen] = '\0';
        strncpy(strIm, cPtr + 1, strImLen);
        strIm[strImLen] = '\0';
        std::cout << "loading " << strRe << " and " << strIm << std::endl;
        initMaskReIm(strRe, strIm);
                draw("initial FPM");
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

void fpmCMCForPupToLyot::initMaskReIm(const char *filenameRe, const char *filenameIm) {
    
    load_cube(filenameRe, complexMaskMatRe);
    load_cube(filenameIm, complexMaskMatIm);
    complexMaskCube.set_size(size(complexMaskMatRe));
    complexMaskMatAmp.set_size(size(complexMaskMatRe));
    complexMaskMatPh.set_size(size(complexMaskMatRe));
    for (int sl=0; sl<complexMaskMatRe.n_slices; ++sl) {
        complexMaskCube.slice(sl).set_real(complexMaskMatRe.slice(sl));
        complexMaskCube.slice(sl).set_imag(complexMaskMatIm.slice(sl));
        complexMaskMatAmp.slice(sl) = abs(complexMaskCube.slice(sl));
        complexMaskMatPh.slice(sl) = arg(complexMaskCube.slice(sl));
    }
//    save_mat("donutfpm.fits", complexMaskCube.slice(0));
}

void fpmCMCForPupToLyot::set_geometry(fpmPupToLyot *p2l, efield* E, double *lambda, double *lambdaFocalLength) {
    // matlab:
    // mask.dx = mask.x(2) - mask.x(1);
    // mask.dxprime = 1.098e-06; (fpmMatScale)
    // mask.resample = mask.dx / mask.dxprime;
    // mask.x = mask.x / mask.resample;
//    p2l->maskGeom.set_xy_m1(complexMaskCube.n_rows, complexMaskCube.n_cols, (E->arrayGeometry.physicalSize/2.)/(E->arrayGeometry.pixelSizeX/fpmMatScale));
    p2l->maskGeom.set_xy(complexMaskCube, fpmMatScale);

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
//    if (globalCoronagraph->get_calibration_state() & p2l->disableForCalibration)
        p2l->fpmMatAmpCalib.slice(sl) = exp(i1*arma::zeros<arma::mat>(size(complexMaskMatAmp.slice(0))));
//    else
        p2l->fpmMatAmp.slice(sl) = complexMaskCube.slice(maskIndex);
    
    // for the CMC we pass mask - 1 to zoomFFT
    if (doBabinet) {
        p2l->fpmMatAmpCalib.slice(sl) -= 1;
        p2l->fpmMatAmp.slice(sl) -= 1;
    }
}

//void fpmCMCForPupToLyot::apply_babinet(fpmPupToLyot *p2l, int sl){
//    if (doBabinet)
//        p2l->fftMHatTimesfftPaddedE += p2l->paddedE.slice(sl);
//}

void fpmCMCForPupToLyot::draw(const char *title) {
    char str[200];
    //    sprintf(str, "%s fpmMat %s", title, "CMC mask");
    //    draw_mat(fpmMat, str, "gray");
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

fpmIntHexCMCForPupToLyot::fpmIntHexCMCForPupToLyot(fpmPupToLyot *p2l, initCommandSet*& cmdBlock) {
    
//    std::cout << "in fpmIntHexCMCForPupToLyot constructor" << std::endl;
    std::vector<initCommandSet*> subBlocks = cmdBlock->find_command_blocks();
    for (int i=0; i<subBlocks.size(); i++) {
//        std::cout << "processing " << subBlocks[i]->commandList[0]->getCmdStr() << std::endl;
        // define the FPM
        if (!strcmp(subBlocks[i]->commandList[0]->getCmdStr(), "complexHexMaskFPM")) {
//            std::cout << "zoomFactor = " << zoomFactor << std::endl;
            assert(zoomFactor > 0.0);
            hexFPM = new complexHexMaskFPM(subBlocks[i], zoomFactor);
        }
        else if (!strcmp(subBlocks[i]->commandList[0]->getCmdStr(), "fpmIntHexCMCForPupToLyot")) {
            for (int c=0; c<subBlocks[i]->commandList.size(); c++) {
                set(p2l, subBlocks[i]->commandList[c]->getCmdStr(),
                    subBlocks[i]->commandList[c]->getArgStr());
            }
        }
//        std::cout << "zoomFactor = " << zoomFactor << std::endl;
    }
    set_babinet(true); // turn on babinet for this mask
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
//    p2l->maskGeom.set_xy_m1(hexFPM->interpNSubPix.n_rows, hexFPM->interpNSubPix.n_cols, (E->arrayGeometry.physicalSize/2.)/(E->arrayGeometry.pixelSizeX/hexFPM->pixelScale));
    p2l->maskGeom.set_xy(hexFPM->interpNSubPix, hexFPM->pixelScale);
//    p2l->maskGeom.print("FPM mask geometry");
    //    E->print("in set_geometry");
    double cFRatio = p2l->focalRatio*2*E->beamRadiusPhysical; // focal length at this point
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
    
    p2l->fpmMatAmpCalib.slice(sl) = arma::fliplr(exp(i1*arma::zeros<arma::mat>(size(hexFPM->interpNSubPix.slice(0)))));
    p2l->fpmMatAmp.slice(sl) = arma::fliplr(hexFPM->make_complex_intpolated_mask(lambda, 0.0));
//    arma::cx_mat matlabMask;
//    load_mat("donutPiaaLuvoir/matlabMask", matlabMask);
//    p2l->fpmMatAmp.slice(sl) = matlabMask;
    
    // for the CMC we pass mask - 1 to zoomFFT
    if (doBabinet) {
        p2l->fpmMatAmpCalib.slice(sl) -= 1;
        p2l->fpmMatAmp.slice(sl) -= 1;
    }
}

void fpmIntHexCMCForPupToLyot::get_optimization_data(const char *dataName, void *data) {
    hexFPM->get_optimization_data(dataName, data);
}

void fpmIntHexCMCForPupToLyot::set_optimization_data(const char *dataName, void *data) {
    if (!strcmp(dataName, "setBabinet"))
        doBabinet = *(bool *)data;
    else
        hexFPM->set_optimization_data(dataName, data);
}

//void fpmIntHexCMCForPupToLyot::apply_babinet(fpmPupToLyot *p2l, int sl){
//    if (doBabinet)
//        p2l->fftMHatTimesfftPaddedE += p2l->paddedE.slice(sl);
//}

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
//    double flD = globalTelescope->get("primaryfRatio")*referenceLambda;
    double flD = globalTelescope->compute_loD(referenceLambda);
    double fovFld = 2*outerRadiusFld;
//    p2l->maskGeom.set_xy_offset_m1(fpmMat.n_cols, fpmMat.n_rows, 2*outerRadiusFld*flD/fpmMat.n_rows);
//    p2l->maskGeom.set_xy_m1(fpmMat.n_rows, fpmMat.n_cols, outerRadiusFld*flD);
    p2l->maskGeom.set_xy(fpmMat, 2*outerRadiusFld*flD/fpmMat.n_rows);
//    p2l->maskGeom.print("puplToLyot mask geometry");
    
    //    E->print("in set_geometry");
    for (int i=0; i<E->E[0][0]->n_slices; i++) {
        lambda[i] = E->lambdaData[i].lambda;
        lambdaFocalLength[i] = p2l->fRatioSign*lambda[i]*globalTelescope->get("primaryfLength");
    }
}

void fpmBinaryForPupToLyot::set_fpmMatAmp(fpmPupToLyot *p2l, double lambda, int sl) {
    std::complex<double> r1(1, 0);
    
//    if (globalCoronagraph->get_calibration_state() & p2l->disableForCalibration)
        p2l->fpmMatAmpCalib.slice(sl) = r1*calibMat;
//    else
        p2l->fpmMatAmp.slice(sl) = r1*fpmMat;
    
    // matlab: mask.M_hat = zoomFFT_realunits(mask.x, mask.y, mask.M, pupil_ext.x, pupil_ext.y, pupil.f, lambda);
}

void fpmBinaryForPupToLyot::draw(const char *title) {
    char str[200];
    sprintf(str, "%s fpmMat %s", title, "opaque mask");
    draw_mat(fpmMat, str, "gray");
}

