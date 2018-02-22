//
//  csim.cpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//


#include <iostream>
#include <string>
#include <cstring>
#include "test.hpp"
#include "efield.hpp"
#include "coronagraph.hpp"
#include "telescope.hpp"
#include "csim_lib.hpp"

int main(int argn, char **argv) {
    std::complex<double> r1(1, 0); // for conversion from real to complex
    std::complex<double> i1(0, 1); // for conversion from real to complex
    
    if (argn == 2) {
        if (!strcmp(argv[1], "?")) {
            std::cout << "to display details, pass one of the following:" << std::endl;
            std::cout << "\tfits" << std::endl;
            std::cout << "\tfftShift" << std::endl;
            std::cout << "\tfresnelPropagationAS" << std::endl;
            std::cout << "\tczt" << std::endl;
            std::cout << "\tallUpPIAACMC" << std::endl;
            
            
            return 0;
        }
    }

    //////////////////////////////////////////
    // test fits io on mat
    //////////////////////////////////////////
    {
        arma::mat m = arma::randu(4,6);
        arma::mat mr = arma::randu(4,6);
        save_mat("test_mat.fits", m);
        load_mat("test_mat.fits", mr);
        double saveMatError = max(max(m - mr));
        if (saveMatError > 1e-16) {
            std::cout << "!!!!! save mat result does not match test data" << std::endl;
            std::cout << "max error: " << saveMatError << std::endl;
        } else
            std::cout << "=================> save mat test pass" << std::endl;
        if (argn == 2) {
            if (!strcmp(argv[1], "fits")) {
                m.print("m:");
                mr.print("mr:");
            }
        }
    }
    
    //////////////////////////////////////////
    // test fits io on cx_mat
    //////////////////////////////////////////
    {
        arma::cx_mat m = arma::randu<arma::cx_mat>(4,6);
        arma::cx_mat mr = arma::randu<arma::cx_mat>(4,6);
        save_mat("test_cx_mat.fits", m);
        load_mat("test_cx_mat.fits", mr);
        double saveMatError = max(max(abs(m - mr)));
        if (saveMatError > 1e-16) {
            std::cout << "!!!!! save cx_mat result does not match test data" << std::endl;
            std::cout << "max error: " << saveMatError << std::endl;
        } else
            std::cout << "=================> save cx_mat test pass" << std::endl;
        if (argn == 2) {
            if (!strcmp(argv[1], "fits")) {
                m.print("m:");
                mr.print("mr:");
            }
        }
    }
    
    //////////////////////////////////////////
    // test fftShift
    //////////////////////////////////////////
    {
        arma::mat fftShiftInputMat;
        arma::mat fftShiftRequiredOutputMat;
        load_mat("testData/fftShift/fftShift_in.fits", fftShiftInputMat);
        load_mat("testData/fftShift/fftShift_out.fits", fftShiftRequiredOutputMat);
        arma::cx_cube testIn(fftShiftInputMat.n_rows, fftShiftInputMat.n_cols, 1);
        testIn.slice(0) = r1*fftShiftInputMat;
        arma::cx_cube testOut = fft_shift(testIn);
        double fftShiftError = max(max(abs(testOut.slice(0) - fftShiftRequiredOutputMat)));
        if (fftShiftError > 1e-16) {
            std::cout << "!!!!! fftShift result does not match test data" << std::endl;
            std::cout << "max error: " << fftShiftError << std::endl;
        } else
            std::cout << "=================> fftShift test pass" << std::endl;
        
        if (argn == 2) {
            if (!strcmp(argv[1], "fftShift")) {
                draw_mat(fftShiftInputMat, "fft shift input");
                draw_mat(fftShiftRequiredOutputMat, "fft shift required output");
                draw_mat(abs(testOut.slice(0)), "fft shift test output");
            }
        }
    }
    
    //////////////////////////////////////////
    // test fresnelPropagationAS
    //////////////////////////////////////////
    {
        arma::mat fresnelPropagationASInputMat;
        arma::mat fresnelPropagationASRequiredOutputMat;
        load_mat("testData/fresnelPropagateAS/fresnelPropagateAS_in.fits", fresnelPropagationASInputMat);
        load_mat("testData/fresnelPropagateAS/fresnelPropagateAS_out.fits", fresnelPropagationASRequiredOutputMat);
        arma::cx_cube testIn(fresnelPropagationASInputMat.n_rows, fresnelPropagationASInputMat.n_cols, 1);
        testIn.slice(0) = r1*fresnelPropagationASInputMat;
        FILE *fid = fopen("testData/fresnelPropagateAS/fresnel_test_params.txt", "r");
        double lambda;
        double apRad;
        double z;
        fscanf(fid, "%lf %lf %lf", &lambda, &apRad, &z);
        fclose(fid);
        std::cout << "read fresnel test params lambda = " << lambda << " apRad = " << apRad << " z = " << z << std::endl;
        
        testIn.set_size(fresnelPropagationASInputMat.n_rows, fresnelPropagationASInputMat.n_cols, 1);
        testIn.slice(0) = r1*fresnelPropagationASInputMat;
        fresnelPropagateAS frenelProp;
        arma::cx_cube testOut = testIn;
        frenelProp.execute(testOut, &lambda, apRad, z);
        arma::vec fresnelPropagateASError = vectorise(abs(abs(testOut.slice(0)) - fresnelPropagationASRequiredOutputMat));
        if (max(fresnelPropagateASError/vectorise(fresnelPropagationASRequiredOutputMat)) > 1e-16) {
            std::cout << "!!!!! fresnelPropagateAS result does not match test data" << std::endl;
            std::cout << "max relative error: " << max(fresnelPropagateASError/vectorise(fresnelPropagationASRequiredOutputMat)) << std::endl;
            std::cout << "median relative error: " << median(fresnelPropagateASError/vectorise(fresnelPropagationASRequiredOutputMat)) << std::endl;
            std::cout << "mean relative error: " << mean(fresnelPropagateASError/vectorise(fresnelPropagationASRequiredOutputMat)) << std::endl;
            std::cout << "max absolute error: " << max(fresnelPropagateASError) << std::endl;
            std::cout << "median absolute error: " << median(fresnelPropagateASError) << std::endl;
            std::cout << "mean absolute error: " << mean(fresnelPropagateASError) << std::endl;
        } else
            std::cout << "=================> fresnelPropagateAS test pass" << std::endl;
        
        
        if (argn == 2) {
            if (!strcmp(argv[1], "fresnelPropagationAS")) {
                draw_mat(fresnelPropagationASInputMat, "fresnelPropagationAS input");
                draw_mat(log10(fresnelPropagationASRequiredOutputMat), " log fresnelPropagationAS required output");
                draw_mat(log10(abs(testOut.slice(0))), "log fresnelPropagationAS test output");
                draw_mat(log10(fresnelPropagationASRequiredOutputMat - abs(testOut.slice(0))),
                         "log fresnelPropagationAS difference");
            }
        }
    }

    //////////////////////////////////////////
    // test czt
    //////////////////////////////////////////
    {
        arma::mat cztInputReMat;
        arma::mat cztInputImMat;
        arma::mat cztRequiredOutputReMat;
        arma::mat cztRequiredOutputImMat;
        load_mat("testData/czt/czt_in_re.fits", cztInputReMat);
        load_mat("testData/czt/czt_in_im.fits", cztInputImMat);
        load_mat("testData/czt/czt_out_re.fits", cztRequiredOutputReMat);
        load_mat("testData/czt/czt_out_im.fits", cztRequiredOutputImMat);
        fflush(stdout);
        arma::cx_mat testIn(cztInputReMat.n_rows, cztInputReMat.n_cols);
        testIn = r1*cztInputReMat + i1*cztInputImMat;
        arma::cx_mat cztRequiredOutputMat(cztRequiredOutputReMat.n_rows, cztRequiredOutputReMat.n_cols);
        
        FILE *fid = fopen("testData/czt/czt_test_params.txt", "r");
        int k;
        double w;
        double a;
        fscanf(fid, "%d %lf %lf", &k, &w, &a);
        fclose(fid);
        printf("read czt test params k = %d, w = %0.20f, a = %0.20f\n", k, w, a);

        czt testCzt;
        arma::cx_mat testOut = testCzt.execute(testIn, k, w, a);
        std::cout << "testOut = " << abs(testOut(arma::span(0,4),arma::span(0,4))) << std::endl;
        save_mat("cztTestOut.fits", testOut);
/*
        double wr, wc;
        double ar, ac;
        fscanf(fid, "%d %lf %lf %lf %lf", &k, &wr, &wc, &ar, &ac);
        fclose(fid);
        std::complex<double> w(wr, wc);
        std::complex<double> a(ar, ac);
        std::cout << "read czt test params k = " << k << " w = " << w << " a = " << a << std::endl;
        testIn.set_size(cztInputMat.n_rows, cztInputMat.n_cols);
        testIn = r1*cztInputMat;
        czt testCzt;
        arma::cx_mat testOut = testCzt.execute(testIn, k, w, a);
 */
        arma::vec cztError = vectorise(abs(abs(testOut) - abs(cztRequiredOutputMat)));
        arma::vec cztReError = vectorise(abs(real(testOut) - cztRequiredOutputReMat));
        arma::vec cztImError = vectorise(abs(imag(testOut) - cztRequiredOutputImMat));
        arma::vec cztRequiredAbsValue = vectorise(abs(cztRequiredOutputMat));

        if (max(cztError/vectorise(abs(cztRequiredOutputMat))) > 2e-15) {
            std::cout << "!!!!! czt result does not match test data" << std::endl;
            std::cout << "max relative error: " << max(cztError/cztRequiredAbsValue) << std::endl;
            std::cout << "median relative error: " << median(cztError/cztRequiredAbsValue) << std::endl;
            std::cout << "mean relative error: " << mean(cztError/cztRequiredAbsValue) << std::endl;
            std::cout << "max absolute error: " << max(cztError) << std::endl;
            std::cout << "median absolute error: " << median(cztError) << std::endl;
            std::cout << "mean absolute error: " << mean(cztError) << std::endl;
            std::cout << "==============" << std::endl;
            std::cout << "Real part: max absolute error: " << max(cztReError) << std::endl;
            std::cout << "Real part: median absolute error: " << median(cztReError) << std::endl;
            std::cout << "Real part: mean absolute error: " << mean(cztReError) << std::endl;
            std::cout << "==============" << std::endl;
            std::cout << "Imaginary part: max absolute error: " << max(cztImError) << std::endl;
            std::cout << "Imaginary part: median absolute error: " << median(cztImError) << std::endl;
            std::cout << "Imaginary part: mean absolute error: " << mean(cztImError) << std::endl;
        } else
            std::cout << "=================> czt test pass" << std::endl;
        
        if (argn == 2) {
            arma::cx_mat dMat = testOut;
            if (!strcmp(argv[1], "czt")) {
                draw_mat(abs(testIn), "czt input");
                draw_mat(log10(abs(cztRequiredOutputMat)), "log czt required output");
                draw_mat(log10(abs(dMat)), "log czt test output");
                draw_mat(log10((abs(cztRequiredOutputMat) - abs(testOut))/abs(cztRequiredOutputMat)), "log czt difference");
                
                draw_mat(log10(cztRequiredOutputReMat), "log czt required real part output");
                draw_mat(log10(real(testOut)), "log czt test real part output");
                draw_mat(log10(cztRequiredOutputImMat), "log czt required imaginary part output");
                draw_mat(log10(imag(testOut)), "log czt test real part output");
                
            }
        }
    }
    
    
    //////////////////////////////////////////
    // all-up test PIAACMC coronagraph
    //////////////////////////////////////////
    {
        arma::mat allUpRequiredOutputMat;
        load_mat("testData/allUpPIAACMC/all_up_test_out.fits", allUpRequiredOutputMat);
        
        std::system("src/csim testData/allUpPIAACMC/PIAACMC_all_up_test_input.txt");
        
        arma::cube allUpRunOutputCube;
        load_cube("testData/allUpPIAACMC/E_s0_p0_science_CCD_out_amp.fits", allUpRunOutputCube);
        
        arma::mat allUpRunOutputMat = allUpRunOutputCube.slice(0);
        arma::vec allUpError = vectorise(abs(allUpRunOutputMat - allUpRequiredOutputMat));
        if (max(allUpError > 2e-15)) {
            std::cout << "!!!!! all up PIAACMC result does not match test data" << std::endl;
            std::cout << "max relative error: " << max(allUpError/vectorise(allUpRequiredOutputMat+1e-16)) << std::endl;
            std::cout << "median relative error: " << median(allUpError/vectorise(allUpRequiredOutputMat+1e-16)) << std::endl;
            std::cout << "mean relative error: " << mean(allUpError/vectorise(allUpRequiredOutputMat+1e-16)) << std::endl;
            std::cout << "max absolute error: " << max(allUpError) << std::endl;
            std::cout << "median absolute error: " << median(allUpError) << std::endl;
            std::cout << "mean absolute error: " << mean(allUpError) << std::endl;
        } else
            std::cout << "=================> all up PIAACMC test pass" << std::endl;
        
        if (argn == 2) {
            if (!strcmp(argv[1], "allUpPIAACMC")) {
                draw_mat(log10(allUpRequiredOutputMat), "log allUpPIAACMC required output");
                draw_mat(log10(allUpRunOutputMat), "log allUpPIAACMC test output");
                draw_mat(log10(abs(allUpRequiredOutputMat - allUpRunOutputMat)), "log allUpPIAACMC difference");
            }
        }
    }
    return 0;

/*
    // test fitsio
    void *data = NULL;
    int nDims = 0;
    long *dims = NULL;
    int dataType = 0;
    read_fits_array("pup_1024.fits", &data, nDims, &dims, dataType);
    std::cout << "data = " << data << std::endl;
    std::cout << "read fits array : nDims = " << nDims << std::endl;
    std::cout << "dims = " << dims << std::endl;
    for (int i=0; i<nDims; i++) {
        std::cout << "axis " << i << ": " << dims[i] << std::endl;
    }
    
    write_fits_array("!test.fits", data, nDims, dims, dataType);
*/    

/*
 // test gnuplot
    void *data = NULL;
    int nDims = 0;
    long *dims = NULL;
    int dataType = 0;
    read_fits_array("pup_1024.fits", &data, nDims, &dims, dataType);
    std::cout << "data = " << data << std::endl;
    std::cout << "read fits array : nDims = " << nDims << std::endl;
    std::cout << "dims = " << dims << std::endl;
    for (int i=0; i<nDims; i++) {
    std::cout << "axis " << i << ": " << dims[i] << std::endl;
    }
    
    FILE *gnuDataFile = fopen("gnuData.dat", "w");
    for (int j = 0; j < dims[1]; j++) {
        for (int i = 0; i < dims[0]; i++) {
            fprintf(gnuDataFile, "%lf ", *((float *)data + j*dims[0] + i));
        }
        fprintf(gnuDataFile, "\n");
    }
    fclose(gnuDataFile);
    
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    fprintf(gnuplotPipe, "set title \"test\" \n");
    fprintf(gnuplotPipe, "set size ratio 1 \n");
//    fprintf(gnuplotPipe, "set palette rgb 7,5,15 \n"); // default
//    fprintf(gnuplotPipe, "set palette rgb 21,22,23 \n"); // hot
    fprintf(gnuplotPipe, "set palette rgb 33,13,10 \n"); // rainbow
    fprintf(gnuplotPipe, "plot 'gnuData.dat' matrix with image\n");
    pclose(gnuplotPipe);
*/
}
