//
//  csim_fft.cpp
//  csim
//
//  Created by steve on 4/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include <iostream>
#include <assert.h>

#include "csim_lib.hpp"
#include <math.h>
#include <omp.h>

static int fftPlanType = FFTW_ESTIMATE;
//static int fftPlanType = FFTW_PATIENT;
//static int fftPlanType = FFTW_EXHAUSTIVE;

void init_fft_lib(bool computeFftWisdon) {
    fftw_import_wisdom_from_filename("fftw_wisdom.dat");
    
    if (computeFftWisdon) {
        fftPlanType = FFTW_EXHAUSTIVE;
        std::cout << "computing FFTW widsom" << std::endl;
    } else
        fftPlanType = FFTW_ESTIMATE;
    
    fftw_init_threads();
    fftw_plan_with_nthreads(1);
}

void save_fft_wisdom(void) {
    fftw_export_wisdom_to_filename("fftw_wisdom.dat");
}

void fft::init(arma::cx_cube& in, int length) {
    if (length > 0 & in.n_rows != length) {
        in.reshape(length, in.n_cols, in.n_slices);
    }
    
    if (fftPlan.nRows != in.n_rows || fftPlan.nCols != in.n_cols || fftPlan.nSlices != in.n_slices) {
#pragma omp critical (fftPlanning)  // let's be safe
        {
            if (fftPlan.nRows != 0)
                fftw_destroy_plan(fftPlan.plan);
            fftPlan.plan = fftw_plan_dft_3d(in.n_cols, in.n_rows, in.n_slices, (double(*)[2])&in(0,0,0), (double(*)[2])&in(0,0,0), FFTW_FORWARD, FFTW_WISDOM_ONLY);
            if (fftPlan.plan == NULL) {
                std::cout << "no wisdom, estimating" << std::endl;
                arma::cx_cube tmp = in; // fftw_plan_dft_3d may mess with in
                fftPlan.plan = fftw_plan_dft_3d(in.n_cols, in.n_rows, in.n_slices, (double(*)[2])&in(0,0,0), (double(*)[2])&in(0,0,0), FFTW_FORWARD, fftPlanType);
                in = tmp;
            }
            fftPlan.nRows = in.n_rows;
            fftPlan.nCols = in.n_cols;
            fftPlan.nSlices = in.n_slices;
        }
    }
}

void fft::execute(arma::cx_cube& in, int length) {
    init(in, length);
    fftw_execute_dft(fftPlan.plan, (double(*)[2])&in(0,0,0), (double(*)[2])&in(0,0,0));
}

void ifft::init(arma::cx_cube& in, int length) {
    if (length > 0 & in.n_rows != length) {
        in.reshape(length, in.n_cols, in.n_slices);
    }
    
    if (fftPlan.nRows != in.n_rows || fftPlan.nCols != in.n_cols || fftPlan.nSlices != in.n_slices) {
#pragma omp critical (fftPlanning)  // let's be safe
        {
            if (fftPlan.nRows != 0)
                fftw_destroy_plan(fftPlan.plan);
            fftPlan.plan = fftw_plan_dft_3d(in.n_cols, in.n_rows, in.n_slices, (double(*)[2])&in(0,0,0), (double(*)[2])&in(0,0,0), FFTW_BACKWARD, FFTW_WISDOM_ONLY);
            if (fftPlan.plan == NULL) {
                std::cout << "no wisdom, estimating" << std::endl;
                arma::cx_cube tmp = in; // fftw_plan_dft_3d may mess with in
                fftPlan.plan = fftw_plan_dft_3d(in.n_cols, in.n_rows, in.n_slices, (double(*)[2])&in(0,0,0), (double(*)[2])&in(0,0,0), FFTW_BACKWARD, fftPlanType);
                in = tmp;
            }
            fftPlan.nRows = in.n_rows;
            fftPlan.nCols = in.n_cols;
            fftPlan.nSlices = in.n_slices;
        }
    }
}

void ifft::execute(arma::cx_cube& in, int length) {
    init(in, length);
    fftw_execute_dft(fftPlan.plan, (double(*)[2])&in(0,0,0), (double(*)[2])&in(0,0,0));
    in = in/in.n_elem;
}


void fft::init(arma::cx_mat& in, int length) {
    
    if (length > 0 & in.n_rows != length) {
        arma::cx_mat tmp = in;
        in.zeros(length, in.n_cols);
        in(arma::span(0,tmp.n_rows-1), arma::span(0,tmp.n_cols-1)) = tmp;
    }
    
    if (fftPlan.nRows != in.n_rows || fftPlan.nCols != in.n_cols || fftPlan.nCols != 0) {
#pragma omp critical (fftPlanning)  // let's be safe
        {
            if (fftPlan.nRows != 0)
                fftw_destroy_plan(fftPlan.plan);
            fftPlan.plan = fftw_plan_dft_2d(in.n_cols, in.n_rows, (double(*)[2])&in(0,0), (double(*)[2])&in(0,0), FFTW_FORWARD, FFTW_WISDOM_ONLY);
            if (fftPlan.plan == NULL) {
                std::cout << "no wisdom, estimating" << std::endl;
                arma::cx_mat tmp = in; // fftw_plan_dft_2d may mess with in
                fftPlan.plan = fftw_plan_dft_2d(in.n_cols, in.n_rows, (double(*)[2])&in(0,0), (double(*)[2])&in(0,0), FFTW_FORWARD, fftPlanType);
                in = tmp;
            }
            fftPlan.nRows = in.n_rows;
            fftPlan.nCols = in.n_cols;
            fftPlan.nSlices = 0;
        }
    }
}

void fft::execute(arma::cx_mat& in, int length) {
    init(in, length);
    fftw_execute_dft(fftPlan.plan, (double(*)[2])&in(0,0), (double(*)[2])&in(0,0));
}

void fft::init(arma::cx_mat& in, arma::cx_mat& out, int length) {
    
    if (length > 0 & in.n_rows != length) {
        arma::cx_mat tmp = in;
        in.zeros(length, in.n_cols);
        in(arma::span(0,tmp.n_rows-1), arma::span(0,tmp.n_cols-1)) = tmp;
    }
    
    if (fftPlan.nRows != in.n_rows || fftPlan.nCols != in.n_cols || fftPlan.nCols != 0) {
#pragma omp critical (fftPlanning)  // let's be safe
        {
            if (fftPlan.nRows != 0)
                fftw_destroy_plan(fftPlan.plan);
            fftPlan.plan = fftw_plan_dft_2d(in.n_cols, in.n_rows, (double(*)[2])&in(0,0), (double(*)[2])&out(0,0), FFTW_FORWARD, FFTW_WISDOM_ONLY);
            if (fftPlan.plan == NULL) {
                std::cout << "no wisdom, estimating" << std::endl;
                arma::cx_mat tmp = in; // fftw_plan_dft_2d may mess with in
                fftPlan.plan = fftw_plan_dft_2d(in.n_cols, in.n_rows, (double(*)[2])&in(0,0), (double(*)[2])&out(0,0), FFTW_FORWARD, fftPlanType);
                in = tmp;
            }
            fftPlan.nRows = in.n_rows;
            fftPlan.nCols = in.n_cols;
            fftPlan.nSlices = 0;
        }
    }
}


void fft::execute(arma::cx_mat& in, arma::cx_mat& out, int length) {
    init(in, out, length);
    fftw_execute_dft(fftPlan.plan, (double(*)[2])&in(0,0), (double(*)[2])&out(0,0));
}

void ifft::init(arma::cx_mat& in, int length) {
    
    if (length > 0 & in.n_rows != length) {
        arma::cx_mat tmp = in;
        in.resize(length, in.n_cols);
        in(arma::span(0,tmp.n_rows-1), arma::span(0,tmp.n_cols-1)) = tmp;
    }
    
    if (fftPlan.nRows != in.n_rows || fftPlan.nCols != in.n_cols || fftPlan.nSlices != 0) {
#pragma omp critical (fftPlanning)  // let's be safe
        {
            if (fftPlan.nRows != 0)
                fftw_destroy_plan(fftPlan.plan);
            fftPlan.plan = fftw_plan_dft_2d(in.n_cols, in.n_rows, (double(*)[2])&in(0,0), (double(*)[2])&in(0,0), FFTW_BACKWARD, FFTW_WISDOM_ONLY);
            if (fftPlan.plan == NULL) {
                std::cout << "no wisdom, estimating" << std::endl;
                arma::cx_mat tmp = in; // fftw_plan_dft_2d may mess with in
                fftPlan.plan = fftw_plan_dft_2d(in.n_cols, in.n_rows, (double(*)[2])&in(0,0), (double(*)[2])&in(0,0), FFTW_BACKWARD, fftPlanType);
                in = tmp;
            }
            fftPlan.nRows = in.n_rows;
            fftPlan.nCols = in.n_cols;
            fftPlan.nSlices = 0;
        }
    }
}

void ifft::execute(arma::cx_mat& in, int length) {
    init(in, length);
    fftw_execute_dft(fftPlan.plan, (double(*)[2])&in(0,0), (double(*)[2])&in(0,0));
    in /= in.n_elem;
}

void ifft::init(arma::cx_mat& in, arma::cx_mat& out, int length) {
    
    if (length > 0 & in.n_rows != length) {
        arma::cx_mat tmp = in;
        in.resize(length, in.n_cols);
        in(arma::span(0,tmp.n_rows-1), arma::span(0,tmp.n_cols-1)) = tmp;
    }
    
    if (fftPlan.nRows != in.n_rows || fftPlan.nCols != in.n_cols || fftPlan.nSlices != 0) {
#pragma omp critical (fftPlanning)  // let's be safe
        {
            if (fftPlan.nRows != 0)
                fftw_destroy_plan(fftPlan.plan);
            fftPlan.plan = fftw_plan_dft_2d(in.n_cols, in.n_rows, (double(*)[2])&in(0,0), (double(*)[2])&out(0,0), FFTW_BACKWARD, FFTW_WISDOM_ONLY);
            if (fftPlan.plan == NULL) {
                std::cout << "no wisdom, estimating" << std::endl;
                arma::cx_mat tmp = in; // fftw_plan_dft_2d may mess with in
                fftPlan.plan = fftw_plan_dft_2d(in.n_cols, in.n_rows, (double(*)[2])&in(0,0), (double(*)[2])&out(0,0), FFTW_BACKWARD, fftPlanType);
                in = tmp;
            }
            fftPlan.nRows = in.n_rows;
            fftPlan.nCols = in.n_cols;
            fftPlan.nSlices = 0;
        }
    }
}

void ifft::execute(arma::cx_mat& in, arma::cx_mat& out, int length) {
    init(in, out, length);
    fftw_execute_dft(fftPlan.plan, (double(*)[2])&in(0,0), (double(*)[2])&out(0,0));
    out = out/out.n_elem;
}

void fft::init(arma::cx_vec& in, arma::cx_vec& out, int length) {
    
    if (length > 0 & in.n_elem != length) {
        arma::cx_vec tmp = in;
        in.zeros(length);
        in(arma::span(0,tmp.n_elem-1)) = tmp;
    }
    
    if (fftPlan.nRows != in.n_elem || fftPlan.nCols != 0 || fftPlan.nSlices != 0) {
#pragma omp critical (fftPlanning)  // let's be safe
        {
            if (fftPlan.nRows != 0)
                fftw_destroy_plan(fftPlan.plan);
            fftPlan.plan = fftw_plan_dft_1d(in.n_elem, (double(*)[2])&in(0), (double(*)[2])&out(0), FFTW_FORWARD, FFTW_WISDOM_ONLY);
            if (fftPlan.plan == NULL) {
                std::cout << "no wisdom, estimating" << std::endl;
                arma::cx_vec tmp = in; // fftw_plan_dft_1d may mess with in
                fftPlan.plan = fftw_plan_dft_1d(in.n_elem, (double(*)[2])&in(0), (double(*)[2])&out(0), FFTW_FORWARD, fftPlanType);
                in = tmp;
            }
            fftPlan.nRows = in.n_elem;
            fftPlan.nCols = 0;
            fftPlan.nSlices = 0;
        }
    }
}

void fft::execute(arma::cx_vec& in, arma::cx_vec& out, int length) {
    init(in, out, length);
    fftw_execute_dft(fftPlan.plan, (double(*)[2])&in(0), (double(*)[2])&out(0));
}

void ifft::init(arma::cx_vec& in, arma::cx_vec& out, int length) {
    
    if (length > 0 & in.n_elem != length) {
        arma::cx_vec tmp = in;
        in.zeros(length);
        in(arma::span(0,tmp.n_elem-1)) = tmp;
    }
    
    if (fftPlan.nRows != in.n_elem || fftPlan.nCols != 0 || fftPlan.nSlices != 0) {
#pragma omp critical (fftPlanning)  // let's be safe
        {
            if (fftPlan.nRows != 0)
                fftw_destroy_plan(fftPlan.plan);
            fftPlan.plan = fftw_plan_dft_1d(in.n_elem, (double(*)[2])&in(0), (double(*)[2])&out(0), FFTW_BACKWARD, FFTW_WISDOM_ONLY);
            if (fftPlan.plan == NULL) {
                std::cout << "no wisdom, estimating" << std::endl;
                arma::cx_vec tmp = in; // fftw_plan_dft_1d may mess with in
                fftPlan.plan = fftw_plan_dft_1d(in.n_elem, (double(*)[2])&in(0), (double(*)[2])&out(0), FFTW_BACKWARD, fftPlanType);
                in = tmp;
            }
            fftPlan.nRows = in.n_elem;
            fftPlan.nCols = 0;
            fftPlan.nSlices = 0;
        }
    }
}

void ifft::execute(arma::cx_vec& in, arma::cx_vec& out, int length) {
    init(in, out, length);
    fftw_execute_dft(fftPlan.plan, (double(*)[2])&in(0), (double(*)[2])&out(0));
    out /= out.n_elem;
}

//////////////////////////////////////

arma::cx_cube fft_shift(arma::cx_cube& in) {
    arma::cx_cube shiftE = in;
    
    for (int k=0; k<shiftE.n_slices; k++) {
        shiftE.slice(k) = shift(shiftE.slice(k), floor(shiftE.slice(k).n_rows/2), 0);
        shiftE.slice(k) = shift(shiftE.slice(k), floor(shiftE.slice(k).n_cols/2), 1);
    }
    return shiftE;
}

arma::mat fft_shift(arma::mat& in) {
    arma::mat shiftE = in;
    
    shiftE = shift(shiftE, floor(shiftE.n_rows/2), 0);
    shiftE = shift(shiftE, floor(shiftE.n_cols/2), 1);
    return shiftE;
}

arma::cx_mat fft_shift(arma::cx_mat& in) {
    arma::cx_mat shiftE = in;
    
    shiftE = shift(shiftE, floor(shiftE.n_rows/2), 0);
    shiftE = shift(shiftE, floor(shiftE.n_cols/2), 1);
    return shiftE;
}

//////////////////////////////////////

void fresnelPropagateAS::execute(arma::cx_cube& in, double *lambda, double apRad, double z) {
    std::complex<double> r1(1, 0); // for conversion from real to complex
    std::complex<double> i1(0, 1); // for conversion from real to complex
    int N = in.n_rows>in.n_cols?in.n_rows:in.n_cols;
    
    // Computation of the transfer function spectrum frequency squared f2 = fx^2 + fy^2
    H = arma::zeros<arma::mat>(in.n_rows, in.n_cols);
    vH = arma::regspace<arma::vec>(-N/2,N/2-1);
    H.each_row() = vH.t();
    H = square(H) + square(H.t());
    H = fft_shift(H);
    
//    myFft.execute(in);
    for (int sl=0; sl<in.n_slices; sl++) {
        double k = 2.*M_PI/lambda[sl];
        arma::cx_mat cH = exp(i1*k*z * (sqrt( 1. - 1./4. * pow(lambda[sl]/apRad,2.) * H ))); // Goodman 4-20
        myFft.execute(in.slice(sl));
        in.slice(sl) = in.slice(sl)%cH;
        myIfft.execute(in.slice(sl));
    }
//    myIfft.execute(in);
}

//////////////////////////////////////

arma::cx_cube czt::execute(arma::cx_cube& in, int outN, arma::vec pRatio, arma::vec start) {
    return execute(in, outN, pRatio, start, pRatio, start);
}

arma::cx_cube czt::execute(arma::cx_cube& in, int outN, arma::vec pRatioX, arma::vec startX, arma::vec pRatioY, arma::vec startY) {
    arma::cx_cube out(outN, outN, in.n_slices);
    
    for (int s=0; s<in.n_slices; s++) {
//        std::cout << "czt::execute cube: slice " << s << std::endl;
        out.slice(s) = execute(in.slice(s), outN, pRatioX[s], startX[s], pRatioY[s], startY[s]);
    }
    return out;
}

arma::cx_mat czt::execute(arma::cx_mat& in, int outN, double pRatio, double start) {
    return execute(in, outN, pRatio, start, pRatio, start);
}

arma::cx_mat czt::execute(arma::cx_mat& in, int outN, double pRatioX, double startX, double pRatioY, double startY) {
    arma::cx_mat out(outN, outN);
//    std::cout << "czt2d: in: " << size(in) << " outN: " << outN << " pRatioX: " << pRatioX << " startX: " << startX << " pRatioY: " << pRatioY << " startY: " << startY << std::endl;
    
    czt2d(in, in.n_rows, outN, pRatioX, startX, pRatioY, startY, out);
    return out;
}

int czt::czt2d(arma::cx_mat& in, int m, int k, double wx, double ax, double wy, double ay, arma::cx_mat& out) {
    // symmetric case: square input, same {k,w,a} in both dims
    // trivial mod for different {m,k,w,a} in each dim
//    std::cout << "czt2d: in: " << size(in) << " m: " << m << " k: " << k << std::endl;
//    printf("wx = %0.15f\n", wx);
//    printf("ax = %0.15f\n", ax);
//    printf("wy = %0.15f\n", wy);
//    printf("ay = %0.15f\n", ay);
    
    // matlab: out = czt_ph(czt_ph(Ein, length(eta), Wy, Ay).',length(xi), Wx, Ax).';
    arma::cx_mat out1;
    czt2d(in,k,wy,ay,out1, false);
    out1 = out1.st(); // non-conjugate transpose
    
    czt2d(out1,k,wx,ax,out, false);
    out = out.st(); // non-conjugate transpose

    return 0;
}

int czt::czt2d(arma::cx_mat& in, int k, double w, double a, arma::cx_mat& out, bool diagnostics) {
    std::complex<double> i1(0, 1);
    int nfft;
    arma::vec ww;
    arma::vec aa;
    arma::vec kk;
    arma::cx_vec caa;
    arma::cx_vec iww;
    arma::cx_mat y;
    arma::cx_mat fy;
    arma::cx_vec fv;
    arma::cx_mat g;
    fftw_plan py2fy = NULL;
    fftw_plan piww2fv = NULL;
    fftw_plan ipfy2g = NULL;
    arma::wall_clock timer;

    
    int m = in.n_rows;
    int n = in.n_cols;
    
    // matlab: nfft = 2^nextpow2(m+k-1);
    nfft = pow(2,ceil(log2(m+k-1)));  // next largest power of 2, for ffts
//    std::cout << "m = " << m << ", n = " << n << ", k = " << k << ", nfft = " << nfft << std::endl;
//    printf("w = %0.30f\n", w);
//    printf("a = %0.30f\n", a);
    
    // set mm=m = max(m,k)
    int mm = 0;
    if (m>=k)
        mm=m;
    else
        mm=k;

//    timer.tic();
#pragma omp critical (fftPlanning)  // let's be safe
    {
        y = arma::zeros<arma::cx_mat>(nfft, n);
        fy = arma::cx_mat(nfft, n);
        fv = arma::cx_vec(nfft);
        g = arma::cx_mat(nfft, n);
        ww.set_size(m+mm-1);
        iww = arma::zeros<arma::cx_vec>(nfft);
        kk = arma::vec(ww.n_elem);
        if (diagnostics)
            std::cout << "finished allocation" << std::endl;
        
        py2fy = fftw_plan_dft_1d(nfft,(double(*)[2])&y(0,0),(double(*)[2])&fy(0,0),FFTW_FORWARD,FFTW_WISDOM_ONLY);
        if (py2fy == NULL) {
            std::cout << "czt py2fy: no wisdom, estimating" << std::endl;
            arma::cx_mat tmp = y; // fftw_plan_dft_1d may mess with in
            py2fy = fftw_plan_dft_1d(nfft,(double(*)[2])&y(0,0),(double(*)[2])&fy(0,0),FFTW_FORWARD,fftPlanType);
            y = tmp;
        }
        piww2fv = fftw_plan_dft_1d(nfft,(double(*)[2])&iww(0),(double(*)[2])&fv(0),FFTW_FORWARD,FFTW_WISDOM_ONLY);
        if (piww2fv == NULL) {
            std::cout << "czt piww2fv: no wisdom, estimating" << std::endl;
            arma::cx_vec tmp = iww; // fftw_plan_dft_3d may mess with in
            piww2fv = fftw_plan_dft_1d(nfft,(double(*)[2])&iww(0),(double(*)[2])&fv(0),FFTW_FORWARD,fftPlanType);
            iww = tmp;
        }
        ipfy2g = fftw_plan_dft_1d(nfft,(double(*)[2])&fy(0,0),(double(*)[2])&g(0,0),FFTW_BACKWARD,FFTW_WISDOM_ONLY);
        if (ipfy2g == NULL) {
            std::cout << "czt ipfy2g: no wisdom, estimating" << std::endl;
            arma::cx_mat tmp = fy; // fftw_plan_dft_3d may mess with in
            ipfy2g = fftw_plan_dft_1d(nfft,(double(*)[2])&fy(0,0),(double(*)[2])&g(0,0),FFTW_BACKWARD,fftPlanType);
            fy = tmp;
        }
    }
//    std::cout << "=========== setting up mem and fft time: " << timer.toc() << " seconds" << std::endl;
    
//    timer.tic();
    // matlab: kk = ( (-m+1):max(k-1,m-1) ).';
    // matlab: kk2 = (kk .^ 2) ./ 2;
    // matlab: ww = w .* (kk2);
    for (double i=0.0;i<m+mm-1;i+=1.0){
        // argument will go from -m+1 to m+mm-2-m+1=mm-1=max(m-1,k-1)
        kk(i) = i-m+1.0;
        ww(i) = w*(kk(i)*kk(i)/2.0);
    }
//    std::cout << "length of ww: " << ww.n_elem << std::endl;

    aa.set_size(m);
    // matlab: nn = (0:(m-1))';
    // matlab: aa = a .* ( -nn );
    for (int i=0;i<m;++i)
        aa(i) = a*((double) -i);
    // matlab: aa = exp(i*(aa + ww(m+nn)));
    caa.set_size(m);
    arma::vec dww(m);
    arma::vec aww(m);
    for (int i=0;i<m;++i) {
        dww(i) = ww(m-1+i);
        aww(i) = aa(i) + ww(m-1+i);
        caa(i) = exp(i1*(aa(i) + ww(m-1+i)));
    }
//    std::cout << "=========== setting up w and a time: " << timer.toc() << " seconds" << std::endl;
    
    if (diagnostics) {
        for (int i=0; i<5; i++)
            printf("ww(%d) = %0.15f\n", i, ww(i));
        for (int i=ww.n_elem - 5; i<ww.n_elem; i++)
            printf("ww(%d) = %0.15f\n", i, ww(i));
        for (int i=0; i<5; i++)
            printf("aa(%d) = %0.15f\n", i, aa(i));
        for (int i=aa.n_elem - 5; i<aa.n_elem; i++)
            printf("aa(%d) = %0.15f\n", i, aa(i));
        for (int i=0; i<5; i++)
            printf("caa(%d) = %0.15f + %0.15fi\n", i, real(caa(i)), imag(caa(i)));
        for (int i=aa.n_elem - 5; i<aa.n_elem; i++)
            printf("caa(%d) = %0.15f + %0.15fi\n", i, real(caa(i)), imag(caa(i)));
        
        save_mat("ww.fits", ww);
        save_mat("dww.fits", dww);
        save_mat("aww.fits", aww);
        save_mat("aa.fits", aa);
        save_mat("caa.fits", caa);
        save_mat("kk.fits", kk);
    }
//    y.zeros(nfft, n); // zero pad
    
    // matlab: fv = fft( exp(-i * ww(1:(k-1+m))), nfft );
    // fft2.execute(iww, fv);
//    timer.tic();
    iww(arma::span(0, (double)(m-1+k-1))) = exp(-i1*ww(arma::span(0,(double)(m-1+k-1)))); // chirp filter
//    std::cout << "length of iww: " << iww.n_elem << std::endl;

    fftw_execute_dft(piww2fv, (double(*)[2])&iww(0), (double(*)[2])&fv(0));
    if (diagnostics) {
        std::cout << "fv:" << std::endl;
        fv(arma::span(0, 4)).print();
        std::cout << "finished fft2" << std::endl;
        save_mat("testfv.fits", fv);
        save_mat("testiww.fits", iww);
    }
//    std::cout << "=========== computing fv time: " << timer.toc() << " seconds" << std::endl;
    
    
    // matlab: y = x .* aa(:,ones(1,n));
    // matlab: fy = fft(  y, nfft );
    // matlab: fy = fy .* fv(:,ones(1, n));
    // matlab: g  = ifft( fy );
//    timer.tic();
#pragma omp parallel
    {
#pragma omp for
        for (int c=0; c<y.n_cols; ++c) {
            y(arma::span(0,m-1),c) = in(arma::span::all,c) % caa;
            fftw_execute_dft(py2fy, (double(*)[2])&y(0,c), (double(*)[2])&fy(0,c));
        }
    }
    if (diagnostics) {
        save_mat("testfy1.fits", fy);
    }
#pragma omp parallel
    {
#pragma omp for
        for (int c=0; c<y.n_cols; ++c) {
            fy(arma::span(0,fy.n_rows-1),c) %= fv;
            fftw_execute_dft(ipfy2g, (double(*)[2])&fy(0,c), (double(*)[2])&g(0,c));
        }
    }
    g /= g.n_rows;
//    std::cout << "=========== computing g time: " << timer.toc() << " seconds" << std::endl;
    
    if (diagnostics) {
        save_mat("testin.fits", in);
        save_mat("testy.fits", y);
    }
    if (diagnostics) {
        std::cout << "fy after fft2 mult:" << std::endl;
        fy(arma::span(0, 4), arma::span(0,4)).print();
        std::cout << "finished post fft2 mult" << std::endl;
        save_mat("testfy2.fits", fy);
    }
//    std::cout << "shape of g before trim: " << size(g) << ", n = " << n << std::endl;
//    timer.tic();
//    g = g(arma::span::all, arma::span(0,n-1));
//    std::cout << "shape of g after trim: " << size(g) << std::endl;
//    std::cout << "=========== trimming g time: " << timer.toc() << " seconds" << std::endl;
    if (diagnostics) {
        std::cout << "shape of g: " << size(g) << std::endl;
        std::cout << "g:" << std::endl;
        g(arma::span(0, 4), arma::span(0,4)).print();
        save_mat("testg.fits", g);
    }
    
    // matlab: g = g( m:(m+k-1), : ) .* exp(i*ww( m:(m+k-1),ones(1, n)) );
//    arma::cx_vec expww = exp(i1*ww(arma::span(m-1,m+k-2)));
//    timer.tic();
    arma::cx_vec expww(k);
    #pragma omp parallel
    {
        #pragma omp for
        for (int i=0; i<k; ++i) {
            expww(i) = exp(i1*ww(m-1+i));
        }
    }
//    std::cout << "=========== computing expww time: " << timer.toc() << " seconds" << std::endl;
//    std::cout << "created expww, expww = " << size(expww) << std::endl;
    out.set_size(size(g(arma::span(m-1,m+k-2), arma::span::all)));
//    std::cout << "created out, out = " << size(out) << std::endl;
//    timer.tic();
    #pragma omp parallel
    {
        #pragma omp for
        for (int c=0; c<out.n_cols; ++c) {
            out(arma::span::all,c) = g(arma::span(m-1,m+k-2), c) % expww;
        }
    }
//    std::cout << "=========== computing out time: " << timer.toc() << " seconds" << std::endl;

    if (diagnostics) {
        std::cout << "shape of out: " << size(out) << std::endl;
        std::cout << "out:" << std::endl;
        out(arma::span(0, 4), arma::span(0,4)).print();
        save_mat("testout.fits", out);
    }
    return 0;
}

//////////////////////////////////////

// init geometry with physical units
void zoomFft::init(arrayGeom agIn, arrayGeom agOut, double *lambdaFocalLength, int nLambdas) {
    std::complex<double> i1(0, 1);
    Nx = agOut.pixelX.n_elem;
    Ny = agOut.pixelY.n_elem;
    inPixSizeX = agIn.pixelSizeX;
    inPixSizeY = agIn.pixelSizeY;
    
    lFl.set_size(nLambdas);
    Ax.set_size(nLambdas);
    Wx.set_size(nLambdas);
    Ay.set_size(nLambdas);
    Wy.set_size(nLambdas);
    
    for (int i=0; i<nLambdas; i++) {
        lFl[i] = lambdaFocalLength[i];
        Ax[i] = 2*M_PI*agIn.pixelSizeX*agOut.pixelX[0]/lFl[i];
        Wx[i] = -2*M_PI*agIn.pixelSizeX*agOut.pixelSizeX/lFl[i];
        if (agIn.pixelSizeX != agIn.pixelSizeY | agOut.pixelY[0] != agOut.pixelSizeY) {
            Ay[i] = 2*M_PI*agIn.pixelSizeY*agOut.pixelY[0]/lFl[i];
            Wy[i] = -2*M_PI*agIn.pixelSizeY*agOut.pixelSizeY/lFl[i];
        } else {
            Ay[i] = Ax[i];
            Wy[i] = Wx[i];
        }
        //        std::cout << "Ax Ay Wx Wy" << Ax << " " << Ay << " " << Wx << " " << Wy << std::endl;
    }
    //    agIn.print("agIn");
    //    agOut.print("agOut");
    mxPmy = agIn.pixelX[0]*arma::repmat(agOut.pixelX.t(), agOut.pixelY.n_elem, 1) + agIn.pixelY[0]*arma::repmat(agOut.pixelY, 1, agOut.pixelX.n_elem);
//    std::cout << "size of mxPmy: " << arma::size(mxPmy) << std::endl;
    /*
     for (int i=0; i<10; i++) {
     for (int j=0; j<10; j++) {
     printf("%0.20f ", mxPmy(i,j));
     }
     printf("\n");
     }
     */
}

// init geometry with normalized units
void zoomFft::init(int nx, int ny, double zoomFactorIn, int dirIn, int nLambdas) {
    std::complex<double> i1(0, 1);
    
    std::cout << "in zoomFft init" << std::endl;
    std::cout << "nx = " << nx << ", ny = " << ny << std::endl;
    
    dir = dirIn;
    zoomFactor = zoomFactorIn;
    Nx = nx;
    Ny = ny;
    
    lFl.set_size(nLambdas);
    Ax.set_size(nLambdas);
    Wx.set_size(nLambdas);
    Ay.set_size(nLambdas);
    Wy.set_size(nLambdas);
    
    double dx = 1.0/nx;
    double dy = 1.0/ny;
    
    // xIn = i/nx - 0.5
    arma::vec xIn;
    arma::vec yIn;
    arma::vec xOut;
    arma::vec yOut;
    if (dir == 1) {
        xIn = arma::linspace<arma::vec>(-0.5, 0.5-dx, nx);
        yIn = arma::linspace<arma::vec>(-0.5, 0.5-dy, ny);
        xOut = nx*xIn/zoomFactor;
        yOut = ny*yIn/zoomFactor;
    } else {
        xOut = arma::linspace<arma::vec>(-0.5, 0.5-dx, nx);
        yOut = arma::linspace<arma::vec>(-0.5, 0.5-dy, ny);
        xIn = nx*xOut/zoomFactor;
        yIn = ny*yOut/zoomFactor;
    }
    
    inPixSizeX = xIn[1] - xIn[0];
    inPixSizeY = yIn[1] - yIn[0];
    double outPixSizeX = xOut[1] - xOut[0];
    double outPixSizeY = yOut[1] - yOut[0];
    
    Ax[0] = dir*2*M_PI*inPixSizeX*xOut[0];
    Wx[0] = -dir*2*M_PI*inPixSizeX*outPixSizeX;
    if (nx != ny) {
        Ay[0] = dir*2*M_PI*inPixSizeY*yOut[0];
        Wy[0] = -dir*2*M_PI*inPixSizeY*outPixSizeY;
    } else {
        Ay[0] = Ax[0];
        Wy[0] = Wx[0];
    }
//    std::cout << "Ax Ay Wx Wy " << Ax[0] << " " << Ay[0] << " " << Wx[0] << " " << Wy[0] << std::endl;
    
    for (int i=0; i<Ax.n_elem; i++) {
        lFl[i] = 1.0*dir;
        Ax[i] = Ax[0];
        Ay[i] = Ay[0];
        Wx[i] = Wx[0];
        Wy[i] = Wy[0];
    }
    
//    std::cout << "Ax Ay Wx Wy " << Ax << " " << Ay << " " << Wx << " " << Wy << std::endl;
    
    mxPmy = xIn[0]*arma::repmat(xOut.t(), yOut.n_elem, 1) + yIn[0]*arma::repmat(yOut, 1, xOut.n_elem);
    /*
     for (int i=0; i<10; i++) {
     for (int j=0; j<10; j++) {
     printf("%0.20f ", mxPmy(i,j));
     }
     printf("\n");
     }
     */
}


arma::cx_cube zoomFft::execute(arma::cx_cube& in) {
    std::complex<double> i1(0, 1);
    
    assert(Nx > 0); // make sure we're initialized
    
    arma::cx_cube cztOut = propCzt.execute(in, Nx, Wx, Ax, Wy, Ay);
    for (int sl=0; sl<cztOut.n_slices; sl++) {
//        cztOut.slice(sl) %= exp(-2*M_PI*i1*mxPmy/lFl[sl])*inPixSizeX*inPixSizeY/(i1*lFl[sl]);
        std::complex<double> f1 = -2*M_PI*i1/lFl[sl];
        std::complex<double> f2 = inPixSizeX*inPixSizeY/(i1*lFl[sl]);
#pragma omp parallel
        {
#pragma omp for
            for (int c=0; c<cztOut.n_cols; ++c) {
                for (int r=0; r<cztOut.n_rows; ++r)
                     cztOut(r, c, sl) *= exp(f1*mxPmy(r, c))*f2;
            }
        }
    }
    return cztOut;
}

arma::cx_mat zoomFft::execute(arma::cx_mat& in, int waveIndex) {
    std::complex<double> i1(0, 1);
    
    assert(Nx > 0); // make sure we're initialized

    arma::cx_mat cztOut = propCzt.execute(in, Nx, Wx[waveIndex], Ax[waveIndex], Wy[waveIndex], Ay[waveIndex]);
    std::complex<double> f1 = -2*M_PI*i1/lFl[waveIndex];
    std::complex<double> f2 = inPixSizeX*inPixSizeY/(i1*lFl[waveIndex]);
#pragma omp parallel
    {
#pragma omp for
        for (int c=0; c<cztOut.n_cols; ++c) {
            for (int r=0; r<cztOut.n_rows; ++r)
                cztOut(r, c) *= exp(f1*mxPmy(r, c))*f2;
        }
    }

    return cztOut;
}

//////////////////////////////////////
//////////////////////////////////////

// zoomDft impements Olivier's fft_Dft
//
// For the complex input values z, computes pointwise
// the complex output value
// val_j = (sum_k{z_k*exp(2*dir*M_PI*i*(xout_j*xin_k + yout_j*yin_k))})/Zfactor
// where xin_j = j/xsize - 0.5 and xout = xsize*xin/Zfactor
// here j, k are linear indices into the 2D arrays xin etc.

// zoomFactor is zoom factor
// dir = -1 for FT, 1 for inverse FT
void zoomDft::init(int nx, int ny, double zoomFactorIn, int dirIn) {

    std::cout << "in zoomDft init" << std::endl;
    std::cout << "nx = " << nx << ", ny = " << ny << std::endl;
    
    dir = dirIn;
    zoomFactor = zoomFactorIn;
    
    xIn.set_size(nx, nx);
    xOut.set_size(nx, nx);
    yIn.set_size(ny, ny);
    yOut.set_size(ny, ny);
    
    double dx = 1.0/nx;
    double dy = 1.0/ny;
    
    // xIn = i/nx - 0.5
    arma::vec x = arma::linspace<arma::vec>(-0.5, 0.5-dx, nx);
    arma::vec y = arma::linspace<arma::vec>(-0.5, 0.5-dy, nx);
    
    xIn = arma::repmat(x.t(), ny, 1);
    yIn = arma::repmat(y, 1, nx);
    std::cout << "xIn = " << xIn.size() << ", yIn = " << yIn.size() << std::endl;
    std::cout << "xIn = " << xIn.n_rows << "x" << yIn.n_cols << std::endl;
    
    // xOut = nx*(i/nx - 0.5)/Zfactor = nx*xIn/Zfactor
    xOut = nx*xIn/zoomFactor;
    yOut = ny*yIn/zoomFactor;
    
}

arma::cx_cube  zoomDft::execute(arma::cx_cube& in) {
    arma::cx_cube out(in.n_rows, in.n_cols, in.n_slices);
    
    for (int s=0; s<in.n_slices; s++) {
        std::cout << "zoomDft::execute cube: slice " << s << std::endl;
        out.slice(s) = execute(in.slice(s));
    }
    
    return out;
}

arma::cx_mat zoomDft::execute(arma::cx_mat& in) {
    std::complex<double> i1(0, 1);
    arma::cx_mat out(in.n_rows, in.n_cols);
    
#pragma omp parallel for
    for (int i=0; i<in.n_rows; i++) {
        if (i%20 == 0) {
            std::cout << i << " ";
            fflush(stdout);
        }
        for (int j=0; j<in.n_cols; j++) {
            out(i,j) = arma::sum(arma::sum(in%exp(2*dir*M_PI*i1*(xOut(i,j)*xIn + yOut(i,j)*yIn))))/zoomFactor;
        }
    }
    if (dir == 1)
        out = out/(out.n_rows*out.n_cols);
    return out;
}


