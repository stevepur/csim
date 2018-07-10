//
//  csim_fft.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef csim_fft_hpp
#define csim_fft_hpp

#include <fftw3.h>
#include "armadillo"
#include "csim_geom.hpp"

void init_fft_lib(bool computeFftWisdom = false);
void save_fft_wisdom(void);

class fftPlanData {
public:
    fftw_plan plan = NULL;
    int nRows = 0;
    int nCols = 0;
    int nSlices = 0;
    
    ~fftPlanData() {if (plan != NULL) fftw_destroy_plan(plan); };
};

class fft {
    fftPlanData fftPlan;
public:
    fft() {};
    void init(arma::cx_cube& in, int length = 0);
    void init(arma::cx_mat& in, int length = 0);
    void init(arma::cx_mat& in, arma::cx_mat& out, int length = 0);
    void init(arma::cx_vec& in, arma::cx_vec& out, int length = 0);
    void execute(arma::cx_cube& in, int length = 0);
    void execute(arma::cx_mat& in, int length = 0);
    void execute(arma::cx_mat& in, arma::cx_mat& out, int length = 0);
    void execute(arma::cx_vec& in, arma::cx_vec& out, int length = 0);
};

class ifft {
    fftPlanData fftPlan;
public:
    ifft() {};
    void init(arma::cx_cube& in, int length = 0);
    void init(arma::cx_mat& in, int length = 0);
    void init(arma::cx_mat& in, arma::cx_mat& out, int length = 0);
    void init(arma::cx_vec& in, arma::cx_vec& out, int length = 0);
    void execute(arma::cx_cube& in, int length = 0);
    void execute(arma::cx_mat& in, int length = 0);
    void execute(arma::cx_mat& in, arma::cx_mat& out, int length = 0);
    void execute(arma::cx_vec& in, arma::cx_vec& out, int length = 0);
};

arma::cx_cube fft_shift(arma::cx_cube& in);
arma::mat fft_shift(arma::mat& in);
arma::cx_mat fft_shift(arma::cx_mat& in);

class fresnelPropagateAS {
    fft myFft;
    ifft myIfft;
    
    arma::mat H;
    arma::mat vH;
    
public:
    fresnelPropagateAS() {};
    void execute(arma::cx_cube& in, double *lambda, double apRad, double z);
};

class czt {
//    fft fft1;
//    fft fft2;
//    ifft ifft1;
    
    
    int n;
    int nfft;

public:
    czt(){};
    arma::cx_cube execute(arma::cx_cube& in, int outN, arma::vec pRatio, arma::vec start);
    arma::cx_cube execute(arma::cx_cube& in, int outN, arma::vec pRatioX, arma::vec startX, arma::vec pRatioY, arma::vec startY);
    arma::cx_mat execute(arma::cx_mat& in, int outN, double pRatio, double start);
    arma::cx_mat execute(arma::cx_mat& in, int outN, double pRatioX, double startX, double pRatioY, double startY);
    int czt2d(arma::cx_mat& in, int m, int k, double wx, double ax, double wy, double ay, arma::cx_mat& out);
    int czt2d(arma::cx_mat& in, int k, double w, double a, arma::cx_mat& out, bool diagnostics = false);
};

class zoomFft {
    czt propCzt;
    
    int Nx = 0;
    int Ny = 0;
    double inPixSizeX = 0;
    double inPixSizeY = 0;
    
    arma::vec Ax;
    arma::vec Wx;
    arma::vec Ay;
    arma::vec Wy;
    
    arma::mat mxPmy;
    arma::mat my;
    
    arma::vec lFl;
    
    int dir;
    double zoomFactor;
    
public:
    zoomFft() {};
    void init(arrayGeom agIn, arrayGeom agOut, double *lambdaFocalLength, int nLambdas);
    void init(int nx, int ny, double zoomFactorIn, int dirIn, int nLambdas);
    arma::cx_cube execute(arma::cx_cube& in);
    arma::cx_mat execute(arma::cx_mat& in, int waveIndex = 0);
};

/*
class czt {
    //    fft fft1;
    //    fft fft2;
    //    ifft ifft1;
    
    
    int n;
    int nfft;
    
public:
    czt(){};
    arma::cx_cube execute(arma::cx_cube& in, int outN, arma::cx_vec pRatio, arma::cx_vec start);
    arma::cx_cube execute(arma::cx_cube& in, int outN, arma::cx_vec pRatioX, arma::cx_vec startX, arma::cx_vec pRatioY, arma::cx_vec startY);
    arma::cx_mat execute(arma::cx_mat& in, int outN, std::complex<double> pRatio, std::complex<double> start);
    arma::cx_mat execute(arma::cx_mat& in, int outN, std::complex<double> pRatioX, std::complex<double> startX, std::complex<double> pRatioY, std::complex<double> startY);
    int czt2d(arma::cx_mat& in, int m, int k, std::complex<double> wx, std::complex<double> ax, std::complex<double> wy, std::complex<double> ay, arma::cx_mat& out);
    int czt1d(arma::cx_vec& x, int m, int k, std::complex<double> w, std::complex<double> a, arma::cx_vec& out, bool doPrint = false);
};

 class zoomFft {
    czt propCzt;
    
    int Nx = 0;
    int Ny = 0;
    double inPixSizeX = 0;
    double inPixSizeY = 0;
    
    arma::cx_vec Ax;
    arma::cx_vec Wx;
    arma::cx_vec Ay;
    arma::cx_vec Wy;
    
    arma::mat mxPmy;
    arma::mat my;
    
    arma::vec lFl;

    int dir;
    double zoomFactor;

public:
    zoomFft() {};
    void init(arrayGeom agIn, arrayGeom agOut, double *lambdaFocalLength, int nLambdas);
    void init(int nx, int ny, double zoomFactorIn, int dirIn, int nLambdas);
    arma::cx_cube execute(arma::cx_cube& in);
    arma::cx_mat execute(arma::cx_mat& in, int waveIndex = 0);
};
 */

class zoomDft {
    arma::mat xIn;
    arma::mat xOut;
    arma::mat yIn;
    arma::mat yOut;
    
    int dir;
    double zoomFactor;
    
public:
    zoomDft() {};
    void init(int nx, int ny, double zoomFactor, int dir);
    arma::cx_cube execute(arma::cx_cube& in);
    arma::cx_mat execute(arma::cx_mat& in);
};

#endif /* csim_fft_hpp */
