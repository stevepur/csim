//
//  celem.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef fpmPupToLyot_hpp
#define fpmPupToLyot_hpp

#include "celem.hpp"
#include "complexHexMaskFPM.hpp"
#include "../lib/csim_lib.hpp"
#include "armadillo"

class fpmPupToLyot;

// helper class for handling different FPM mask types
// sets data in fpmPupToLyot
class fpmForPupilToLyot {
protected:
    
public:
    fpmForPupilToLyot() {}
    fpmForPupilToLyot(fpmPupToLyot *p2l, initCommandSet*& cmdBlock) {}
    
    virtual void set_geometry(fpmPupToLyot *prop, efield* E, double *lambda, double *lambdaFocalLength) {}
    virtual void set_fpmMatAmp(fpmPupToLyot *p2l, double lambda, int sl) {}
    virtual void apply_babinet(fpmPupToLyot *p2l) {}
    virtual void draw(const char *title = ""){}
    virtual void get_optimization_data(const char *dataName, void *data) {}
    virtual void set_optimization_data(const char *dataName, void *data) {}
    
};

////////

class fpmCMCForPupToLyot : public fpmForPupilToLyot {
    arma::cube complexMaskMatAmp;
    arma::cube complexMaskMatPh;
    arma::cx_cube complexMaskCube;
    int maskIndex = -1;
    double fpmMatScale = 0;
    
public:
    fpmCMCForPupToLyot(fpmPupToLyot *p2l, initCommandSet*& cmdBlock);
    
    void set(fpmPupToLyot *p2l, std::string fieldName, const char *arg);
    void initMask(const char *filenameAmp, const char *filenamePh);
    void set_geometry(fpmPupToLyot *prop, efield* E, double *lambda, double *lambdaFocalLength);
    void set_fpmMatAmp(fpmPupToLyot *p2l, double lambda, int sl);
    void apply_babinet(fpmPupToLyot *p2l);
    void draw(const char *title = "");
};

////////

class fpmIntHexCMCForPupToLyot : public fpmForPupilToLyot {
    complexHexMaskFPM *hexFPM;
    double zoomFactor = -1;
    
public:
    fpmIntHexCMCForPupToLyot(fpmPupToLyot *p2l, initCommandSet*& cmdBlock);
    
    void set(fpmPupToLyot *p2l, std::string fieldName, const char *arg);
    void initMask(const char *filenameAmp, const char *filenamePh) {}
    void set_geometry(fpmPupToLyot *prop, efield* E, double *lambda, double *lambdaFocalLength);
    void set_fpmMatAmp(fpmPupToLyot *p2l, double lambda, int sl);
    
    void get_optimization_data(const char *dataName, void *data);
    void set_optimization_data(const char *dataName, void *data);
    
    void apply_babinet(fpmPupToLyot *p2l);
    void draw(const char *title = "");
};

////////

class fpmBinaryForPupToLyot : public fpmForPupilToLyot {
    arma::mat fpmMat;
    arma::mat calibMat;
    double referenceLambda = 0;
    double innerRadiusFld = 0;
    double outerRadiusFld = 0;
    
public:
    fpmBinaryForPupToLyot(fpmPupToLyot *p2l, initCommandSet*& cmdBlock);
    
    void set(fpmPupToLyot *p2l, std::string fieldName, const char *arg);
    void set_geometry(fpmPupToLyot *prop, efield* E, double *lambda, double *lambdaFocalLength);
    void set_fpmMatAmp(fpmPupToLyot *p2l, double lambda, int sl);
    void apply_babinet(fpmPupToLyot *p2l) {}
    void draw(const char *title = "");
};


//////////////////////////////////////////

class fpmPupToLyot : public celem {
    friend class fpmForPupilToLyot;
    friend class fpmCMCForPupToLyot;
    friend class fpmIntHexCMCForPupToLyot;
    friend class fpmBinaryForPupToLyot;

    arma::cx_mat fpmMatAmp;
    arma::cx_mat maskHat;
    arma::cx_mat fftMaskHat;
    
    arma::cx_mat paddedE;
    arma::cx_mat fftPaddedE;
    arma::cx_mat fftMHatTimesfftPaddedE;
    
    zoomFft propZoomFft;

    arrayGeom maskGeom;
    arrayGeom paddedGeom;
    
    fpmForPupilToLyot *mask;
    
    fft maskFft;
    fft paddedEFft;
    ifft myIfft;
    
    double focalRatio = -1; // focal length at this point in the coronagraph
    double fRatioSign = 1; // sign convention for fratio in zoomFFT
    double zoomFftSign = 1; // sign convention for zoomFFT

    
public:
    fpmPupToLyot();
    fpmPupToLyot(const char *inName);
    fpmPupToLyot(initCommandSet*& cmdBlock);
    
    void init(initCommandSet*& cmdBlock);
    void init_block(initCommandSet*& cmdBlock);
    
    void set(std::string fieldName, const char *arg);
    
    efield* execute(efield* E, celem* prev, celem* next, double time);
    
    void get_optimization_data(const char *dataName, void *data);
    void set_optimization_data(const char *dataName, void *data);

    void draw(const char *title = "");
};

#endif /* fpmPupToLyot_hpp */
