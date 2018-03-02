//
//  complexHexMaskFPM.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef complexHexMaskFPM_hpp
#define complexHexMaskFPM_hpp

#include "celem.hpp"
#include "armadillo"

class complexHexMaskFPM : public celem {
    friend class fpmIntHexCMCForPupToLyot;
    
    double hexStep = -1;
    int maskNRings = -1;
    int fpmArraySize = -1;
    double fpmScaleFactor = -1;
    double hexGap = -1;
    double lambdaRef = -1;
    double fpmFRatio = -1;
    int nSubPix = -1;
    double fpmRadius = -1;
    
    double fpmScale = -1;
    double zoomFactor = 1.0;
    
    std::string maskFilenameRoot;
    
    arma::vec hexX;
    arma::vec hexY;
    arma::ivec hexNum;
    arma::ivec hexRing;
    
    arma::vec fpmSags;
    arma::mat fpmDesignArray;
    arma::mat fpmDesignSagArray;
    
    arma::cube interpNSubPix;
    arma::cube interpNSubPixToUse;
    arma::cube interpSagVals;
    arma::cube interpHexNum;
    
    int useOnlyThisHex = -1;

    
public:
    complexHexMaskFPM();
    complexHexMaskFPM(initCommandSet*& cmdBlock, double zoomFactor);
    
    void init(initCommandSet*& cmdBlock);
    void set(std::string fieldName, const char *arg);
    void init_mask(void);
    void load_sags(const char *filename);
    double compute_fpmScale(int);
    
    arma::cx_mat make_complex_intpolated_mask(double lambda, double time);
    void execute(efield* E, arma::cx_mat& tcmat, int sl, double time);
    
    void get_optimization_data(const char *dataName, void *data);
    void set_optimization_data(const char *dataName, void *data);

    void compute_hex_centers(void);
    void make_hex_array(void);
    void make_interpolation_data(void);
    
    double get_fpmScale(void) { return fpmScale; }
    void set_zoomFactor(double z) {  zoomFactor = z; }
    
    void print(const char *hdr = "");
    void draw(const char *title = "");
};

#endif /* binary_complexHexMaskFPM_hpp */
