//
//  deformableMirror.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef deformableMirror_hpp
#define deformableMirror_hpp

#include "celem.hpp"
#include "armadillo"
#include "../lib/csim_lib.hpp"

class deformableMirror : public celem {
    arma::mat mirrorMat;
    arma::mat influenceFunction;
    arma::mat actuatorMat;
    
    arma::mat actuatorUpsample;
    arma::cx_mat FTactuatorUpsample;
    arma::cx_mat FTinfluenceFunction;

    fft *actuatorUpsampleFft = NULL;
    fft *influenceFunctionFft = NULL;
    ifft *surfIfft = NULL;

    double mirrorSign = 1;
    double sigma = 1;
    double nActuatorRows = 0;
    double nActuatorCols = 0;
    
    double optr0 = 0;
    double optc0 = 0;
    double optStride = 1;
    double optNRows = 32;
    double optNCols = 32;

public:
    deformableMirror();
    deformableMirror(const char *inName);
    deformableMirror(initCommandSet*& cmdBlock);
    
    void init(initCommandSet*& cmdBlock);
    void read_mirror_file(const char *filename);
    void read_actuator_file(const char *filename);

    void set(std::string fieldName, const char *arg);
    void get_optimization_data(const char *dataName, void *data);
    void set_optimization_data(const char *dataName, void *data);
    void save_optimization_data(const char *dataName, char *outputDirectory = NULL);
    void compute_influence_function(void);
    void compute_deformable_mirror_surface(void);

    efield* execute(efield* E, celem* prev, celem* next, double time);
    
    void draw(const char *title = "");
};

#endif /* deformableMirror_hpp */
