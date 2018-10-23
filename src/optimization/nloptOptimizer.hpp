//
//  nloptOptimizer.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef nloptOptimizer_hpp
#define nloptOptimizer_hpp

#include "../lib/csim_lib.hpp"
#include "armadillo"
#include "nlopt.hpp"
#include "../tools/regionContrast.hpp"
#include "../data/efield.hpp"

typedef struct {
    double a, b;
} my_constraint_data;

class nloptOptimizer {

    FILE *bestValFid = NULL;
    FILE *histValFid = NULL;
    
    regionContrast *region = NULL;

    efield *calibEfield = NULL;
    efield *fullEfield = NULL;
    double calibMaxIntensity = 0;
    arma::uvec regionPixelIndex;
    
    arma::vec initX;
    arma::vec startVec;
    arma::vec finalVec;
    arma::vec previousGradient;
    arma::vec previousGradientX;
    double previousGradientContrast = 0;
    double lastContrast = 0;
    double gradientDx = 1e-8;
    double gradientErrorThreshold = 1e-8;

    char *componentName = NULL;
    char *dataName = NULL;

    double globalLowerBound;
    double globalUpperBound;
    double initLowerBound;
    double initUpperBound;
    int optVecSize = 0;

    arma::wall_clock optimizationTimer;

    nlopt::algorithm optMethod = nlopt::LN_NELDERMEAD;
    int nIterations = 0;
    int maxIterations = 1000000;
    double bestVal = 1e15;
    double stopObjectiveValue = 1e-9;
    double xRelTolerance = 1e-5;
    double fRelTolerance = 1e-5;
    bool stop_criterion(void);

    int outputInterval = 100;
    char *outputDirectory = NULL;
    bool saveState = true;

public:

    nloptOptimizer();
    nloptOptimizer(initCommandSet*& cmdBlocks);
    
    static double eval_contrast(const std::vector<double> &x, std::vector<double> &grad, void *parentPointer);

    void init(initCommandSet*& cmdBlocks);
    void init_block(initCommandSet*& cmdBlocks);
    void set(std::string fieldName, const char *arg);
    
    void optimize(void);
    void reset_fullEfield(void);
    
    FILE *get_bestValFid(void) { return bestValFid; }
    FILE *get_histValFid(void) { return histValFid; }
    void increment_nIterations(void) { nIterations++; }
    int get_nIterations(void) { return nIterations; }
    double get_bestVal(void) { return bestVal; }
    double set_bestVal(double val) { bestVal = val; }
    char *get_ouput_directory(void) { return outputDirectory; }
    void print(const char *hdr = NULL);
};


#endif /* nloptOptimizer_hpp */
