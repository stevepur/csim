//
//  linearOptimizer.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef linearOptimizer_hpp
#define linearOptimizer_hpp

#include "../lib/csim_lib.hpp"
#include "armadillo"
#include "nlopt.hpp"
#include "../tools/regionContrast.hpp"
#include "../data/efield.hpp"

class linearOptimizer {

    FILE *bestValFid = NULL;
    FILE *histValFid = NULL;
    
    regionContrast *region = NULL;

    efield *calibEfield = NULL;
    efield *fullEfield = NULL;
    efield *jacEfield = NULL;
    double calibMaxIntensity = 0;
    arma::uvec regionPixelIndex;
    
    arma::vec optVec;
    arma::vec deltaOptVec;
    arma::vec lastDeltaOptVec;
    arma::vec startVec;
    arma::vec finalVec;
    arma::mat objMat;
    arma::mat lastObjMat;
    arma::vec objVec;
    
    double contrast = 0;
    double lastObjMatContrast = 0;
    double lastContrast = 0;
    double jacobianDx = 1e-8;
    double regularizationMu = 0.1;

    char *componentName = NULL;
    char *dataName = NULL;

    double globalLowerBound;
    double globalUpperBound;
    int optVecSize = 0;
    double stepSizeScale = 0.1;
    int nJacRows = 0;
    
    arma::wall_clock optimizationTimer;

    int nIterations = 0;
    int maxIterations = 20;
    double bestVal = 1e15;
    double stopContrastValue = 1e-9;
    double contrastRelTolerance = 1e-5;
    double deltaOptVecTolerance = 1e-9;
    bool stop_criterion(void);
    
    bool optVecNotChanging = false;
    bool contrastNotChanging = false;
    bool contrastReachedStopValue = false;

    int outputInterval = 100;
    char *outputDirectory = NULL;
    bool saveState = true;

public:

    linearOptimizer();
    linearOptimizer(initCommandSet*& cmdBlocks);
    
    static double eval_mu_contrast(const std::vector<double> &x, std::vector<double> &grad, void *parentPointer);

    void init(initCommandSet*& cmdBlocks);
    void init_block(initCommandSet*& cmdBlocks);
    void set(std::string fieldName, const char *arg);
    
    void optimize(void);
    bool compute_step(void);
    void compute_contrast(void);
    void compute_jacobian(void);
    void compute_objective_vector(void);
    void reset_Efield(efield *& E);
    bool stopping_criteria(void);

    FILE *get_bestValFid(void) { return bestValFid; }
    FILE *get_histValFid(void) { return histValFid; }
    int get_nIterations(void) { return nIterations; }
    double get_bestVal(void) { return bestVal; }
    double set_bestVal(double val) { bestVal = val; }
    char *get_ouput_directory(void) { return outputDirectory; }
    void print(const char *hdr = NULL);
};


#endif /* linearOptimizer_hpp */
