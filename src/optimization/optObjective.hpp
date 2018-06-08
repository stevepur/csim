//
//  optObjective_hpp.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef optObjective_hpp
#define optObjective_hpp

#include "../data/responseData.hpp"

#define OBJECTIVE_MEAN 0
#define OBJECTIVE_MAX 1

class optObjective {
    
    
public:
    int objectiveType = OBJECTIVE_MEAN;

    optObjective() {}
    optObjective(initCommandSet*& cmdBlocks) {}
    
    virtual void init(initCommandSet*& cmdBlock) {}
    virtual void set(std::string fieldName, const char *arg) {}
    
    virtual void execute(arma::mat& optData, arma::rowvec& objVal) {}
    
    virtual int get_opt_vector_size(void) {}
    virtual void print(const char *hdr = NULL) {}
};

class optHexRespObjective : public optObjective {
    responseData *response = NULL;
    arma::cx_mat rotatedPhases;
    double phaseSign = 1.0;
    
public:
    
    optHexRespObjective();
    optHexRespObjective(initCommandSet*& cmdBlocks);
    
    void init(initCommandSet*& cmdBlock);
    void set(std::string fieldName, const char *arg);
    
    void execute(arma::mat& optData, arma::rowvec& objVal);
    int get_opt_vector_size(void);
    void print(const char *hdr = NULL);
};

#endif /* optObjective_hpp */
