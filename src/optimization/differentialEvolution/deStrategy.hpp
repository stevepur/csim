//
//  deStrategy_hpp.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef deStrategy_hpp
#define deStrategy_hpp

#include "armadillo"

class deStrategy {
    
public:
    const char *name;

    deStrategy() {}

    virtual void execute(arma::mat& population, arma::mat& candidatePopulation, arma::mat& originalPopulation, arma::vec& bestVector, double fWeight) {}
    void print(const char *hdr = NULL);
    void makePopulationPermutations(arma::mat& population, arma::cube& popPermutations);
};

class deRand1 : public deStrategy {
    
public:
    
    deRand1();

    void execute(arma::mat& population, arma::mat& candidatePopulation, arma::mat& originalPopulation, arma::vec& bestVector, double fWeight);
};

class deLocalToBest1 : public deStrategy {
    
public:
    
    deLocalToBest1();
    
    void execute(arma::mat& population, arma::mat& candidatePopulation, arma::mat& originalPopulation, arma::vec& bestVector, double fWeight);
};

class eitherOr : public deStrategy {
    
public:
    
    eitherOr();
    
    void execute(arma::mat& population, arma::mat& candidatePopulation, arma::mat& originalPopulation, arma::vec& bestVector, double fWeight);
};

#endif /* deStrategy_hpp */
