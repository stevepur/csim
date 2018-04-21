//
//  differentialEvolutionOptimizer.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef differentialEvolutionOptimizer_hpp
#define differentialEvolutionOptimizer_hpp

#include "armadillo"
#include "../optObjective.hpp"
#include "deStrategy.hpp"

#define DE_INIT_RANDU 0
#define DE_INIT_RANDN 1

class differentialEvolutionOptimizer {
    arma::mat population;
    arma::mat originalPopulation;
    arma::mat mutatedPopulation;
    arma::mat candidatePopulation;
    
    arma::vec populationMean;
    arma::vec popDistanceFromMean;
    double popMeanDistanceFromMean = 0.0;
    double popMaxDistanceFromMean = 0.0;
    double bestVal = 1e15;
    arma::vec bestVector;
    
    double globalBestVal = 1e15;
    arma::vec globalBestVector;
    
    arma::vec globalLowerBound;
    arma::vec globalUpperBound;
    arma::vec initLowerBound;
    arma::vec initUpperBound;
    int initType = DE_INIT_RANDU;
    double initSpread = 0.1;
    double reinitCompression = 1.0;
    double reinitSpread = 0.0;
    
    int populationSize = 0;
    int optVecSize = 0;
    int iteration = 0;
    
    int maxIterations = 1000000;
    double stopObjectiveValue = 0.0;
    double stopPopulationSpread = 0.0;
    
    bool doReinit = false;
    double reinitChangefracThreshold = 0.0;
    int maxSmallChange = 5;
    int maxRetries = 5;

    double fWeight = 0.8;
    double crossoverProb = 1.0;
    
    optObjective *objective = NULL;
    arma::rowvec objectiveValue;
    arma::rowvec candidateObjectiveValue;
    
    deStrategy *strategy = NULL;
    deStrategy *startStrategy = NULL;
    deStrategy *baseStrategy = NULL;
    deStrategy *endStrategy = NULL;
    deStrategy *reinitStrategy = NULL;
    double baseStrategySpreadSwitch = 0.0;
    double endStrategySpreadSwitch = 0.0;
    
    int outputInterval = 100;
    char *outputDirectory = NULL;
    bool saveState = true;
    
    FILE *bestValFid = NULL;
    FILE *histValFid = NULL;
    
    deStrategy *make_strategy(const char *arg);
    bool stop_criterion(void);
    void init_population(arma::vec initOrigin, int randomType = DE_INIT_RANDU, double spread = 1.0);
    void find_best_member(void);
    void do_crossover(void);
    void enforce_boundaries(void);
    void update_population(void);
    void compute_population_statistics(void);
    void display_population_statistics(arma::mat& p);
    
public:

    differentialEvolutionOptimizer();
    differentialEvolutionOptimizer(initCommandSet*& cmdBlocks);
    
    void init(initCommandSet*& cmdBlocks);
    void init_block(initCommandSet*& cmdBlock);
    void set(std::string fieldName, const char *arg);
    
    void optimize(void);
    void print(const char *hdr = NULL);
};


#endif /* differentialEvolutionOptimizer_hpp */
