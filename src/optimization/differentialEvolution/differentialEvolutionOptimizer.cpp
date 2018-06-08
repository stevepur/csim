//
//  differentialEvolutionOptimizer.cpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>

#include "armadillo"
#include "differentialEvolutionOptimizer.hpp"
#include "../lib/csim_lib.hpp"

differentialEvolutionOptimizer::differentialEvolutionOptimizer() {
}

differentialEvolutionOptimizer::differentialEvolutionOptimizer(initCommandSet*& cmdBlocks) {
    init(cmdBlocks);
}

void differentialEvolutionOptimizer::init(initCommandSet*& cmdBlock) {
    
    std::vector<initCommandSet*> subBlocks = cmdBlock->find_command_blocks();
    
    for (int j=subBlocks.size()-1; j>=0; j--) {
        if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "optHexRespObjective")) {
            objective = new optHexRespObjective(subBlocks[j]);
            
            optVecSize = objective->get_opt_vector_size() - 1; // the response matrix has an extra column that is not optimized
        }
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "differentialEvolutionOptimizer")) {
            init_block(subBlocks[j]);
        }
    }
    if (startStrategy != NULL)
        strategy = startStrategy;
    else if (baseStrategy != NULL)
        strategy = baseStrategy;
    else if (endStrategy != NULL)
        strategy = endStrategy;
    
    if (reinitStrategy == NULL)
        reinitStrategy = baseStrategy;
    
    popDistanceFromMean.set_size(populationSize);
    
    print("initial differential evoluation optimizer");

    assert(strategy);
    
}


void differentialEvolutionOptimizer::init_block(initCommandSet*& cmdBlock) {

    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
}

void differentialEvolutionOptimizer::set(std::string fieldName, const char *arg) {
    if (fieldName == "differentialEvolutionOptimizer")
        ;
    // set strategies
    else if (fieldName == "startStrategy") {
        startStrategy = make_strategy(arg);
    }
    else if (fieldName == "baseStrategy") {
        baseStrategy = make_strategy(arg);
    }
    else if (fieldName == "endStrategy") {
        endStrategy = make_strategy(arg);
    }
    else if (fieldName == "reinitStrategy") {
        reinitStrategy = make_strategy(arg);
    }
    // initialize population
    else if (fieldName == "populationSize") {
        // arg is a single int
        populationSize = atoi(arg);
    }
    else if (fieldName == "initPopulationMean") {
        // arg is a string
        if (!strcmp(arg, "zero")) {
            populationMean = arma::zeros<arma::vec>(optVecSize);
        }
        else {
            char *filename = new char[strlen(arg)+1];
            strcpy(filename, arg);
            load_vec(filename, populationMean);
        }
    }
    else if (fieldName == "initSpread") {
        // arg is a single float
        initSpread = atof(arg);
    }
    else if (fieldName == "reinitCompression") {
        // arg is a single float
        reinitCompression = atof(arg);
    }
    else if (fieldName == "reinitSpread") {
        // arg is a single float
        reinitSpread = atof(arg);
    }
    // set bounds
    else if (fieldName == "globalBounds") {
        float v1, v2;
        // arg is a comma separated pair of doubles
        sscanf(arg, "%f, %f", &v1, &v2);
        globalLowerBound = v1*arma::ones<arma::vec>(optVecSize);
        globalUpperBound = v2*arma::ones<arma::vec>(optVecSize);
    }
    else if (fieldName == "initBounds") {
        float v1, v2;
        // arg is a comma separated pair of doubles
        sscanf(arg, "%f, %f", &v1, &v2);
        initLowerBound = v1*arma::ones<arma::vec>(optVecSize);
        initUpperBound = v2*arma::ones<arma::vec>(optVecSize);
    }
    else if (fieldName == "initType") {
        // arg is a string
        if (!strcmp(arg, "uniformRandom")) {
            initType = DE_INIT_RANDU;
        }
        else if (!strcmp(arg, "normalRandom")) {
            initType = DE_INIT_RANDN;
        }
        else
            std::cout << "!!! differentialEvolutionOptimizer initType bad arg: " << arg << std::endl;
    }
    // set stopping criteria
    else if (fieldName == "maxIterations") {
        // arg is a single int
        maxIterations = atoi(arg);
    }
    else if (fieldName == "stopObjectiveValue") {
        // arg is a single float
        stopObjectiveValue = atof(arg);
    }
    else if (fieldName == "stopPopulationSpread") {
        // arg is a single float
        stopPopulationSpread = atof(arg);
    }
    // set retry/reinit criteria
    else if (fieldName == "reinitChangefracThreshold") {
        // arg is a single float
        reinitChangefracThreshold = atof(arg);
    }
    else if (fieldName == "maxSmallChange") {
        // arg is a single int
        maxSmallChange = atoi(arg);
    }
    else if (fieldName == "maxRetries") {
        // arg is a single int
        maxRetries = atoi(arg);
    }
    // set differential evoluation parameters
    else if (fieldName == "fWeight") {
        // arg is a single float
        fWeight = atof(arg);
    }
    else if (fieldName == "crossoverProb") {
        // arg is a single float
        crossoverProb = atof(arg);
    }
    // set strategy population spread swtich thresholds
    else if (fieldName == "baseStrategySpreadSwitch") {
        // arg is a single float
        baseStrategySpreadSwitch = atof(arg);
    }
    else if (fieldName == "endStrategySpreadSwitch") {
        // arg is a single float
        endStrategySpreadSwitch = atof(arg);
    }
    // set output parameters
    else if (fieldName == "outputInterval") {
        // arg is a single int
        outputInterval = atoi(arg);
    }
    else if (fieldName == "outputDirectory") {
        // arg is a string
        outputDirectory = new char[strlen(arg)+1];
        strcpy(outputDirectory, arg);
    }
    else if (fieldName == "saveState") {
        // arg is a single int
        saveState = atoi(arg);
    }
    else
        std::cout << "!!! differentialEvolutionOptimizer bad set field name: " << fieldName << std::endl;
    
}

deStrategy *differentialEvolutionOptimizer::make_strategy(const char *arg) {
    if (!strcmp(arg, "deRand1")) {
        return new deRand1;
    }
    else if (!strcmp(arg, "deLocalToBest1")) {
        return new deLocalToBest1;
    }
    else if (!strcmp(arg, "eitherOr")) {
        return new eitherOr;
    }
    else {
        std::cout << "!!! make_strategy bad arg: " << arg << std::endl;
        assert(NULL);
    }
    return NULL;
}

void differentialEvolutionOptimizer::optimize(void) {
    arma::wall_clock timer;
    int retryCount = 0;
    int smallChangeCount = 0;
    double preiousBestVal = 0;
    bool selectedBaseStrategy = false;
    bool selectedEndStrategy = false;
    
    std::cout << "============= starting optimization =============" << std::endl;
    std::string fname = (std::string)outputDirectory + "/optBestVal.txt";
    bestValFid = fopen(fname.c_str(), "w");
    
    fname = (std::string)outputDirectory + "/optValHistory.txt";
    histValFid = fopen(fname.c_str(), "w");
    
    assert(bestValFid);

    init_population(populationMean, initType);
    
    objective->execute(population, objectiveValue);
    
    find_best_member();
    compute_population_statistics();
    
    if (saveState) {
        char outStr[1000];
        sprintf(outStr, "%s/initial_population.fits", outputDirectory, iteration);
        save_mat(outStr, population);
    }
    timer.tic();
    arma::mat previousPopulation = population;
    while((!stop_criterion() && retryCount <= maxRetries) || iteration == 0) {
        
        if (smallChangeCount >= maxSmallChange || doReinit || popMeanDistanceFromMean < reinitSpread) {
            std::cout << "reininting" << std::endl;
            bestVal = 1e15;
            strategy = reinitStrategy;
            init_population(globalBestVector, DE_INIT_RANDU, reinitCompression);
            objective->execute(population, objectiveValue);
            find_best_member();
            smallChangeCount = 0;
            doReinit = false;
            retryCount++;
        }
    
        strategy->execute(population, mutatedPopulation, originalPopulation, bestVector, fWeight);
//        std::cout << "after execution population stats" << std::endl;
//        display_population_statistics(population);
//        std::cout << "after execution mutated population stats" << std::endl;
//        display_population_statistics(mutatedPopulation);
        
        do_crossover();
//        std::cout << "after crossover candidatePopulation stats" << std::endl;
//        display_population_statistics(candidatePopulation);

        enforce_boundaries();
//        std::cout << "after boundaries candidatePopulation stats" << std::endl;
//        display_population_statistics(candidatePopulation);
        
        update_population();
//        std::cout << "after update population stats" << std::endl;
//        display_population_statistics(population);
        compute_population_statistics();
        

        find_best_member();
        
        if (iteration % outputInterval == 0 || iteration == 0) {
            compute_population_statistics();

            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << "iteration " << iteration << ": bestVal = " << bestVal << ", smallChangeCount = " << smallChangeCount << ", time " << timer.toc() << " seconds, " << strategy->name << std::endl;
            std::cout << "popMeanDistanceFromMean = " << popMeanDistanceFromMean << ", popMaxDistanceFromMean = " << popMaxDistanceFromMean << std::endl;
            
            if (popMeanDistanceFromMean < baseStrategySpreadSwitch & !selectedBaseStrategy) {
                strategy = baseStrategy;
                selectedBaseStrategy = true;
                doReinit = true;
            }
            if (popMeanDistanceFromMean < endStrategySpreadSwitch & !selectedEndStrategy) {
                strategy = endStrategy;
                selectedEndStrategy = true;
                doReinit = true;
            }
            
            if (fabs((bestVal - preiousBestVal)/bestVal) < reinitChangefracThreshold)
                smallChangeCount++;
            else
                smallChangeCount = 0;
            preiousBestVal = bestVal;
            
            if (all(all(population == previousPopulation)))
                std::cout << "population unchanged" << std::endl;
            previousPopulation = population;
            
            fprintf(histValFid, "%d %g\n", iteration, bestVal);
            fflush(histValFid);
            
            if (saveState) {
                char outStr[1000];
                sprintf(outStr, "%s/population_it_%d.fits", outputDirectory, iteration);
                save_mat(outStr, population);
            }
            
            timer.tic();
        }
    
        iteration++;
    }
    fclose(bestValFid);
    fclose(histValFid);
    std::cout << "iteration = " << iteration << ", bestVal = " << bestVal << ", popMaxDistanceFromMean = " << popMaxDistanceFromMean << ", retryCount = " << retryCount << std::endl;
    
}

bool differentialEvolutionOptimizer::stop_criterion(void) {
    return ((iteration > maxIterations) || (bestVal < stopObjectiveValue) || (popMaxDistanceFromMean < stopPopulationSpread));
}

void differentialEvolutionOptimizer::init_population(arma::vec initOrigin, int randomType, double spread) {
    // population is optVecSize x populationSize
    // populationMean is as column vector, need to turn into row
    switch (randomType) {
        case DE_INIT_RANDU:
            population = repmat(initOrigin, 1, populationSize) + spread*initSpread*(arma::randu(optVecSize, populationSize) - 0.5) % repmat((initUpperBound - initLowerBound), 1, populationSize);
            break;
            
        case DE_INIT_RANDN:
            population = repmat(initOrigin, 1, populationSize) + spread*initSpread*arma::randn(optVecSize, populationSize) % repmat((initUpperBound - initLowerBound), 1, populationSize);
            break;
            
        default:
            assert(NULL);
    }
    if (0) {
        arma::mat im = repmat(initOrigin, 1, populationSize);
        std::cout << "input sags:" << std::endl;
        std::cout << im(arma::span(0,9), arma::span(0,9)) << std::endl;
        std::cout << "init population:" << std::endl;
        std::cout << population(arma::span(0,9), arma::span(0,9)) << std::endl;
    }
}

void differentialEvolutionOptimizer::do_crossover(void) {
    if (crossoverProb < 1.0) {
        arma::umat crossTrue = arma::randu(size(population)) < crossoverProb;
        arma::umat crossFalse = crossTrue < 0.5;
        
        candidatePopulation = population % crossFalse + mutatedPopulation % crossTrue;
    }
    else
        candidatePopulation = mutatedPopulation;
}

void differentialEvolutionOptimizer::find_best_member(void) {
    // population is optVecSize x populationSize
    // populationMean is as column vector, need to turn into row
//    std::cout << "starting find_best_member" << std::endl;
//    std::cout << "mean(objectiveValue) = " << mean(objectiveValue) << std::endl;
    int bestIndex = objectiveValue.index_min();
    if (objectiveValue(bestIndex) < bestVal) {
        bestVal = objectiveValue(bestIndex);
        bestVector = population(arma::span::all, bestIndex);
        
        if (bestVal < globalBestVal) {
            globalBestVal = bestVal;
            globalBestVector = bestVector;
            
            std::string fname = (std::string)outputDirectory + "/bestOptVector.fits";
            save_vec(fname.c_str(), globalBestVector);
            
            fprintf(bestValFid, "%d %g\n", iteration, globalBestVal);
            fflush(bestValFid);
        }
    }
//    std::cout << "completed find_best_member" << std::endl;
}

void differentialEvolutionOptimizer::enforce_boundaries(void) {
    // population is optVecSize x populationSize
    // populationMean is as column vector, need to turn into row
    
//    std::cout << "before enforce_boundaries: max(abs(pop)) = " << max(max(population)) << std::endl;
//    std::cout << "mean(populationMean) = " << mean(populationMean) << ", max(populationMean) = " << max(populationMean) << ", min(populationMean) = " << min(populationMean) << std::endl;

    arma::vec upperBound = populationMean + initUpperBound;
    arma::vec lowerBound = populationMean + initLowerBound;
    
    #pragma omp parallel for
    for (int i=0; i<upperBound.n_elem; i++) {
        if (upperBound(i) > globalUpperBound(i))
            upperBound(i) = globalUpperBound(i);
        if (lowerBound(i) < globalLowerBound(i))
            lowerBound(i) = globalLowerBound(i);
    }
//    std::cout << "upperBound Range: " << min(upperBound) << " to " << max(upperBound) << ", lowerBound Range: " << min(lowerBound) << " to " << max(lowerBound) << std::endl;

    // enforce boundaries with a bounceback
    #pragma omp parallel for
    for (int c=0; c<population.n_cols; c++) {
        for (int r=0; r<population.n_rows; r++) {
            if (candidatePopulation(r,c) < lowerBound(r)) {
                double rnd = (double)rand() / ((double)RAND_MAX);
                candidatePopulation(r,c) = lowerBound(r) + rnd*(originalPopulation(r,c) - lowerBound(r));
            }
            else if (candidatePopulation(r,c) > upperBound(r)) {
                double rnd = (double)rand() / ((double)RAND_MAX);
                candidatePopulation(r,c) = upperBound(r) + rnd*(originalPopulation(r,c) - upperBound(r));
            }
        }
    }
    if (0) {
        std::cout << "0: lowerBound = " << lowerBound(0) << ", upperBound = " << upperBound(0) << std::endl;
        std::cout << "after enforce_boundaries: 0: " << min(candidatePopulation(0, arma::span::all)) << " to " << max(candidatePopulation(0, arma::span::all)) << std::endl;
        std::cout << "1: lowerBound = " << lowerBound(1) << ", upperBound = " << upperBound(1) << std::endl;
        std::cout << "after enforce_boundaries: 1: " << min(candidatePopulation(1, arma::span::all)) << " to " << max(candidatePopulation(1, arma::span::all)) << std::endl;
    }
}

void differentialEvolutionOptimizer::update_population(void) {
    objective->execute(candidatePopulation, candidateObjectiveValue);
    
    bool foundBetter = false;
    for (int i=0; i<populationSize; i++) {
        if (candidateObjectiveValue(i) < objectiveValue(i)) {
            population(arma::span::all, i) = candidatePopulation(arma::span::all, i);
            objectiveValue(i) = candidateObjectiveValue(i);
            foundBetter = true;
        }
    }
//    if (!foundBetter) {
//        std::cout << "candidate population did not find a better value" << std::endl;
//    }
//    else
//        std::cout << "candidate population found a better value" << std::endl;

}


void differentialEvolutionOptimizer::compute_population_statistics(void) {
    //    std::cout << "starting compute_population_statistics" << std::endl;
    //    std::cout << "size of arma::mean(population, 1): " << size(arma::mean(population, 1)) << std::endl;
    populationMean = arma::mean(population, 1);
    //    std::cout << "populationMean: " << size(populationMean) << std::endl;
    for (int i=0; i<population.n_cols; i++) {
        popDistanceFromMean(i) = arma::norm(population(arma::span::all, i) - populationMean);
    }
    //    std::cout << "popDistanceFromMean: " << size(popDistanceFromMean) << std::endl;
    popMeanDistanceFromMean = arma::mean(popDistanceFromMean);
    popMaxDistanceFromMean = arma::max(popDistanceFromMean);
    //    std::cout << "completed compute_population_statistics" << std::endl;
}

void differentialEvolutionOptimizer::display_population_statistics(arma::mat& p) {
    arma::vec pMean = arma::mean(p, 1);

    arma::vec pDistanceFromMean(populationSize);
    #pragma omp parallel for
    for (int i=0; i<population.n_cols; i++) {
        pDistanceFromMean(i) = arma::norm(p(arma::span::all, i) - pMean);
    }
    std::cout << "mean distance from mean = " << arma::mean(pDistanceFromMean) << std::endl;
    std::cout << "max distance from mean = " << arma::max(pDistanceFromMean) << std::endl;
}

void differentialEvolutionOptimizer::print(const char *hdr) {
    if (hdr != NULL)
        std::cout << hdr << ":" << std::endl;
    std::cout << "globalBestVal: " << globalBestVal << std::endl;
    std::cout << "populationSize: " << populationSize << std::endl;
    std::cout << "optVecSize: " << optVecSize << std::endl;
    std::cout << "global bounds: " << arma::min(globalLowerBound) << " to " << arma::max(globalUpperBound) << std::endl;
    std::cout << "init bounds: " << arma::min(initLowerBound) << " to " << arma::max(initUpperBound) << std::endl;
    std::cout << "initType: " << initType << std::endl;
    std::cout << "initSpread: " << initSpread << std::endl;
    std::cout << "reinitCompression: " << reinitCompression << std::endl;
    std::cout << "reinitSpread: " << reinitSpread << std::endl;
    std::cout << "maxIterations: " << maxIterations << std::endl;
    std::cout << "stopObjectiveValue: " << stopObjectiveValue << std::endl;
    std::cout << "stopPopulationSpread: " << stopPopulationSpread << std::endl;
    std::cout << "reinitChangefracThreshold: " << reinitChangefracThreshold << std::endl;
    std::cout << "maxSmallChange: " << maxSmallChange << std::endl;
    std::cout << "maxRetries: " << maxRetries << std::endl;
    std::cout << "fWeight: " << fWeight << std::endl;
    std::cout << "crossoverProb: " << crossoverProb << std::endl;
    std::cout << "baseStrategySpreadSwitch: " << baseStrategySpreadSwitch << std::endl;
    std::cout << "endStrategySpreadSwitch: " << endStrategySpreadSwitch << std::endl;
    std::cout << "outputInterval: " << outputInterval << std::endl;
    std::cout << "outputDirectory: " << outputDirectory << std::endl;
    std::cout << "saveState: " << saveState << std::endl;
    fflush(stdout);
    objective->print(">>> objective");
    strategy->print(">>> strategy");
    startStrategy->print(">>> startStrategy");
    baseStrategy->print(">>> baseStrategy");
    endStrategy->print(">>> endStrategy");
}
