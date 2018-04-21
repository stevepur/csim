//
//  optObjective_hpp.cpp
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
#include "deStrategy.hpp"

void deStrategy::print(const char *hdr) {
    if (hdr != NULL)
        std::cout << hdr << ":" << std::endl;
    std::cout << name << std::endl;
}

void deStrategy::makePopulationPermutations(arma::mat& population, arma::cube& popPermutations) {
    
    // make the offset vector for each
    arma::uvec ii2(popPermutations.n_slices);
    for (int i=0; i<ii2.n_rows; i++)
        ii2(i) = i;
    arma::uvec si2 = shuffle(ii2);
//    std::cout << "ii2" << ii2.t() << std::endl;
//    std::cout << "si2" << si2.t() << std::endl;
    
    
    // make the index vector that will be permuted
    arma::uvec ii1(population.n_cols);
    for (int i=0; i<ii1.n_rows; i++)
        ii1(i) = i;
    arma::uvec si1 = shuffle(ii1);
//    std::cout << "si1" << si1.t() << std::endl;
    
    
    // make the permuted populations via offsets into si1
    for (int i=0; i<popPermutations.n_slices; i++) {
//        std::cout << "shifted si1" << si1.t() << std::endl;
        popPermutations.slice(i) = population.cols(si1);
        if (i < popPermutations.n_slices - 1) // don't need to do the last time
            si1 = si1(shift(ii1, si2(i)+1));
    }
    
//    for (int i=0; i<popPermutations.n_slices; i++) {
//        arma::urowvec tv = si1(shift(ii1, si2(i))).t();
//        std::cout << "si1(shift(ii1, si2(i))): " << tv(arma::span(0,9)) << std::endl;
//    }
    if (0) {
        std::cout << "population:" << std::endl;
        std::cout << population(arma::span(0,9), arma::span(0,9)) << std::endl;
        
        for (int i=0; i<popPermutations.n_slices; i++) {
//            std::cout << "si1(shift(ii1, si2(i))): " << si1(shift(ii1, si2(i))).t() << std::endl;
            std::cout << "popPermutations slice " << i << ":" << std::endl;
            std::cout << popPermutations(arma::span(0,9), arma::span(0,9), arma::span(i)) << std::endl;
        }
    }

}

//////////////////////////////////////////////////////////////////

deRand1::deRand1() {
    name = "DE/Rand/1";
}

void deRand1::execute(arma::mat& population, arma::mat& mutatedPopulation, arma::mat& originalPopulation, arma::vec& bestVector, double fWeight) {
    
    arma::cube popPermutations(population.n_rows, population.n_cols, 3);
    makePopulationPermutations(population, popPermutations);
    
    mutatedPopulation = popPermutations.slice(2) + fWeight*(popPermutations.slice(0) - popPermutations.slice(1));
    if (0) {
        std::cout << "deRand1 population:" << std::endl;
        std::cout << population(arma::span(0,9), arma::span(0,9)) << std::endl;

        arma::mat tmat;
        tmat = popPermutations.slice(0);
        std::cout << "deRand1 permuted mat 0:" << std::endl;
        std::cout << tmat(arma::span(0,9), arma::span(0,9)) << std::endl;
        tmat = popPermutations.slice(1);
        std::cout << "deRand1 permuted mat 1:" << std::endl;
        std::cout << tmat(arma::span(0,9), arma::span(0,9)) << std::endl;
        tmat = popPermutations.slice(2);
        std::cout << "deRand1 permuted mat 2:" << std::endl;
        std::cout << tmat(arma::span(0,9), arma::span(0,9)) << std::endl;
        
        std::cout << "deRand1 mutated population fWeight = " << fWeight << ":" << std::endl;
        std::cout << mutatedPopulation(arma::span(0,9), arma::span(0,9)) << std::endl;
    }
    
    originalPopulation = popPermutations.slice(2);
}

//////////////////////////////////////////////////////////////////

deLocalToBest1::deLocalToBest1() {
    name = "DE/local-to-best/1";
}

void deLocalToBest1::execute(arma::mat& population, arma::mat& mutatedPopulation, arma::mat& originalPopulation, arma::vec& bestVector, double fWeight) {
    
    arma::cube popPermutations(population.n_rows, population.n_cols, 2);
    makePopulationPermutations(population, popPermutations);

    mutatedPopulation = population + fWeight*(repmat(bestVector, 1, population.n_cols) - population) + fWeight*(popPermutations.slice(0) - popPermutations.slice(1));
    
    originalPopulation = population;
//    arma::cube popPermutations(2, population.n_elem);
//    
//    mutatedPopulation = population + fWeight*(repmat(bestVector, population.n_rows, 1) - population) + fWeight*(popPermutations(0,arma::span::all) - popPermutations(1,arma::span::all));
//    
//    originalPopulation = population;
}

//////////////////////////////////////////////////////////////////

eitherOr::eitherOr() {
    name = "DE/Either-Or";
}

void eitherOr::execute(arma::mat& population, arma::mat& mutatedPopulation, arma::mat& originalPopulation, arma::vec& bestVector, double fWeight) {
    
    arma::cube popPermutations(population.n_rows, population.n_cols, 3);
    makePopulationPermutations(population, popPermutations);
    
    if (all(all(popPermutations.slice(0) == popPermutations.slice(1))))
        std::cout << "population permutation 0 == population permutation 1!!!" << std::endl;
    if (all(all(popPermutations.slice(1) == popPermutations.slice(2))))
        std::cout << "population permutation 1 == population permutation 2!!!" << std::endl;
    if (all(all(popPermutations.slice(0) == popPermutations.slice(2))))
        std::cout << "population permutation 0 == population permutation 2!!!" << std::endl;
    

    if ((double)rand() / ((double)RAND_MAX) < 0.5) {
        mutatedPopulation = popPermutations.slice(2) + fWeight*(popPermutations.slice(0) - popPermutations.slice(1));
    }
    else {
        mutatedPopulation = popPermutations.slice(2) + 0.5*(fWeight + 1.0)*(popPermutations.slice(0) + popPermutations.slice(1) - 2*popPermutations.slice(2));
    }
    
    originalPopulation = popPermutations.slice(2);
//    arma::mat popPermutations(3, population.n_elem);
//    makePopulationPermutations(population, popPermutations);
    
//    if ((double)rand() / ((double)RAND_MAX) < 0.5)
//        mutatedPopulation = popPermutations(2,arma::span::all) + fWeight*(popPermutations(0,arma::span::all) - popPermutations(1,arma::span::all));
//    else
//        mutatedPopulation = popPermutations(2,arma::span::all) + 0.5*(fWeight + 1.0)*(popPermutations(0,arma::span::all) + popPermutations(1,arma::span::all) - 2*popPermutations(2,arma::span::all));
//    
//    originalPopulation = popPermutations(2,arma::span::all);
}


