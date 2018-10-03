//
//  csim_utils_cpp
//  csim
//
//  Created by steve on 4/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include <iostream>
#include <assert.h>

#include "csim_utils.hpp"

void arma_vec_to_std_vec(arma::vec& vIn, std::vector<double>& vOut) {
    vOut.resize(vIn.n_elem);
    for (int i=0; i<vIn.n_elem; i++)
        vOut[i] = vIn[i];
}

void std_vec_to_arma_vec(const std::vector<double>& vIn, arma::vec& vOut) {
    vOut.set_size(vIn.size());
    for (int i=0; i<vOut.n_elem; i++)
        vOut[i] = vIn[i];
}

