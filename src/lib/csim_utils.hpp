//
//  csim_utils_hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef csim_utils_hpp
#define csim_utils_hpp

#include "armadillo"

void arma_vec_to_std_vec(arma::vec& vIn, std::vector<double>& vOut);
void std_vec_to_arma_vec(const std::vector<double>& vIn, arma::vec& vOut);


#endif /* csim_utils_hpp */
