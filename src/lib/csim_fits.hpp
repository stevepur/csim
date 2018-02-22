//
//  csim_fits.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef csim_fits_hpp
#define csim_fits_hpp

#include "fitsio.h"
#include "armadillo"

int read_fits_array(const char *fn, void **data, int& ndims, long **dims, int& dataType, int verbose = 0);
void write_fits_array(const char *fn, void *data, int ndims, long *dims, int dataType, int verbose = 0);
int load_cube(const char *filename, arma::cube& fillCube, int verbose = 0);
int load_cube(const char *filename, arma::icube& fillCube, int verbose = 0);
void save_cube(const char *filename, arma::cube& saveCube, int verbose = 0);
void save_cube(const char *filename, arma::icube& saveCube, int verbose = 0);
void load_mat(const char *filename, arma::mat& fillMat, int verbose = 0);
void load_mat(const char *filename, arma::cx_mat& fillMat, int verbose = 0);
void load_vec(const char *filename, arma::vec& fillVec, int verbose = 0);
void save_mat(const char *filename, arma::mat& saveMat, int verbose = 0);
void save_mat(const char *filename, arma::cx_mat& saveMat, const char *saveType = "reIm", int verbose = 0);

#endif /* csim_parser_hpp */
