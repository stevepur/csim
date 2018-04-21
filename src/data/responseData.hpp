//
//  responseData.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef responseData_hpp
#define responseData_hpp

#include "armadillo"
#include "../lib/csim_parser.hpp"
#undef ARMA_BLAS_UNDERSCORE


typedef std::vector<arma::cx_cube *> polarizations;

class responseData {
public:
    char *name = NULL;
    int verbose;
    char *outputDirectory = NULL;
    std::vector<polarizations> M;
    arma::vec wavelengths;
    double calibIntensity = 1.0;
    
    responseData();
    responseData(int nRows, int nColumns, int nSlices, int nPolarizations, int nSources);
    responseData(responseData& in);
    responseData(char *dirName);
    responseData(initCommandSet*& cmdBlock);
    ~responseData(void);

    responseData& operator=(const responseData& in);

    void init(int nRows, int nColumns, int nSlices, int nPolarizations, int nSources);
    void init(initCommandSet*& cmdBlock);
    
    void set(std::string fieldName, const char *arg);
    void set_size(int nRows, int nColumns, int nSlices);
    void save(void);
    void load(void);
    void read_problem_params(void);
    void print(const char *header = NULL);
};

#endif /* responseData_hpp */
