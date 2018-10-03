//
//  regionContrast.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef regionContrast_hpp
#define regionContrast_hpp

#include "../data/efield.hpp"

class regionContrast {
    double referenceLambda = -1.0;
    double loD = 0;
    double xlim1 = -1e16;
    double xlim2 = 1e16;
    double radius1 = 1.0;
    double radius2 = 2.0;
    double angle1 = 0.0;
    double angle2 = 30.0;


    char *filename = NULL;
    bool draw = false;
    
public:

    regionContrast();
    regionContrast(initCommandSet*& cmdBlocks);
    
    void init(initCommandSet*& cmdBlocks);
    void set(std::string fieldName, const char *arg);
    void set_loD(void);
    
    void execute(efield *inE, double time);
    void get_region_pixels(efield *E, arma::uvec& pixelIndex);
    void get_region_pixels(efield *E, arma::uvec& pixelIndex, arma::vec& pixelX, arma::vec& pixelY);
    double compute_intensity(efield *E, arma::mat& intensitySum, int calibrationState);
    void compute_contrast(void);
    void print(const char *hdr = NULL);
};


#endif /* regionContrast_hpp */
