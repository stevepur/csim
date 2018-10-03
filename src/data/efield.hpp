//
//  efield.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef efield_hpp
#define efield_hpp

#include "armadillo"
#include "../lib/csim_parser.hpp"
#include "../lib/csim_geom.hpp"
#undef ARMA_BLAS_UNDERSCORE

#define DEFAULT_FIELD_SIZE_X 1024
#define DEFAULT_FIELD_SIZE_Y 1024
#define DEFAULT_N_LAMBDA 1
#define DEFAULT_CONSTANT_VAL 1

class lambdaDataClass {
    public:
    double lambda = 0;
    double focalLengthLambdaOverD = 0;
    
    lambdaDataClass();
    lambdaDataClass(double wavelength);
    
    void set_wavelength(double wavelength);
    double get_wavelength(void);
};

class pointSourceClass {
    public:
    double flux = 1;
    double tipX = 0;
    double tiltY = 0;
    
    pointSourceClass();
    pointSourceClass(double f=1.0, double tx=0.0, double ty=0.0);
    void print(const char *hdr = "");
};

typedef std::vector<arma::cx_cube *> polarizations;

class efield {
public:
    char *name = NULL;
    int verbose;
    char *outputDirectory = NULL;
    std::vector<polarizations> E;
    std::vector<pointSourceClass> pointSourceList;
    double referenceWavelength = 0;
    double wavelengthStart = 0;
    double wavelengthEnd = 0;
    double deltaWavelength = 0;
    int nInitedSources = 0;
    
    lambdaDataClass *lambdaData;
    
    ////////////////////////////////////////////////////////
    // optical parameters that can change as the efield moves through the coronagraph
    //
    // for convenience, call physical measurement units "m"
    // call pixel units "p"
    //
    // beamRadius (beamrad in matlab code) is the physical radius of the beam in m
    double beamRadiusPhysical = 0;
    // pixelScale (pscale in matlab code) is the physical size of a pixel = m/p
    double pixelScale = 0;
    // size of the optical array (N in matlab code).  Must match the size of the E field
    int arraySizePixels = 0;
    // beam radius in pixels (Nbeam in matlab code)
    double beamSizePixels = 0;
    // physical half-size of the array in m (prad in matlab code)
    double arrayRadiusPhysical = 0;
    ////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////
    // array geometry that can change as the efield moves through the coronagraph
    //
    //
    arrayGeom arrayGeometry;
    arrayGeom initArrayGeometry;
    ////////////////////////////////////////////////////////
    
    
    efield();
    efield(efield& in);
    efield(char *reFilename, char *imFilename);
    efield(int nRows, int nColumns, int nLambda = DEFAULT_N_LAMBDA);
    efield(initCommandSet*& cmdBlock);
    ~efield(void);

    efield& operator=(const efield& in);

    void init(char *reFilename, char *imFilename);
    void init(double initValue = DEFAULT_CONSTANT_VAL);
    void init(initCommandSet*& cmdBlock);
    
    void set(std::string fieldName, const char *arg);
    void add_point_source(double flux=1.0, double tipX=0.0, double tiltY=0.0);
    void init_point_sources(void);
    std::complex<double> getV(int r, int c, int lambda = DEFAULT_N_LAMBDA, int p = 0, int s = 0);
    void setV(std::complex<double> val, int r, int c, int lambda = DEFAULT_N_LAMBDA, int p = 0, int s = 0);
    void set_size(int nRows, int nColumns, int nLambda = DEFAULT_N_LAMBDA);
    void set_wavelength_range(const char *arg);
    void set_wavelength_bandwidth(const char *arg);
    void set_wavelength_list(const char *arg);
    void set_optical_parameters(void);
    void set_array_geometry(void);
    void save(const char *coreName);
    void draw(const char *title, const char *drawType = "logamp");
    void print(const char *header = "");
};

extern efield *initialEfield;

#endif /* efield_hpp */
