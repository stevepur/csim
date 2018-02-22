//
//  csim_fits.cpp
//  csim
//
//  Created by steve on 4/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include <iostream>
#include <assert.h>

#include "csim_fits.hpp"

static int bitPixDefs[] = {BYTE_IMG, SHORT_IMG, LONG_IMG, LONGLONG_IMG, FLOAT_IMG, DOUBLE_IMG};
static int dataTypeDefs[] = {TBYTE, TSHORT, TINT, TLONG, TFLOAT, TDOUBLE};

int read_fits_array(const char *fn, void **data, int& ndims, long **dims, int& dataType, int verbose) {
    fitsfile *fptr;
    int status = 0;
    int retStat = 0;
    int bitpix;
    long fpixel[3] = {1, 1, 1};
    int anynul = 0;
    
    if (verbose)
        std::cout << "opening " << fn << std::endl;
    fits_open_file(&fptr, fn, READONLY, &status);
    if (verbose)
        std::cout << "fits open status: " << status << std::endl;
    fits_report_error(stderr, status);
    if (status == FILE_NOT_OPENED)
        return status;

    retStat = fits_get_img_dim(fptr, &ndims, &status);
    *dims = (long *) malloc(ndims*sizeof(long));
    if (verbose)
        std::cout << "dims = " << dims << std::endl;
    
    retStat = fits_get_img_param(fptr, 10, &bitpix, &ndims, *dims, &status);
    if (verbose)
        std::cout << "fits img type bitpix : " << bitpix << " ndims " << ndims << ", retStat " << retStat << " status " << status << std::endl;
    
    for (unsigned int i=1; i<sizeof(bitPixDefs)/sizeof(int); i++) {
        if (bitpix == bitPixDefs[i])  {
            dataType = dataTypeDefs[i];
            break;
        }
    }
    long nElements = 1;
    for (int i=0; i<ndims; i++) {
        if (verbose)
            std::cout << "axis " << i << ": " << (*dims)[i] << std::endl;
        nElements *= (*dims)[i];
    }
    fits_report_error(stderr, status);
    if (verbose)
        std::cout << "nElements = " << nElements << std::endl;
    switch (bitpix) {
        case LONGLONG_IMG:
            if (verbose)
                std::cout << "allocating " << nElements*sizeof(long) << " bytes" << std::endl;
            
            *data = (long *) malloc(nElements*sizeof(long));
            break;
            
        case FLOAT_IMG:
            if (verbose)
                std::cout << "allocating " << nElements*sizeof(float) << " bytes" << std::endl;
            
            *data = (float *) malloc(nElements*sizeof(float));
            break;
            
        case DOUBLE_IMG:
            if (verbose)
                std::cout << "allocating " << nElements*sizeof(double) << " bytes" << std::endl;
            
            *data = (double *) malloc(nElements*sizeof(double));
            break;
    }
    
    fits_read_pix(fptr, dataType, fpixel, nElements, NULL, *data, &anynul, &status);
    if (verbose)
        std::cout << "fits_read_pix status " << status << std::endl;
    fits_report_error(stderr, status);
    
    fits_close_file(fptr, &status);
    
    assert(data);
    return 0;
}

void write_fits_array(const char *fn, void *data, int ndims, long *dims, int dataType, int verbose) {
    fitsfile *fptr;
    int status = 0;
    int retStat = 0;
    int bitpix;
    long *fpixel;
    char filename[200];
    
    
    if (fn[0] != '!') {
        strcpy(filename, "!");
        strcat(filename, fn);
    } else
        strcpy(filename, fn);
    for (unsigned int i=0; i<strlen(filename); i++) {
        if (filename[i] == ' ')
            filename[i] = '_';
    }
    if (verbose)
        std::cout << "creating " << filename << std::endl;
    fits_create_file(&fptr, filename, &status);
    if (verbose)
        std::cout << "fits create status: " << status << std::endl;
    fits_report_error(stderr, status);
    
    for (unsigned int i=1; i<sizeof(bitPixDefs)/sizeof(int); i++) {
        if (dataType == dataTypeDefs[i])  {
            bitpix = bitPixDefs[i];
            break;
        }
    }
    fits_create_img(fptr, bitpix, ndims, dims, &status);
    if (verbose)
        std::cout << "fits create img status: " << status << std::endl;
    fits_report_error(stderr, status);
    
    fpixel = (long *) malloc(ndims*sizeof(long));
    
    if (verbose)
        std::cout << "ndims = " << ndims << std::endl;
    if (verbose)
        std::cout << "dims = " << dims << ": " << dims[0] << "," << dims[1] << "," << dims[2] << std::endl;
    long nElements = 1;
    for (int i=0; i<ndims; i++) {
        nElements *= dims[i];
        fpixel[i] = 1;
    }
    if (verbose)
        std::cout << "nElements = " << nElements << std::endl;
    if (verbose)
        std::cout << "data = " << data << std::endl;
    
    assert(fptr);
    assert(fpixel);
    assert(nElements);
    assert(data);
    
    fits_write_pix(fptr, dataType, fpixel, nElements, data, &status);
    fits_report_error(stderr, status);
    
    fits_close_file(fptr, &status);
}

void load_vec(const char *filename, arma::vec& fillVec, int verbose) {
    void *data = NULL;
    int nDims = 0;
    long *dims = NULL;
    int dataType = 0;
    
    
    read_fits_array(filename, &data, nDims, &dims, dataType, verbose);
    if (verbose) {
        std::cout << "data = " << data << std::endl;
        std::cout << "read fits array : nDims = " << nDims << std::endl;
        std::cout << "dims = " << dims << std::endl;
        for (int i=0; i<nDims; i++)
            std::cout << "dim " << i << ": " << dims[i] << std::endl;
    }
    if (!(nDims == 1 || dims[1] == 1)) {
        std::cout << "vec file has other than 1 dimension1!!" << std::endl;
        assert(NULL);
    }
    fillVec.set_size(dims[0]);
    switch (dataType) {
        case TFLOAT:
            for (int i = 0; i < dims[0]; i++)
                fillVec(i) = *((float *)data + i);
            break;
            
        case TDOUBLE:
            for (int i = 0; i < dims[0]; i++)
                fillVec(i) = *((double *)data + i);
            break;
    }
    free(data);
}

void load_mat(const char *filename, arma::mat& fillMat, int verbose) {
    void *data = NULL;
    int nDims = 0;
    long *dims = NULL;
    int dataType = 0;
    
    
    read_fits_array(filename, &data, nDims, &dims, dataType, verbose);
    if (verbose) {
        std::cout << "data = " << data << std::endl;
        std::cout << "read fits array : nDims = " << nDims << std::endl;
        std::cout << "dims = " << dims << std::endl;
        for (int i=0; i<nDims; i++)
            std::cout << "dim " << i << ": " << dims[i] << std::endl;
    }
    if (nDims != 2) {
        std::cout << "mat file has other than 2 dimensions!!" << std::endl;
        assert(NULL);
    }
    fillMat.set_size(dims[1], dims[0]); // fits is tranposed?!?
    switch (dataType) {
        case TFLOAT:
            for (int j = 0; j < dims[1]; j++) {
                for (int i = 0; i < dims[0]; i++) {
                    fillMat(j,i) = *((float *)data + j*dims[0] + i);
                }
            }
            break;
            
        case TDOUBLE:
            for (int j = 0; j < dims[1]; j++) {
                for (int i = 0; i < dims[0]; i++) {
                    fillMat(j,i) = *((double *)data + j*dims[0] + i);
                }
            }
            break;
    }
    free(data);
}

void save_cube(const char *filename, arma::cube& saveCube, int verbose) {
    if (verbose)
        std::cout << "saving the cube " << filename << std::endl;
    
    int nDims = 3;
    long dims[3];
    int dataType = TDOUBLE;
    arma::cube ct = saveCube;
    double *data = (double(*))&ct(0,0,0);
    
    // fits is transposed so...
    for (int s=0; s<ct.n_slices; s++) {
        inplace_strans(ct.slice(s));
    }
    
    dims[0] = saveCube.n_rows; // but use original row/column sizes
    dims[1] = saveCube.n_cols;
    dims[2] = saveCube.n_slices;
    
    write_fits_array(filename, data, nDims, dims, dataType, verbose);
}

void save_cube(const char *filename, arma::icube& saveCube, int verbose) {
    if (verbose)
        std::cout << "saving the cube " << filename << std::endl;
    
    int nDims = 3;
    long dims[3];
    int dataType = TLONG;
    arma::icube ct = saveCube;
    long *data = (long(*))&ct(0,0,0);
    
    // fits is transposed so...
    for (int s=0; s<ct.n_slices; s++) {
        inplace_strans(ct.slice(s));
    }
    
    dims[0] = saveCube.n_rows; // but use original row/column sizes
    dims[1] = saveCube.n_cols;
    dims[2] = saveCube.n_slices;
    
    write_fits_array(filename, data, nDims, dims, dataType, verbose);
}


int load_cube(const char *filename, arma::cube& fillCube, int verbose) {
    void *data = NULL;
    int nDims = 0;
    long *dims = NULL;
    int dataType = 0;
    int status;
    
    
    if (status = read_fits_array(filename, &data, nDims, &dims, dataType, verbose))
        return status;
    
    if (verbose) {
        std::cout << "data = " << data << std::endl;
        std::cout << "read fits array : nDims = " << nDims << std::endl;
        std::cout << "dims = " << dims << std::endl;
    }
    if (!(nDims == 3 || nDims == 2)) {
        std::cout << "nDims = " << nDims << std::endl;
        std::cout << "mat file has other than 2 or 3 dimensions!!" << std::endl;
        assert(NULL);
    }
    if (nDims == 2) {
        long tdims[2];
        tdims[0] = dims[0];
        tdims[1] = dims[1];
        delete[] dims;
        dims = new long[3];
        dims[0] = tdims[0];
        dims[1] = tdims[1];
        dims[2] = 1;
    }
    
    fillCube.set_size(dims[1], dims[0], dims[2]); // fits is tranposed?!?
    switch (dataType) {
            
        case TFLOAT:
            for (int s = 0; s < dims[2]; s++) {
                for (int j = 0; j < dims[1]; j++) {
                    for (int i = 0; i < dims[0]; i++) {
                        fillCube(j,i,s) = *((float *)data + s*dims[0]*dims[1] + j*dims[0] + i);
                    }
                }
            }
            break;
            
        case TDOUBLE:
            for (int s = 0; s < dims[2]; s++) {
                for (int j = 0; j < dims[1]; j++) {
                    for (int i = 0; i < dims[0]; i++) {
                        fillCube(j,i,s) = *((double *)data + s*dims[0]*dims[1] + j*dims[0] + i);
                    }
                }
            }
            break;
    }
    free(data);
    return 0;
}


int load_cube(const char *filename, arma::icube& fillCube, int verbose) {
    void *data = NULL;
    int nDims = 0;
    long *dims = NULL;
    int dataType = 0;
    int status;
    
    if (status = read_fits_array(filename, &data, nDims, &dims, dataType, verbose))
        return status;
    
    if (verbose) {
        std::cout << "data = " << data << std::endl;
        std::cout << "read fits array : nDims = " << nDims << std::endl;
        std::cout << "dims = " << dims << std::endl;
    }
    if (!(nDims == 3 || nDims == 2)) {
        std::cout << "nDims = " << nDims << std::endl;
        std::cout << "mat file has other than 2 or 3 dimensions!!" << std::endl;
        assert(NULL);
    }
    if (nDims == 2) {
        long tdims[2];
        tdims[0] = dims[0];
        tdims[1] = dims[1];
        delete[] dims;
        dims = new long[3];
        dims[0] = tdims[0];
        dims[1] = tdims[1];
        dims[2] = 1;
    }
    
    fillCube.set_size(dims[1], dims[0], dims[2]); // fits is tranposed?!?
    switch (dataType) {
        case TLONG:
            for (int s = 0; s < dims[2]; s++) {
                for (int j = 0; j < dims[1]; j++) {
                    for (int i = 0; i < dims[0]; i++) {
                        fillCube(j,i,s) = *((long *)data + s*dims[0]*dims[1] + j*dims[0] + i);
                    }
                }
            }
            break;
    }
    free(data);
    return 0;
}

void load_mat(const char *filename, arma::cx_mat& fillMat, int verbose) {
    arma::mat mRe;
    arma::mat mIm;
    
    std::string fString = filename;
    char nameRoot[200];
    char loadName[200];
    
    int fitPos = fString.find(".fits");
    if (fitPos > 0) {
        fString.copy(nameRoot, fitPos);
        nameRoot[fitPos] = '\0';
    } else
        strcpy(nameRoot, filename);
    if (verbose)
        std::cout << "nameRoot: " << nameRoot << std::endl;
    
    strcpy(loadName, nameRoot);
    strcat(loadName, "_re.fits");
    load_mat(loadName, mRe, verbose);
    strcpy(loadName, nameRoot);
    strcat(loadName, "_im.fits");
    load_mat(loadName, mIm, verbose);
    
    fillMat.set_size(mRe.n_rows, mRe.n_cols);
    fillMat.set_real(mRe);
    fillMat.set_imag(mIm);
}

void save_mat(const char *filename, arma::mat& saveMat, int verbose) {
    if (verbose)
        std::cout << "saving the matrix " << filename << std::endl;
    arma::mat tMat = saveMat.t(); // fits is transposed
    double *data = (double(*))&tMat(0,0);
    int nDims = 2;
    long dims[2];
    int dataType = TDOUBLE;
    
    dims[0] = saveMat.n_cols; // but use original row/column sizes
    dims[1] = saveMat.n_rows;
    
    write_fits_array(filename, data, nDims, dims, dataType, verbose);
}

void save_mat(const char *filename, arma::cx_mat& saveMat, const char *saveType, int verbose) {
    std::string fString = filename;
    char nameRoot[200];
    char saveName[200];
    
    int fitPos = fString.find(".fits");
    if (fitPos > 0) {
        fString.copy(nameRoot, fitPos);
        nameRoot[fitPos] = '\0';
    } else
        strcpy(nameRoot, filename);
    if (verbose)
        std::cout << "nameRoot: " << nameRoot << std::endl;
    
    if (!strcmp(saveType, "reIm")) {
        strcpy(saveName, nameRoot);
        strcat(saveName, "_re.fits");
        arma::mat tmp = real(saveMat);
        save_mat(saveName, tmp, verbose);
        strcpy(saveName, nameRoot);
        strcat(saveName, "_im.fits");
        tmp = imag(saveMat);
        save_mat(saveName, tmp, verbose);
    }
    if (!strcmp(saveType, "amPh")) {
        strcpy(saveName, nameRoot);
        strcat(saveName, "_am.fits");
        arma::mat tmp = abs(saveMat);
        save_mat(saveName, tmp, verbose);
        strcpy(saveName, nameRoot);
        strcat(saveName, "_ph.fits");
        tmp = arg(saveMat);
        save_mat(saveName, tmp, verbose);
    }
    
}

    
    
    
    
