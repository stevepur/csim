bin_PROGRAMS = csim csim_gpu test
csim_SOURCES = csim.cpp
csim_gpu_SOURCES = csim.cpp
test_SOURCES = test.cpp

csim_LDADD = lib/libcsim.la \
coronagraph/libcoronagraph.la \
data/libdata.la \
elements/libelements.la \
optimization/liboptimization.la \
telescope/libtelescope.la \
tools/libtools.la \
/Users/steve/programming_projects/gnu/fftw-3.3.7/libfftw3.la \
/Users/steve/programming_projects/gnu/fftw-3.3.7/threads/libfftw3_omp.la \
/Users/steve/programming_projects/gnu/cfitsio/libcfitsio.a \
-lblas -llapack

csim_CPPFLAGS = \
-I/Users/steve/programming_projects/gnu/cfitsio/include \
-I/Users/steve/programming_projects/gnu/armadillo-7.800.1/include \
-I/Users/steve/programming_projects/gnu/fftw-3.3.7/include \
-Ilib \
-Icoronagraph \
-Idata \
-Ielements \
-Itelescope \
-Ioptimization \
-Ioptimization/differentialEvolution \
-Itools

csim_gpu_LDADD = lib/libcsim.la \
coronagraph/libcoronagraph.la \
data/libdata.la \
elements/libelements.la \
optimization/liboptimization.la \
telescope/libtelescope.la \
tools/libtools.la \
/Users/steve/programming_projects/gnu/fftw-3.3.7/libfftw3.la \
/Users/steve/programming_projects/gnu/fftw-3.3.7/threads/libfftw3_omp.la \
/Users/steve/programming_projects/gnu/cfitsio/libcfitsio.a \
-L/Developer/NVIDIA/CUDA-8.0/lib -lnvblas -lblas -llapack

csim_gpu_CPPFLAGS = \
-I/Users/steve/programming_projects/gnu/cfitsio/include \
-I/Users/steve/programming_projects/gnu/armadillo-7.800.1/include \
-I/Users/steve/programming_projects/gnu/fftw-3.3.7/include \
-Ilib \
-Icoronagraph \
-Idata \
-Ielements \
-Itelescope \
-Ioptimization \
-Ioptimization/differentialEvolution \
-Itools

test_LDADD = lib/libcsim.la \
coronagraph/libcoronagraph.la \
data/libdata.la \
elements/libelements.la \
optimization/liboptimization.la \
telescope/libtelescope.la \
tools/libtools.la \
/Users/steve/programming_projects/gnu/fftw-3.3.7/libfftw3.la \
/Users/steve/programming_projects/gnu/fftw-3.3.7/threads/libfftw3_omp.la \
/Users/steve/programming_projects/gnu/cfitsio/libcfitsio.a \
-lblas -llapack

test_CPPFLAGS = \
-I/Users/steve/programming_projects/gnu/cfitsio/include \
-I/Users/steve/programming_projects/gnu/armadillo-7.800.1/include \
-I/Users/steve/programming_projects/gnu/fftw-3.3.7/include \
-Ilib \
-Icoronagraph \
-Idata \
-Ielements \
-Ioptimization \
-Ioptimization/differentialEvolution \
-Itelescope \
-Itools

SUBDIRS = lib coronagraph data elements optimization telescope tools
