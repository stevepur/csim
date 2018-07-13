bin_PROGRAMS = csim csim_gpu test
csim_SOURCES = csim.cpp
csim_gpu_SOURCES = csim.cpp
test_SOURCES = test.cpp

csim_LDADD = tools/libtools.la \
optimization/liboptimization.la \
coronagraph/libcoronagraph.la \
elements/libelements.la \
telescope/libtelescope.la \
data/libdata.la \
lib/libcsim.la \
$(MY_LIBDIR)/fftw-3.3.7/libfftw3.la \
$(MY_LIBDIR)/fftw-3.3.7/threads/libfftw3_omp.la \
$(MY_LIBDIR)/cfitsio/libcfitsio.a \
-lblas -llapack

csim_CPPFLAGS = \
-I$(MY_LIBDIR)/cfitsio/include \
-I$(MY_LIBDIR)/armadillo-8.500.1/include \
-I$(MY_LIBDIR)/fftw-3.3.7/include \
-Ilib \
-Icoronagraph \
-Idata \
-Ielements \
-Itelescope \
-Ioptimization \
-Ioptimization/differentialEvolution \
-Itools

csim_gpu_LDADD = tools/libtools.la \
optimization/liboptimization.la \
coronagraph/libcoronagraph.la \
elements/libelements.la \
telescope/libtelescope.la \
data/libdata.la \
lib/libcsim.la \
$(MY_LIBDIR)/fftw-3.3.7/libfftw3.la \
$(MY_LIBDIR)/fftw-3.3.7/threads/libfftw3_omp.la \
$(MY_LIBDIR)/cfitsio/libcfitsio.a \
-L/Developer/NVIDIA/CUDA-9.2/lib -lnvblas -lblas -llapack

csim_gpu_CPPFLAGS = \
-I$(MY_LIBDIR)/cfitsio/include \
-I$(MY_LIBDIR)/armadillo-8.500.1/include \
-I$(MY_LIBDIR)/fftw-3.3.7/include \
-Ilib \
-Icoronagraph \
-Idata \
-Ielements \
-Itelescope \
-Ioptimization \
-Ioptimization/differentialEvolution \
-Itools

test_LDADD = tools/libtools.la \
optimization/liboptimization.la \
coronagraph/libcoronagraph.la \
elements/libelements.la \
telescope/libtelescope.la \
data/libdata.la \
lib/libcsim.la \
$(MY_LIBDIR)/fftw-3.3.7/libfftw3.la \
$(MY_LIBDIR)/fftw-3.3.7/threads/libfftw3_omp.la \
$(MY_LIBDIR)/cfitsio/libcfitsio.a \
-lblas -llapack

test_CPPFLAGS = \
-I$(MY_LIBDIR)/cfitsio/include \
-I$(MY_LIBDIR)/armadillo-8.500.1/include \
-I$(MY_LIBDIR)/fftw-3.3.7/include \
-Ilib \
-Icoronagraph \
-Idata \
-Ielements \
-Ioptimization \
-Ioptimization/differentialEvolution \
-Itelescope \
-Itools

SUBDIRS = lib coronagraph data elements optimization telescope tools
