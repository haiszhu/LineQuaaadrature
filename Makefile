# LineQuaaadrature Makefile
# Builds: libLineQuaaadrature.a + LineQuaaadrature_mex.<ext>
#
# Goals:
#   1. linequaaadrature_mod + solidangle_mod verification via MEX
#   2. quad precision support (real(16) / __float128 via -lquadmath)
#
# Usage:
#   make -f Makefile          — show build targets
#   make -f Makefile lib      — build static library
#   make -f Makefile mex      — build MATLAB MEX gateway
#   make -f Makefile test     — build solid-angle r64 test
#   make -f Makefile test_r128 — build solid-angle r128 test
#   make -f Makefile clean    — remove build artifacts

SHELL := /bin/bash

ROOT      := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
ROOT      := $(patsubst %/,%,$(ROOT))
SRC_DIR   := $(ROOT)/src
MATLAB_DIR := $(ROOT)/matlab
BLD_DIR   := $(ROOT)/build

# ---- compilers ----
FC := gfortran-15
CC := gcc-15
MW := ~/mwrap/mwrap
MWFLAGS := -c99complex -i8 -mex

# ---- platform detection ----
UNAME := $(shell uname)
ARCH  := $(shell uname -m)

ifeq ($(UNAME), Darwin)
  MATLAB_ROOT  := $(shell ls -d /Applications/MATLAB_R*.app 2>/dev/null | sort | tail -n1)
  ifeq ($(ARCH), arm64)
    MATLAB_ARCH  := maca64
    MEX_EXT      := mexmaca64
    OPENBLAS_DIR := /opt/homebrew/opt/openblas-singlethread
  else
    MATLAB_ARCH  := maci64
    MEX_EXT      := mexmaci64
    OPENBLAS_DIR := /usr/local/opt/openblas-singlethread
  endif
  MATLAB_INC   := -I$(MATLAB_ROOT)/extern/include
  MATLAB_LIBS  := $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmx.dylib \
                  $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmex.dylib \
                  $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmat.dylib -lm
  MEX_LDFLAGS  := -bundle -Wl,-undefined,dynamic_lookup
else
  MATLAB_ROOT  := $(shell ls -d /usr/local/MATLAB/R* 2>/dev/null | sort | tail -n1)
  MATLAB_ARCH  := glnxa64
  MEX_EXT      := mexa64
  MATLAB_INC   := -I$(MATLAB_ROOT)/extern/include
  MATLAB_LIBS  := -L$(MATLAB_ROOT)/bin/$(MATLAB_ARCH) -lmx -lmex -lmat -lm
  MEX_LDFLAGS  := -shared
endif

OPENBLAS_LIBS := -L$(OPENBLAS_DIR)/lib -lopenblas
HDF5_ROOT ?= /opt/homebrew/opt/hdf5
HDF5_INC := -I$(HDF5_ROOT)/include
HDF5_LIBS := -L$(HDF5_ROOT)/lib -lhdf5

# ---- Fortran flags ----
# -fdefault-integer-8 : 8-byte integers (matches mwrap -i8)
# -freal-4-real-16    : TODO enable for quad precision build
# -lquadmath          : TODO link quad math library for real(16) support
FFLAGS := -g -O3 -fPIC \
           -ffp-contract=off -fno-unsafe-math-optimizations \
           -fallow-argument-mismatch -std=legacy -w \
           -fdefault-integer-8 \
           -frecursive \
           -cpp \
           -march=native -funroll-loops \
           -fopenmp \
           -J$(BLD_DIR) -I$(BLD_DIR) -I$(SRC_DIR)

# ---- sources and objects ----
# TODO: add additional .f90 files to LQ_SOURCES as needed
LQ_SOURCES := $(SRC_DIR)/linequaaadrature_mod.f90 \
              $(SRC_DIR)/linequaaadrature_mex.f90 \
              $(SRC_DIR)/koorn_geom_mod.f90 \
              $(SRC_DIR)/kernel_mod.f90 \
			  $(SRC_DIR)/kernel_mex.f90 \
              $(SRC_DIR)/adaptive_mod.f90 \
			  $(SRC_DIR)/adaptive_mex.f90 \
	              $(SRC_DIR)/solidangle_mod.f90 \
				  $(SRC_DIR)/solidangle_mex.f90 \
	              $(SRC_DIR)/ellipsoid_mesh_mod.f90 \
	              $(SRC_DIR)/lap3d_mod.f90
	              

LQ_OBJECTS := $(patsubst $(SRC_DIR)/%.f90, $(BLD_DIR)/%.o, $(LQ_SOURCES)) \
		              $(BLD_DIR)/hdf5_io.o

LIB := $(BLD_DIR)/libLineQuaaadrature.a

# ---- mwrap-generated gateway ----
MW_SRC  := $(MATLAB_DIR)/LineQuaaadrature.mw
MEX_C   := $(MATLAB_DIR)/LineQuaaadrature_mex.c
MEX_OUT := $(MATLAB_DIR)/LineQuaaadrature_mex.$(MEX_EXT)

# ============================================================

TEST_SRC     := $(ROOT)/test/solid_angle/test_solid_angle.f90
TEST_BIN     := $(BLD_DIR)/test_solid_angle
TEST_R128_SRC := $(ROOT)/test/solid_angle/test_solid_angle_r128.f90
TEST_R128_BIN := $(BLD_DIR)/test_solid_angle_r128

.PHONY: all mex lib test test_r128 clean

all:
	@echo "LineQuaaadrature build targets"
	@echo ""
	@echo "  make -f Makefile mex        build LineQuaaadrature_mex.$(MEX_EXT)"
	@echo "  make -f Makefile lib        build build/libLineQuaaadrature.a"
	@echo "  make -f Makefile test       build build/test_solid_angle"
	@echo "  make -f Makefile test_r128  build build/test_solid_angle_r128"
	@echo "  make -f Makefile clean      remove build artifacts"
	@echo ""

mex: $(MEX_OUT)

lib: $(LIB)

test: $(TEST_BIN)

test_r128: $(TEST_R128_BIN)

$(TEST_BIN): $(LIB) $(TEST_SRC) | $(BLD_DIR)
	$(FC) $(FFLAGS) -fopenmp $(TEST_SRC) -L$(BLD_DIR) -lLineQuaaadrature \
	  $(OPENBLAS_LIBS) -lgfortran -lm -o $(TEST_BIN)

$(TEST_R128_BIN): $(LIB) $(TEST_R128_SRC) | $(BLD_DIR)
	$(FC) $(FFLAGS) $(TEST_R128_SRC) -L$(BLD_DIR) -lLineQuaaadrature \
	  $(OPENBLAS_LIBS) -lgfortran -lm -o $(TEST_R128_BIN)

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

# ---- compile Fortran objects ----
$(BLD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BLD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BLD_DIR)/hdf5_io.o: $(SRC_DIR)/hdf5_io.c | $(BLD_DIR)
	$(CC) -c -fPIC $(HDF5_INC) $< -o $@

$(BLD_DIR)/kernel_mod.o: $(BLD_DIR)/linequaaadrature_mod.o

$(BLD_DIR)/adaptive_mod.o: $(BLD_DIR)/linequaaadrature_mod.o

$(BLD_DIR)/solidangle_mod.o: $(BLD_DIR)/linequaaadrature_mod.o \
                             $(BLD_DIR)/koorn_geom_mod.o \
                             $(BLD_DIR)/kernel_mod.o \
                             $(BLD_DIR)/adaptive_mod.o

$(BLD_DIR)/linequaaadrature_mex.o: $(BLD_DIR)/linequaaadrature_mod.o \
                                   $(BLD_DIR)/adaptive_mod.o \
                                   $(BLD_DIR)/solidangle_mod.o \
                                   $(BLD_DIR)/ellipsoid_mesh_mod.o

$(BLD_DIR)/solidangle_mex.o: $(BLD_DIR)/solidangle_mod.o

$(BLD_DIR)/kernel_mex.o: $(BLD_DIR)/kernel_mod.o $(BLD_DIR)/solidangle_mod.o

$(BLD_DIR)/adaptive_mex.o: $(BLD_DIR)/adaptive_mod.o

# note: module files (.mod) are emitted to BLD_DIR via -J flag
# TODO: add -J$(BLD_DIR) once modules are non-empty

# ---- static library ----
$(LIB): $(LQ_OBJECTS)
	rm -f $@
	ar rcs $@ $^

# ---- mwrap: two-pass generation ----
$(MEX_C): $(MW_SRC) | $(BLD_DIR)
	cd $(MATLAB_DIR) && $(MW) $(MWFLAGS) LineQuaaadrature_mex -mb -list LineQuaaadrature.mw
	cd $(MATLAB_DIR) && $(MW) $(MWFLAGS) LineQuaaadrature_mex -c LineQuaaadrature_mex.c LineQuaaadrature.mw
	perl -pi -e 's/_{2,}/_/g' $(MEX_C)

# ---- link MEX ----
$(MEX_OUT): $(LIB) $(MEX_C)
	$(CC) $(MEX_LDFLAGS) -fPIC \
	  -DMATLAB_MEX_FILE -DMATLAB_DEFAULT_RELEASE=R2018a -DMX_COMPAT_32=0 \
	  $(MATLAB_INC) \
	  $(MEX_C) \
	  -L$(BLD_DIR) -lLineQuaaadrature \
	  $(HDF5_LIBS) \
	  $(MATLAB_LIBS) \
	  $(OPENBLAS_LIBS) \
		  -lgfortran -lquadmath -lm \
		  -o $(MEX_OUT)

# ---- clean ----
clean:
	rm -rf $(BLD_DIR) $(MEX_OUT)
