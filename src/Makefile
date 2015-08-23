# LAMMPS multiple-machine Makefile

SHELL = /bin/bash

# Definitions
ROOT =	pasc
EXE =	$(ROOT)_$@
SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)
TESTDIR = $(PASCAL_SRC_DIR)/../applications/test
# Package variables
# List of all targets

help:
	@echo ''
	@echo '***PASCAL MAKE SYSTEM***'
	@echo 'make clean-all           delete all object files'
	@echo 'make makelib             create Makefile.lib for static library build'
	@echo 'make makeshlib           create Makefile.shlib for shared library build'
	@echo 'make -f Makefile.lib machine      build LAMMPS as static library for machine'
	@echo 'make -f Makefile.shlib machine    build LAMMPS as shared library for machine'
	@echo 'make machine             build LAMMPS where machine is one of:'
	@echo ''
	@files="`ls MAKE/Makefile.*`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''

# Build the code
.DEFAULT:
	@test -f MAKE/Makefile.$@
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@$(SHELL) Make.sh style
#	@$(SHELL) Make.sh models
	@cp MAKE/Makefile.$@ Obj_$@/Makefile
	@cd Obj_$@; \
	$(MAKE) $(MFLAGS) "OBJ = $(OBJ)" "INC = $(INC)" "SHFLAGS =" \
	  "EXE = ../$(EXE)" ../$(EXE)

# build thirdParty
thirdParty:
	@cd $(PASCAL_SUNDIALS_DIR); \
	$(SHELL) ./configure --prefix=\${PASCAL_INST_DIR} --with-cflags=-fPIC
	@cd $(PASCAL_SUNDIALS_DIR); \
	make
	@cd $(PASCAL_SUNDIALS_DIR); \
	make install
	@cd $(PASCAL_THIRDPARTY_DIR)/chemkinReader/src/; \
	make -f Makefile.shlib fedora_fpic; \
	make -f Makefile.lib fedora_fpic

# build test programs
test:
	@cd $(PASCAL_SUNDIALS_DIR)/../../applications/test/CVODE_AdDiff_1D; \
	make

# Remove all object files, static libraries, and libs of thirdparty
clean-all:
	rm -rf Obj_*  lib$(ROOT)*.a
	@cd $(PASCAL_SUNDIALS_DIR); \
	make clean
	rm -rf $(PASCAL_SUNDIALS_DIR)/libtool
	rm -rf $(PASCAL_SUNDIALS_DIR)/src/*/Makefile
	rm -rf $(PASCAL_INST_DIR)/include
	rm -rf $(PASCAL_INST_DIR)/lib64
	rm -rf $(PASCAL_INST_DIR)/*/*
	@cd $(PASCAL_THIRDPARTY_DIR)/chemkinReader/src/; \
	make clean-all

# NOTIMPLEMENTED
makelib:
	@$(SHELL) Make.sh style
	@$(SHELL) Make.sh Makefile.lib

# NOTIMPLEMENTED
makeshlib:
	@$(SHELL) Make.sh style
	@$(SHELL) Make.sh Makefile.shlib

# NOTIMPLEMENTED - install LAMMPS shared lib and Python wrapper for Python usage
install-python:
	@python ../python/install.py
