#######################################################################
#
#   Makefile for DPMJET-III
#
#######################################################################
#choose gnu, intel or g95
COMP_TYPE=gfortran
COMP_MODE=release

# Extension for PYTHON library
LEXT?=.so
LIB_DIR?=lib
# Version
# BUILD_VERSION = 19.1.0
# BUILD_DATE = 2019/09/02
# CPPFLAGS =
#-DFOR_FLUKA 
#######################################################################
#
#   handling of different compilers
#
#######################################################################
ifeq ($(COMP_TYPE),g95)
	FC = g95
	CC = c++
	F2PY_COMP='--fcompiler=g95'
	ifeq ($(COMP_MODE),debug)
		FOPT = -Wall -fbounds-check -Wno-unused -O0 -g \
			   -ftrace=full -ftrapping-math -pedantic -standard=legacy
	else
		FOPT = -O3 -g -fmultiple-save -fexpensive-optimizations -funroll-loops
	endif

	F90FLAGS = -cpp -ffree-form -Wobsolescent -fno-second-underscore

else ifeq ($(COMP_TYPE),gfortran)
		FC = gfortran
		CC = g++
		# Make sure the flags are selected according to OS
		F2PY_COMP = --fcompiler=gnu95 

		ifeq ($(COMP_MODE),debug)
			FOPT = -fPIC -Wall  -Wno-maybe-uninitialized \
				   -fbounds-check -O0 -g -ffpe-trap=invalid,zero,overflow \
				   -Wno-unused-dummy-argument
	   	else
	   		FOPT = -fPIC -mtune=native -Og -Wno-uninitialized -Wno-unused-dummy-argument -std=legacy
	   	endif

	   	F90FLAGS = -xf77-cpp-input -ffree-form -fno-second-underscore

else
	FC = ifort
	CC = icc
	F2PY_COMP='--fcompiler=intelem'
	ifeq ($(COMP_MODE),debug)
		FOPT = -stand f90 -traceback -gen-interfaces -warn all -O0 -g \
			   -fpe:0 -diag-disable 6717 -module ./src/mod -I ./src/mod
	else
		FOPT = -fast -g -module ./src/mod -I ./src/mod -no-ipo
	endif

	F90FLAGS = -fpp -free
endif

#Linker
LD = $(FC)

# Files
PHOJET_INCS = ./include/phojet
PHOJET_SRCS :=$(wildcard ./src/phojet/*.f)
PHOJET_OBJS :=$(PHOJET_SRCS:.f=.o)
DPMJET_INCS = ./include/dpmjet
DPMJET_FLUKA_INCS = ./include/flinclude #./flukaadd
DPMJET_SRCS :=$(wildcard ./src/dpmjet/*.f)
DPMJET_OBJS :=$(DPMJET_SRCS:.f=.o)
PYTHIA_INCS = ./include/pythia
PYTHIA_SRCS :=$(wildcard ./src/pythia/*.f)
PYTHIA_OBJS :=$(PYTHIA_SRCS:.f=.o)

APP_SRCS :=$(wildcard ./src/exe/*.f)
APP_OBJS :=$(APP_SRCS:.f=.o)
APP_EXE :=$(APP_OBJS:.o=)
APP_EXE :=$(subst ./src/exe/,,$(APP_EXE))

DUMMY_SRCS :=$(wildcard ./common/*.f)
DUMMY_OBJS :=$(DUMMY_SRCS:.f=.o)

DPMJET_FUNCS = pho_event dt_init dt_kkinc \
idt_icihad dt_xsglau pycomp dt_initjs dt_rndmst dt_rndm dt_inucas idt_ipdgha dt_evtout
DPMJET_FUNCS += pho_init pho_setpar poevt1 poevt2 pho_pname pho_pmass pho_setmdl \
pho_setpdf pycomp pho_xsect pho_borncs pho_harmci pho_fitout pho_mcini pho_ptcut \
pytune pho_rregpar pho_sregpar pho_prevnt ipho_pdg2id ipho_id2pdg pho_harint \
impy_openlogfile impy_closelogfile pho_harxto pho_harxpt pho_setpcomb \
dt_phoxs dt_xshn dt_flahad dt_title pho_ghhias

INCLU = -I$(PYTHIA_INCS) -I$(PHOJET_INCS) -I$(DPMJET_INCS) -I$(DPMJET_FLUKA_INCS)

pylib = dpmjetIII191$(LEXT)
F2PY = python -m numpy.f2py

all: exe 

.PHONY: pylib
pylib: $(pylib)

$(pylib): lib/libDPMJET.a common/dpmjetIII191.pyf
	$(F2PY) -c $(F2PY_COMP) --opt="$(FOPT)" \
	     $(INCLU) common/dpmjetIII191.pyf -Llib -lDPMJET

.PHONY: install
install: $(pylib)
	cp *$(LEXT) $(LIB_DIR)

.PHONY: exe
exe: $(APP_OBJS) lib/libDPMJET.a
	for exec in $(APP_EXE) ; do \
        $(LD) -o bin/$$exec ./src/exe/$$exec.o -Llib -lDPMJET  ; \
	done

common/dpmjetIII191.pyf:
	cat $(PHOJET_SRCS) $(PYTHIA_SRCS) $(DPMJET_SRCS) $(DUMMY_SRCS) > f2pytemp.f
	gfortran -E -cpp f2pytemp.f > f2py_cpp.f
	$(F2PY) -m dpmjetIII191 -h common/dpmjetIII191.pyf \
	--include-paths $(DPMJET_INCS):$(PHOJET_INCS):$(PYTHIA_INCS):$(DPMJET_FLUKA_INCS) \
	--overwrite-signature only: $(DPMJET_FUNCS) : f2py_cpp.f
	rm -f f2pytemp.f f2py_cpp.f f2pytemp.s

lib/libDPMJET.a:  $(PHOJET_OBJS) $(PYTHIA_OBJS) $(DPMJET_OBJS) $(DUMMY_OBJS)
	ar -cr lib/libDPMJET.a $(DPMJET_OBJS) $(PHOJET_OBJS) $(PYTHIA_OBJS) $(DUMMY_OBJS)

.f.o:
	$(FC) -c -cpp $(CPPFLAGS) $(FOPT) $(INCLU) -o $@ $<   

.PHONY: clean
clean:
	rm -rf lib/* *.so common/*.o *.dSYM 
	rm -f *.o src/pythia/*.o src/phojet/*.o src/dpmjet/*.o src/exe/*.o 
	rm -f *.s src/pythia/*.s src/phojet/*.s src/dpmjet/*.s src/exe/*.s
	rm -f bin/*

.PHONY: distclean
distclean: clean
	rm -rf common/dpmjetIII191.pyf
