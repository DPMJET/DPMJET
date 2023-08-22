#######################################################################
#
#   Makefile for DPMJET-III
#
#######################################################################
#choose gnu, intel or g95
CVendor = "GNU"
Config?=Release

#######################################################################
#
#   compiler
#
#######################################################################

ifeq ($(CVendor),"GNU")
	#  GNU
	FC := $(or $(FC), gfortran)
else
	#  Intel
	FC = ifort
endif

#######################################################################
#
#   compiler options for different platforms
#
#######################################################################
ifeq ($(CVendor),"GNU")
	ifeq ($(Config),Debug)
		# GNU Debug
		OPT = -fPIC -Wall -fbounds-check -O0 -g \
			  -ffpe-trap=invalid,zero,overflow -Wuninitialized -std=legacy
		OPTF90 = -fPIC -Wall -fbounds-check -O0 -g \
			  -ffpe-trap=invalid,zero,overflow -Wuninitialized \
			  -fno-second-underscore
		#OPT = -fPIC -Wall -Wno-uninitialized -Wno-unused-variable -O3 -g -ffpe-trap=invalid,zero,overflow
	else
		# GNU Release
		OPT = -O3 -Wno-uninitialized -fPIC -std=legacy
		OPTF90 = -O3 -Wno-uninitialized -fPIC -fno-second-underscore
	endif
else
	ifeq ($(Config),Debug)
	# Intel Debug (-gen-interfaces -warn interfaces)
		OPT = -check bounds -O0 -g -check pointer -fpe0 -traceback
		OPTF90 = -check bounds -O0 -g -check pointer -fpe0 -traceback \ 
			 -cpp -ffree-form -Wobsolescent -fno-second-underscore
	else
		# Intel Release
		OPTF90 = -fast -fpe0 \ 
		      -cpp -ffree-form -Wobsolescent -fno-second-underscore
		OPT = -fast -fpe0	
	endif
endif

#Linker
LD = $(FC)

# Directories
WORK_DIR = $(CURDIR)
LIB_DIR?="$(WORK_DIR)/lib"
# New line define
define \n


endef
#######################################################################
#
#   Files
#
#######################################################################
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

# Portability (I know that this is insane...)
ifeq ($(OS),Windows_NT)
  DEL_COMMAND = del /q /f
  MKDIR_COMMAND = if not exist "$(LIB_DIR)" mkdir
  COPY_COMMAND = copy /b
  COPY_DUMP = > nul  2>&1
  CAT_COMMAND = type
  EXESUFX = .exe
  PATHSEP2=\\
  PATHSEP=$(strip $(PATHSEP2))
  space := $(null) #
  comma := ,
  
  PHOJET_SRCS_CMMA := $(subst $(space),$(comma),$(strip $(PHOJET_SRCS)))
  DPMJET_SRCS_CMMA := $(subst $(space),$(comma),$(strip $(DPMJET_SRCS)))
  PYTHIA_SRCS_CMMA := $(subst $(space),$(comma),$(strip $(PYTHIA_SRCS)))
  DUMMY_SRCS_CMMA := $(subst $(space),$(comma),$(strip $(DUMMY_SRCS)))
  
  PHOJET_SRCS_CMMA := $(subst /,\,$(strip $(PHOJET_SRCS_CMMA)))
  DPMJET_SRCS_CMMA := $(subst /,\,$(strip $(DPMJET_SRCS_CMMA)))
  PYTHIA_SRCS_CMMA := $(subst /,\,$(strip $(PYTHIA_SRCS_CMMA)))
  DUMMY_SRCS_CMMA := $(subst /,\,$(strip $(DUMMY_SRCS_CMMA)))
  
else
  DEL_COMMAND = rm -rf
  MKDIR_COMMAND = mkdir -p
  COPY_COMMAND = cp
  COPY_DUMP =
  CAT_COMMAND = cat
  EXESUFX = 
  PATHSEP=/
endif

INCLU = -I$(PYTHIA_INCS) -I$(PHOJET_INCS) -I$(DPMJET_INCS) -I$(DPMJET_FLUKA_INCS)

ifneq ($(FLDOTINCL_DIR),)
  CPPFLAGS = -DFOR_FLUKA -DFLDOTINCL
  INCLU += -I$(FLDOTINCL_DIR)
endif

ifneq ($(FLINCINCL_DIR),)
  CPPFLAGS = -DFOR_FLUKA -DFLINCINCL -DPNUTINC
  INCLU += -I$(FLINCINCL_DIR)
endif

all: exe 

.PHONY: exe
exe: $(APP_OBJS) lib/libDPMJET.a
	$(foreach a, $(APP_EXE), $(LD) -o bin/$(a) ./src/exe/$(a).o -Llib -lDPMJET ${\n})

lib/libDPMJET.a:  $(PHOJET_OBJS) $(PYTHIA_OBJS) $(DPMJET_OBJS) $(DUMMY_OBJS)
	if [ ! -d lib ]; then $(MKDIR_COMMAND) lib; fi
	ar -crs lib/libDPMJET.a $(DPMJET_OBJS) $(PHOJET_OBJS) $(PYTHIA_OBJS) $(DUMMY_OBJS)

.f.o:
	$(FC) -c -cpp $(CPPFLAGS) $(OPT) $(INCLU) -o $@ $<

.PHONY: clean
clean:
	$(DEL_COMMAND) lib$(PATHSEP)libDPMJET.a *.so common$(PATHSEP)*.o *.dSYM $(COPY_DUMP) *.o
	$(DEL_COMMAND) *.o src$(PATHSEP)pythia$(PATHSEP)*.o src$(PATHSEP)phojet$(PATHSEP)*.o src$(PATHSEP)dpmjet$(PATHSEP)*.o src$(PATHSEP)exe$(PATHSEP)*.o $(COPY_DUMP)
	$(DEL_COMMAND) *.s src$(PATHSEP)pythia$(PATHSEP)*.s src$(PATHSEP)phojet$(PATHSEP)*.s src$(PATHSEP)dpmjet$(PATHSEP)*.s src$(PATHSEP)exe$(PATHSEP)*.s $(COPY_DUMP)
	$(DEL_COMMAND) $(addprefix bin$(PATHSEP),$(addsuffix $(EXESUFX), $(APP_EXE))) $(COPY_DUMP)

.PHONY: distclean
distclean: clean

# **************************
# Variable printing target
# make print-INCPATH
# **************************
print-%:
	@echo $* = $($*)
