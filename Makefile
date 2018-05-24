# makefile for mcmax (with comments!)
# Tested on Fedora Core 8 with ifort 10.1.015 (20/12/2012)
# Tested on MacOSX 10.6 with ifort 11.1.080 (20/12/2012)
# Tested on MacOSX 10.8 with ifort 14.0.1 (30/12/2013)
# Tested on MacOSX 10.9 with ifort 12.0.0 (20/12/2013)

# Added gfortran compilation (usage "make compiler=gfortran")
# Tested on Scientific Linux 7.4 with gcc/gfortran 4.8.5 (13/10/2017)
# Tested on MacOSX 10.11.6 with gfortran 4.9.4 (13/10/2017)

	
#GITVERSION = $(echo "#define gitversion = \"$(shell git rev-parse HEAD)\"" > gitversion.h)
GITVERSION = $(ls -l)

# compiler= FC, flags = FFlags
# linker= LINKER, flags= LDFLAGS, libraries=LIBS
#################################################################
# Default compiler path (ifort)
ifneq ($(compiler),gfortran)
	FC	      = ifort
	LINKER	      = ifort
	
	# enforce single core compilation with:
  # cl> make multi=false
  ifneq ($(multi),false)
  	MULTICORE = -openmp -fp-model strict
  	ifeq ($(debug),true)
    		MULTICORE = -openmp
  	endif
  endif
  
  # array boundary check
  ifeq ($(debug),true)
    DEBUGGING = -debug -check bounds -ftrapuv -fpe3 -O0 -check all,noarg_temp_created -fp-stack-check
    #DEBUGGING = -debug -ftrapuv -g -check all -fp-stack-check
    #DEBUGGING = -heap-arrays
    #DEBUGGING = -gen-interfaces -warn interfaces
  endif
  
  # Platform specific compilation options
  FLAG_ALL      = -O3 -extend-source -g -traceback -zero -prec-div $(MULTICORE) $(DEBUGGING)
  FLAG_LINUX    = -msse3 #-prefetch
  FLAG_MAC      = -xHOST -static-intel #-opt-prefetch 
  
  
  ifeq ($(shell uname),Linux)
    FFLAGS   = $(FLAG_ALL) $(FLAG_LINUX) -diag-disable vec
    LDFLAGS  = $(FLAG_ALL) $(FLAG_LINUX) -fpp Version.f
    LIBS     = -lm -lfftw3 -lcfitsio -I/sw/include -L/home/sw-astro/cfitsio/lib -L$(HOME)/lib
  else
    FFLAGS  = $(FLAG_ALL) $(FLAG_MAC) -diag-disable 8290,8291
    LDFLAGS = $(FLAG_ALL) $(FLAG_MAC) -fpp Version.f 
    #LIBS	  =  -L/sw/lib -lm -lfftw3 -lcfitsio -I/sw/include
    LIBS	  =  -L/sw/lib -lm -lfftw3 -L/usr/local/lib -lcfitsio -L/opt/local/lib -L$(HOME)/lib
  endif
#################################################################
# GFOTRAN path 
else
	FC	      = gfortran
	LINKER	  = gfortran
	
  	# cl> make multi=false
  # maybe remove this as it does  not compile without openmp anyway
  ifneq ($(multi),false)
  	MULTICORE = -fopenmp
  endif
  
  # array boundary check
  ifeq ($(debug),true)
    DEBUGGING =  -fcheck=all
  endif
  
  # -fno-whole-file deactivates some checks (e.g. explicit interface required for optional arguments)
  FLAG_ALL = -O3 -g -fbacktrace -finit-local-zero -ffixed-line-length-none $(MULTICORE) $(DEBUGGING)
    
  FFLAGS   = $(FLAG_ALL)
  LDFLAGS  = $(FLAG_ALL) -cpp Version.f

# add here library paths if required  
	ifeq ($(shell uname),Linux)
		LIBS   = -lm -lfftw3 -lcfitsio
	else
		LIBS   = -lm -lfftw3 -lcfitsio -L/opt/local/lib
	endif
endif


# use a suffix in file name (i.e. static, test etc.)
# cl> make name=test
ifneq ($(name),)
  SUFFIX = -$(name)
endif

# files to make
OBJS  = Modules.o \
				MCMax3D.o \
				InputOutput.o \
				Init.o \
				SetupStructure.o \
				Radiation.o \
				ComputePart.o \
				ReadParticleFits.o \
				RefIndData.o \
				KuruczData.o \
				SizeDis.o \
				RadiativeTransfer.o \
				Subroutines.o \
				writeFITS.o \
				Traveling.o \
				OutputMCobs.o \
				Raytrace.o \
				Convolution.o \
				Reddening.o \
				SpecialZone.o \
				OutputStats.o \
				ZoneInputOutput.o \
				EvolveStructure.o \
				delaunay_lmap_2d.o

# program name and install location
PROGRAM       = MCMax3D
DEST1	      = ${HOME}/bin
DEST2	      = mmin@michielmin.nl:/volume1/web/michielmin/resources

# make actions 
all:		version $(PROGRAM)
version:;	echo "#define gitversion \"$(shell git rev-parse HEAD)\"" > gitversion.h
clean:;		rm -f $(OBJS) $(PROGRAM)
install:	version $(PROGRAM)
			mv $(PROGRAM) $(DEST1)
test:		version $(PROGRAM)
			mv $(PROGRAM) $(DEST1)/MCMax3D_test
upload:		version $(PROGRAM)
			scp $(PROGRAM) $(DEST2)
			mv $(PROGRAM) $(DEST1)
echo:;		@echo $(SUFFIX)

# how to compile program 
.SUFFIXES : .o .f .f90 .F

.f.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

.F.o:
	$(FC) $(FFLAGS) -c $<

$(PROGRAM):     $(OBJS)
		$(LINKER) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

# recompile everything if Modules.f has changed 
$(OBJS):	Modules.f



