
# Master Makefile for making the executable of the generalized biological model.

.SUFFIXES:
#.SUFFIXES: .F
.SUFFIXES: .o .f90 .F .F90 
	SHELL	= /bin/sh
	DEF_FLAGS     = -P -C -traditional 
#	NETCDFINCDIR	= /usr/local/include
# 	ETCDFLIB       = /usr/local/lib -lnetcdf

         CPP      = /usr/bin/cpp 
         CPPFLAGS = $(DEF_FLAGS) #-DINTEL   
#         FC    = ifort -checkall
         FC       = mpif90 #-static-libcxa -i-static

BIOMODULE	= \
libbiomodule.a(mod_1D.o)			\
libbiomodule.a(mod_phytoplankton.o)		\
libbiomodule.a(mod_zooplankton.o)		\
libbiomodule.a(mod_bacteria.o)			\
libbiomodule.a(mod_detritus.o)			\
libbiomodule.a(mod_DOM.o)			\
libbiomodule.a(mod_nutrient.o)			\
libbiomodule.a(phytoplankton.o)			\
libbiomodule.a(zooplankton.o)			\
libbiomodule.a(bacteria.o)			\
libbiomodule.a(detritus.o)			\
libbiomodule.a(DOM.o)				\
libbiomodule.a(nutrient.o)			\
libbiomodule.a(get_parameter.o)			\
libbiomodule.a(bio_mixing.o)

LIBS	=	libbiomodule.a

libbiomodule.a: $(BIOMODULE)


clean:
	-rm -f lib*.a  *.mod *.o *.f90

.SUFFIXES: .o .f90 .F .F90 
.F.o:
	$(CPP) $(CPPFLAGS) $*.F > $*.f90
	$(FC) -c $(FFLAGS) $*.f90 > $*.o
	rm $*.f90
