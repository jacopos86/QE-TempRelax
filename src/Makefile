# Makefile for TempRelax calculations

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../dft-d3/ \
          $(MOD_FLAG)../../LR_Modules/ \
	  $(MOD_FLAG)../../PHonon/PH/ \
	  $(MOD_FLAG)../../PP/src

TRELOBJS = \
consts.o \
temp_rel_vars.o \
vscr_mod.o \
screening_vars.o \
apaw_pot.o \
ph_lt.o \

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a \
         ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a \
         ../../dft-d3/libdftd3qe.a
PPOBJS = ../../PP/src/libpp.a
PHOBJS = ../../PHonon/PH/libph.a ../../PHonon/PH/libphaux.a

MODULES = $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PPOBJS) $(TRELOBJS)

TLDEPS= trelibs

all : tldeps libtrel.a temp_relax.x\

libtrel.a : $(TRELOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

qexml.x: qexml.o qexml_example.o
	$(LD) $(LDFLAGS) -o $@ qexml_example.o qexml.o $(QEOBJS) $(LIBOBJS)

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

temp_relax.x : temp_relax.o libgw.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		temp_relax.o libtrel.a $(MODULES) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../TempRelax/src/$@ . )

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L libtrel.a

include make.depend
