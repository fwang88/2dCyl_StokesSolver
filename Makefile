all: stokes_test

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

#CLEANFILES = stokes
OBJECTS = main_test.o stokes.o io.o
#CLINKER = mpicc

stokes_test: $(OBJECTS)
	$(CLINKER) $(OBJECTS) -Ofast -g -o $@ ${PETSC_LIB}

#
#
#all: stokes
#
#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules
#
#CLEANFILES = stokes
#OBJECTS = main.o  InitialLevelSet.o RK2DReinit.o Euler_Reinit_2D.o AssembleB.o AssembleForce.o RK2DPeriod.o Euler_mex2D.o timestep.o
#CLINKER = mpicc
#
#stokes: $(OBJECTS)
#	$(CLINKER) $(OBJECTS) -Ofast -o $@ ${PETSC_LIB}
#	$(RM) -f $(OBJECTS)
#clean1:
#	-rm -f OBJECTS
