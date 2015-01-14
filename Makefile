all: stokes

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

CLEANFILES = stokes
OBJECTS = main03.o  InitialLevelSet.o RK2DReinit.o Euler_Reinit_2D.o AssembleB.o AssembleForce.o RK2DPeriod.o Euler_mex2D.o timestep.o
#CLINKER = mpicc

stokes: $(OBJECTS)
	$(CLINKER) $(OBJECTS) -Ofast -o $@ ${PETSC_LIB}