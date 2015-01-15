LISTI="0.1"
LISTJ="0.02"
# 0.5 0.2 0.1 0.05 0.02 0.01"
for i in $LISTI
do
for j in $LISTJ 
do
${PETSC_DIR}/${PETSC_ARCH}/bin/mpirun -np 2 ./stokes --dr=$i --pertb=$j --dt=0.001 --outputstep=0.001 --maxZ=20 --dz=$i --maxR=10 --r0=0.97 --tension=1.0 --initflag=0
echo $i
done
done

#initflag = 0: end; initflag = 1: periodic