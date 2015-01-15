LISTTW="910"
LISTVEL="1"
LISTThigh="1800"
LISTTlow="1400"
dr=0.1
dz=0.1
maxz=100
maxr=10
pertb=0.02
lowtwidth=100
trestart=0
restart=0
initflag=0
temp_profile=0

for Thigh in $LISTThigh
do
for Tlow in $LISTTlow
do
for tw in $LISTTW
do
for vel in $LISTVEL
do
outputstep=`echo "scale=4; $dr/$vel*20.0" | bc`
${PETSC_DIR}/${PETSC_ARCH}/bin/mpirun -np 4 ./stokes --maxR=$maxr --maxZ=$maxz --dr=$dr --dz=$dz --r0=0.97 --tension=1.0 --Twidth=$tw --Tlow=$Tlow --Thigh=$Thigh --vf=$vel --outputstep=$outputstep --pertb=$pertb --LowTWidth=$lowtwidth --trestart=$trestart --restart=$restart --initflag=0 --temp_profile=$temp_profile
done
done
done
done

#initflag = 0: end; initflag = 1: periodic
