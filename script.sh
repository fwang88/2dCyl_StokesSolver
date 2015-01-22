LISTTW="910"
LISTVEL="0.5"
LISTThigh="2200"
LISTTlow="1400"
dr=0.1
dz=0.1
maxz=400
maxr=10
pertb=0.0
lowtwidth=100
trestart=500
restart=1
initflag=0
temp_profile=1
tension=1.5

for Thigh in $LISTThigh
do
for Tlow in $LISTTlow
do
for tw in $LISTTW
do
for vel in $LISTVEL
do
outputstep=10
#outputstep=`echo "scale=4; $dr/$vel*20.0" | bc`
${PETSC_DIR}/${PETSC_ARCH}/bin/mpirun -np 4 ./stokes --r0=0.97 --muRat=20 --tension=$tension --maxR=$maxr --maxZ=$maxz --dr=$dr --dz=$dz --pertb=$pertb --Tlow=$Tlow --Thigh=$Thigh --Twidth=$tw --outputstep=$outputstep --vf=$vel --LowTWidth=$lowtwidth --trestart=$trestart --restart=$restart --initflag=$initflag --temp_profile=$temp_profile &
done
done
done
done

#initflag = 0: end; initflag = 1: periodic
