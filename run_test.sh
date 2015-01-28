maxr=10
maxz=100
dr=0.1
dz=0.1
r0=0.97
tension=6.0
pertb=0.0
period_or_end=0
mui=0.91
muo=1.0
vf=0.5
temp_profile=0
tlow=1400
thigh=2200
twidth=910
lowtwidth=100
restart=1
trestart=226
outputdt=0
nghostlayer=3
epsilon=0.3
reinitstep=4

outputdt=`echo "scale=4; $dr/$vf*5.0" | bc`

${PETSC_DIR}/${PETSC_ARCH}/bin/mpirun -np 1 ./stokes_test --maxr=$maxr --maxz=$maxz --dr=$dr --dz=$dz --r0=$r0 --tension=$tension --pertb=$pertb --period_or_end=$period_or_end --mui=$mui --muo=$muo --vf=$vf --temp_profile=$temp_profile --tlow=$tlow --thigh=$thigh --twidth=$twidth --lowtwidth=$lowtwidth --restart=$restart --trestart=$trestart --outputdt=$outputdt --nghostlayer=3 --epsilon=0.3 --reinitstep=$reinitstep &
