You are running Q-Chem version: 5.0.2

#
# job setting
#
local host:  nid04204
current dir: /scratch1/scratchdirs/ewcss/with_water
input file:  /scratch1/scratchdirs/ewcss/with_water/243_2674637/rct_1_471171.in
output file: /scratch1/scratchdirs/ewcss/with_water/243_2674637/rct_1_471171.out
nprocs     : 1
nthreads   : 24
#
# qchem installation setting
#
QC:          /global/common/edison/software/qchem/5.1beta
QCAUX:       /global/common/edison/software/qchem/5.1beta/qcaux
QCPROG:      /global/common/edison/software/qchem/5.1beta/exe/qcprog.exe
QCPROG_S:    /global/common/edison/software/qchem/5.1beta/exe/qcprog.exe
PARALLEL:    -DPARALLEL
QCMPI:       slurm
#
# qchem directory setting
#
qcrun:       qchem43027
QCSCRATCH:   /dev/shm/qcscratch/
QCLOCALSCR:  /tmp
QCTMPDIR:    /tmp
QCFILEPREF:  /tmp/qchem43027
QCSAVEDIR:   
workdirs:    /tmp/qchem43027
workdir0:    /tmp/qchem43027
partmpdirs =  
#
# parallel setting
#
QCRSH:           ssh
QCMPI:           slurm
QCMPIRUN:        srun
QCMACHINEFILE:   

#
# env setting
#
exported envs:   QC QCAUX QCSCRATCH QCRUNNAME QCFILEPREF QCPROG QCPROG_S GUIFILE
remove work dirs /tmp/qchem43027.0 -- /tmp/qchem43027.0
rm -rf /tmp/qchem43027.0
rm -rf /tmp/qchem43027
