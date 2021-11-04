#!/bin/bash
#PBS -V
#PBS -S /bin/tcsh
#
#### Job name
#PBS -N run_01
#
#### Nodes,Cores,Time
#PBS -l nodes=1:ppn=10,walltime=08:00:00
#
#### Mail Options (email me on: abort,begin,end)
#PBS -M <your netid>@nyu.edu
#PBS -m abe
#
#### Redirect error and output files
#PBS -e localhost:/scratch/<your netid>/HPC/${PBS_JOBNAME}.err${PBS_JOBID}
#PBS -o localhost:/scratch/<your netid>/HPC/${PBS_JOBNAME}.out${PBS_JOBID}

#### load matlab module
module load matlab

#### start MATLAB, everything between EOF's is ran in matlab
#########################################################################
matlab -nodisplay -nodesktop -nosplash << EOF
% Matlab Code
%___________________________________________________________________
variable1 	= 1;
nworkers    = 10;
%___________________________________________________________
parpool(nworkers);
addpath(genpath('/scratch/<your netid>/COMPECON'));
start(variable1);
matlabpool close
EOF
#########################################################################
### MATLAB closed

#### Compile pdf of figures and tables
latex RESULTS.tex
dvips RESULTS.dvi
ps2pdf RESULTS.ps

echo ""
echo "Done at "`date`
