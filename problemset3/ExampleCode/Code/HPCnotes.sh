#___________________________________________________________________________________________________
# 1. LOGGING IN
# From computer
ssh netid@hpc.nyu.edu
# From HPC-Bastion host
ssh netid@mercer.es.its.nyu.edu
cd /scratch/$USER

#___________________________________________________________________________________________________
# 2. TEXT EDITTING
# A very versatile but steep learning curve text editor is vi. Another widely used one is emacs
vi textfile.txt
vi mfile.m
# Most improtant commands - :x (quit and save), :q (quit not save)

#___________________________________________________________________________________________________
# 3. LOADING PROGRAMS
# Check available modules
module avail
# Load a module, eg matlab
module load matlab

#___________________________________________________________________________________________________
# 4. SUBMITTING / MANAGING JOBS
# Start an itneractive session on a node (nodes,cores,time). Any job which isn't just moving files around
# or doing very very simple computations should be done using an interactive or batch processing session
qsub -I -q interactive -l nodes=1:ppn=1,walltime=01:00:00
	# Note: You can have one instance of Cygwin/Terminal open running an interactive session, running, say, Matlab
	#  		Another instance open running a separate session running an editor, say vi, editting your m-files etc
# Submit a batch processing job
qsub RUN.pbs
# Check on the status of all jobs
qstat -u $USER
# Delete a specific job
qdel <number here>
# Hold / release a specific job
qhld <number here>
qrls <number here>
# Delete all running jobs
qstat -u $USER | grep "R" | cut -d' ' -f1 | xargs qdel
# Delete all queued jobs
qstat -u $USER | grep "Q" | cut -d' ' -f1 | xargs qdel
# Gaze into the multicolor wonder of the HPC. Type :q to exit
pbstop
pbstop -u netid 

#___________________________________________________________________________________________________
# 5. BASIC SHELL COMMANDS
# Obviously many many more can be found online, these are the ones I use the most /
# Clear screen
clear
# List files in directory
ls -l
# Copy a file to a directory
cp file ./directory
# Check history to make sure you've done everything in the right order
history
# Create a new sub-directory
mkdir dirname
# Remove a file
rm file
# Remove a directory
rm -r directory
# Note: Use the -r flag whenever moving / copying / removing etc a directory

#________________________________________________________________________________________
# 6. TRANSFERRING FILES
# Example: Transfer RESULTS directory to PC 
# (-r used if transferring a directory, not used when just a file)
# Mercer to HPC (from the directory with RESULTS in it, while logged in to Mercer)
scp -r RESULTS $USER@hpc.nyu.edu:~/.
# HPC to PC (from PC)
scp -r netid@hpc.nyu.edu:~/RESULTS .

#_____________________________________________________________________________________________________
# 7. USING LOOPS TO SCALE CODE
# Shell-programming makes it incredibly easy to manage files / directories
# Here's an example that leads to running / queuing multiple jobs
# Example: 
# You have a set of *.code files that run your code in a directory /scratch/$USER/Project
# You have written a pbs-file to batch submit a job
# Now you want to run 20 instances of this code with a small change to each run 
	# Example: Different starting point for search as in Guvenen's example
# This small change should be representable just by the index of the run, eg run 1, run 2, ..., run 20
# For example in the example pbs file I provided, I would change variable1 to the run number, I would also change the
# name of the job so I can keep track of submitted jobs in the server queue
#######################################################################################################
# 0. Change to the base directory
cd /scratch/$USER/Project
# 1. Make directories
for x in $(seq 1 20); do mkdir run$x;done;
# 2. Copy all code to all directories
for x in $(seq 1 20); do cp *.code run$x;done;
# 3. Copy your base pbs-file to a seperate directory
mkdir RUNfiles
cp RUN.pbs RUNfiles
# 4. Make 20 instances of BASE.pbs
for x in $(seq 1 20); do cp RUNfiles/RUN.pbs RUNfiles/RUN_$x.pbs;done;
# 5. Manually edit each of the files so it knows what run it is (save each file)
for x in $(seq 1 20); do vi RUNfiles/RUN_$x.pbs ;done;
# 6. Distribute the pbs files
for x in $(seq 1 20)); do cp RUNfiles/RUN_$x.pbs run$x;done;
# 7. Submit all jobs
for x in $(seq 1 20); do cd /scratch/$USER/Project/run$x; qsub RUN_$x.pbs;done;
# Once this is set up the only code that should require editting is in the ../Project folder
# This can then copied out to all other sub-directories and ran from there

#_____________________________________________________________________________________________________
# 8. SEARCH AND REPLACE TEXT USING PERL
# This is really handy for editting a set of pbs-files
module load perl
perl -pi -e 's/texttoreplace/texttoreplaceit/g'
# In a loop over files (e.g. changing memory and walltime in *all* pbs files)
for x in $(seq 1 20); do perl -pi -e 's/ppn=10,mem=10GB,walltime=08:00:00/ppn=10,mem=40GB,walltime=10:00:00/g' RUN$x.pbs;done;

#___________________________________________________________________________________________________
# 9. CHANGES TO SHELL ENVIRONMENT AND ALIASES
# The environment you are operating in (the shell) can be editted, with new commands added
# This can make your life a lot easier when coding in the shell
# The file which needs to be editted is the bashrc file, to do this, when you log in to Mercer
vi .bashrc
# Any edits made to .bashrc will only take effect if saved, and then only upon exitting and re-entering Mercer
# Here is what my .bashrc looks like,
##################################################################################
# .bashrc
# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

# User specific aliases and functions

# Prompt color
export PS1="\e[0;32m[\u@\h \W]\$ \e[m"

# Aliases
alias ls='ls -l -hF --color=tty'
alias rm='rm -i'
alias ..='cd ..'
alias scratch='cd $SCRATCH'
alias home='cd ~/'
alias qstats='clear;qstat -u $USER;'
alias mercer='ssh sm4125@mercer.es.its.nyu.edu'
alias showqs='showq -u $USER'
alias HPC='cd $SCRATCH/HPC'
alias cls='clear;ls;'
##################################################################################
# The alias lines create new commands, shortcuts for things you do all the time
# E.g. cls clears the screen and lists all files
# E.g. ls replaces the standard ls with one that creates a vertical list with colors














