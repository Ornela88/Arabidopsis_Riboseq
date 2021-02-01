#!/bin/bash
#!
#! Example SLURM job script for Darwin (Sandy Bridge, ConnectX3)
#! Last updated: Sat Apr 18 13:05:53 BST 2015
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J Map59
#! Which project should be charged:
#SBATCH -A CHIARUGI-CCLD-SL2-CPU
#! How many whole nodes should be allocated?
#SBATCH --nodes=2
#! How many (MPI) tasks will there be in total? (<= nodes*16)
#SBATCH --ntasks=30
#! How much wallclock time will be required?
#SBATCH --time=5:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ornela.maloku@gmail.com
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! Do not change:
#SBATCH -p clincloud-himem

#! sbatch directives end here (put any additional directives above this line)

#! Notes:
#! Charging is determined by core number*walltime.
#! The --ntasks value refers to the number of tasks to be launched by SLURM only. This
#! usually equates to the number of MPI tasks launched. Reduce this from nodes*16 if
#! demanded by memory requirements, or if OMP_NUM_THREADS>1.
#! Each task is allocated 1 core by default, and each core is allocated 3994MB. If this
#! is insufficient, also specify --cpus-per-task and/or --mem (the latter specifies
#! MB per node).

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#! ############################################################
#! Modify the settings below to specify the application's environment, location 
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load default-impi                   # REQUIRED - loads the basic environment

#! Insert additional module load commands after this line if needed:
module load star/2.5.0a
module load subread-1.6.2-gcc-5.4.0-7zywp5u

#! Full path to application executable: 
application=""

#! Run options for the application:
options=""

#! Work directory (i.e. where the job will run):
:

cd /rds/project/dc702/rds-dc702-bio2_core_rds/Personal_folders_ext/hpcmalo1/

STAR --runMode genomeGenerate --runThreadN 6 --genomeDir star_index --genomeFastaFiles Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --sjdbGTFfile Arabidopsis_thaliana.TAIR10.45.gtf
