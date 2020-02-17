#!/bin/bash
#PBS -N XXX
#PBS -j XXX
#PBS -o XXX
#PBS -q XXX
#PBS -l nodes=5:ppn=YYY,mem=MMM
#PBS -r n

cd $PBS_O_WORKDIR

# do not forget to export
# mpi library

#---------------------------------------------------
# !!! HOSTLIST !!!
#     get hostnames of the nodes
#     and prepare host_list for run
cat ${PBS_NODEFILE} > all_hosts

#     replace YYY below, 
#     with the given ppn value above
#     i.e. if 
#          PBS -l nodes=5:ppn=20,mem=MMM
#          then 
#          /bin/sed -n '1~20p' all_hosts > host_list
/bin/sed -n '1~YYYp' all_hosts > host_list
#---------------------------------------------------

#---------------------------------------------------
# !!! RUN !!!

#     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     !!! WARNING                                !!!
#     !!! please prefer to                       !!!
#     !!! run this test on 5 nodes               !!!
#     !!!                                        !!!
#     !!! number of beads are 7, but             !!!
#     !!! frozen end.                            !!!
#     !!! number of beads to optimize is 5, thus !!!
#     !!! do not use more than 5 nodes           !!!
#     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#     replace NNN with number of nodes
#     replace YYY with ppn value

mpirun -np 5 --hostfile host_list lichem.MPI -n YYY -x react.xyz -c connect.inp -r regions.inp -o qsm_mpi.xyz -l qsm_mpi.log

#---------------------------------------------------

