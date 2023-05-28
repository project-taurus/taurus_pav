#!/bin/bash

# For debugging purposes
#set -x

#################
#  Environment  #
#################

inter="usdb"  

omp="no"
nthreads=4

mpi="no"
nprocess=4

heredir=$(pwd)
heredir=$(pwd)
outdir=$heredir/out
wrkdir=$heredir/wrk
exedir=$heredir/../../exe
auxdir=$heredir/data
intdir=$heredir/../../int

code=taurus_pav.exe

if [ ! -d $outdir ]; then mkdir $outdir; fi 
if [ ! -d $wrkdir ]; then mkdir $wrkdir; fi 

#################
#  Calculation  #
#################

cd $wrkdir 

cp $exedir/$code .
cp $auxdir/template_input.txt input.txt
cp $auxdir/wavefunction_mg25_usdb_minimum.bin left_wf.bin
cp $auxdir/wavefunction_mg25_usdb_minimum.bin right_wf.bin
cp $intdir/$inter.sho .

# Running the script      
# For OpenMP calculations
if [ $omp = "yes" ] && [ $nthreads -ne 0 ]; then
  export OMP_NUM_THREADS=$nthreads
fi

# Runs the code (MPI or not)
if [ $mpi = "yes" ]; then
  mpirun -np $nprocess ./$code < input.txt > results
else
  ./$code < input.txt > results
fi 
    
echo "The results can be found in the files: results, out/*txt"

# Clean up
rm -f $code $inter.sho $inter.red input.txt left_wf.bin right_wf.bin

mv results $heredir/
mv *bin $outdir/
