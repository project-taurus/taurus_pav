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
cp $intdir/$inter.sho .

i=0

for deformation1 in beta_0.250_gamma_20 beta_0.280_gamma_12 beta_0.300_gamma_05
do

  let i=$i+1
  j=0

  for deformation2 in beta_0.250_gamma_20 beta_0.280_gamma_12 beta_0.300_gamma_05
  do

    let j=$j+1
    if [ $j -lt $i ]; then continue; fi

    echo "Left state: $deformation1 --- Right state: $deformation2"

    fileref1=wavefunction_mg24_usdb_${deformation1}
    fileref2=wavefunction_mg24_usdb_${deformation2}
 
    fileelm=matelem_mg25_usdb_LEFT_${deformation1}_RIGHT_${deformation2}
  
    cp $auxdir/$fileref1.bin left_wf.bin
    cp $auxdir/$fileref2.bin right_wf.bin
  
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
  
    mv results $heredir/results_$fileelm 
    mv projmatelem_states.bin $outdir/${fileelm}_projmatelem_states.bin
    mv projmatelem_M1.bin     $outdir/${fileelm}_projmatelem_M1.bin
    mv projmatelem_E2.bin     $outdir/${fileelm}_projmatelem_E2.bin
      
    rm -f left_wf.bin right_wf.bin
     
  done
done
    
echo "The results can be found in the files: results_*, out/*txt"

# Clean up
rm -f $code $inter.sho $inter.red input.txt 
