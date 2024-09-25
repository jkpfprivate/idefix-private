#!/bin/bash

module purge
module load intelmpi/18.2
module load fftw/3.3.8
module load python/3.6.8
module load gcc/9.1.0
module load texlive/2021
module list

source activate cartopy

lmin=0
lmax=17

for (( l=$lmin; l<$lmax; l++ ))
do
  for (( m=0; m<=$l; m++ ))
  do
    if ! test -f "python/sliceplot_nofcb/MOLL_Ylm"$l$m".png"; then
      echo "sliceplot_nofcb/MOLL_Ylm"$l$m".png does not exist"
#    else
#      echo "sliceplot_nofcb/MOLL_Ylm"$l$m".png exists"
    fi
    if ! test -f "python/sliceplot_nofcb/MOLL_Slm"$l$m"_th.png"; then
      echo "sliceplot_nofcb/MOLL_Slm"$l$m"_th.png does not exist"
    fi
    if ! test -f "python/sliceplot_nofcb/MOLL_Slm"$l$m"_phi.png"; then
      echo "sliceplot_nofcb/MOLL_Slm"$l$m"_phi.png does not exist"
    fi
    if ! test -f "python/sliceplot_nofcb/MOLL_Tlm"$l$m"_th.png"; then
      echo "sliceplot_nofcb/MOLL_Tlm"$l$m"_th.png does not exist"
    fi
    if ! test -f "python/sliceplot_nofcb/MOLL_Tlm"$l$m"_phi.png"; then
      echo "sliceplot_nofcb/MOLL_Tlm"$l$m"_phi.png does not exist"
    fi
  done
done

