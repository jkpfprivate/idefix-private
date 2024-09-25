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
lmax=4

for (( l=$lmin; l<$lmax; l++ ))
do
  for (( m=0; m<=$l; m++ ))
  do
    echo "Hi $l$m"
    sed -i "67s/.*/d.Vc(VX1,k,j,i) = d.Ylm_r("$l","$m",k,j);/"    setup.cpp
    sed -i "68s/.*/d.Vc(VX2,k,j,i) = d.Slm_th("$l","$m",k,j);/"   setup.cpp
    sed -i "69s/.*/d.Vc(VX3,k,j,i) = d.Slm_phi("$l","$m",k,j);/"  setup.cpp
    make -j 8
    ./idefix
    cp data.0000.vtk "data.0000.Sl"$l"m"$m".vtk"
    cd python
    python3 slicePlotTmp.py -vtk 0 -var VX1 VX2 VX3 -r -proj moll -save
    mv sliceplot_nofcb/RAD_MOLL_VX1_000.png "sliceplot_nofcb/MOLL_Ylm"$l$m".png"
    mv sliceplot_nofcb/RAD_MOLL_VX2_000.png "sliceplot_nofcb/MOLL_Slm"$l$m"_th.png"
    mv sliceplot_nofcb/RAD_MOLL_VX3_000.png "sliceplot_nofcb/MOLL_Slm"$l$m"_phi.png"
    cd ..
    sed -i "68s/.*/d.Vc(VX2,k,j,i) = d.Tlm_th("$l","$m",k,j);/"   setup.cpp
    sed -i "69s/.*/d.Vc(VX3,k,j,i) = d.Tlm_phi("$l","$m",k,j);/"  setup.cpp
    make -j 8
    ./idefix
    cp data.0000.vtk "data.0000.Tlm"$l$m".vtk"
    cd python
    python3 slicePlotTmp.py -vtk 0 -var VX2 VX3 -r -proj moll -save
    mv sliceplot_nofcb/RAD_MOLL_VX2_000.png "sliceplot_nofcb/MOLL_Tlm"$l$m"_th.png"
    mv sliceplot_nofcb/RAD_MOLL_VX3_000.png "sliceplot_nofcb/MOLL_Tlm"$l$m"_phi.png"
    cd ..
  done
done

