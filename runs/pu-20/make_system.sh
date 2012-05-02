#!/bin/bash

python ../../ibi/bead_spring_system.py  \
   --num_chains 20 \
   --block HHSSSSSSSSSSSSSS \
   --num-blocks 14 \
   --filename coarse-system.lammps \
   --bond_length "SS=4.8666,HH=11.099,HS=9.245" \
   --density 1.0698 \
   --bead_masses "S=72.10776,H=253.261245" 


#bond SS r0 4.866691   k 7.193422
#bond HH r0 11.099044  k 1.056467
#bond HS r0 9.2456160311    k 5.34556100981

# bead masses
#mass S 72.10776
#mass H 253.261245

# known pairs for mixed.
#pair SS pair.table.SS
#pair HH pair.table.HH
#rdf files


