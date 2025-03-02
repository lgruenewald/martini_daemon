#!/usr/bin/bash

source /usr/local/gromacs-2021.7-double/bin/GMXRC
gmx_d grompp -f ../md.mdp -c system.gro -p system.top -o run.tpr
gmx_d mdrun -deffnm run -rerun system.gro -nt 1 -v
echo pot | gmx_d energy -f run.edr
echo 0 | gmx_d traj -f run.trr -s run.tpr -of forces
