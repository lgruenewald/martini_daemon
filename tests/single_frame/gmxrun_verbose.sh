#!/usr/bin/bash

if [ -e respos.gro ]; then
  gmx_d grompp -f ../md.mdp -c system.gro -p system.top -r respos.gro -o run.tpr
else
  gmx_d grompp -f ../md.mdp -c system.gro -p system.top -o run.tpr
fi
gmx_d mdrun -deffnm run -rerun system.gro -nt 1 -v
echo pot | gmx_d energy -f run.edr
echo 0 | gmx_d traj -f run.trr -s run.tpr -of forces
