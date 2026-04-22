#!/usr/bin/bash

if [ -e respos.gro ]; then
  gmx_d grompp -f ../md.mdp -c system.gro -p gromacs.top -o run.tpr -r respos.gro -maxwarn 1
else
  gmx_d grompp -f ../md.mdp -c system.gro -p gromacs.top -o run.tpr -maxwarn 1
fi
gmx_d mdrun -deffnm run -rerun system.gro -nt 1 -v
echo pot | gmx_d energy -f run.edr
echo 0 | gmx_d traj -f run.trr -s run.tpr -of forces
rm mdout.mdp
rm log.txt
rm run*
rm \#*
