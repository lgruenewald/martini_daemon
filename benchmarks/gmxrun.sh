#!/usr/bin/bash

source /usr/local/gromacs-2021.7/bin/GMXRC
sed "s/XSTEPS/$3/" ../template.mdp > run.mdp
gmx grompp -f run.mdp -c $2 -p $1 -o run.tpr >>gro.log &>> gro.log
gmx mdrun -deffnm run -nt 10 -v >>gro.log &>> gro.log
rm mdout.mdp
rm run*
