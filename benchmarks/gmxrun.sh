#!/usr/bin/bash

# args: topol, gro, nsteps, pcoupl_type, pcoupl_pressure

source /usr/local/gromacs-2023.4/bin/GMXRC
sed "s/XSTEPS/$3/;s/PCOUPL/$4/;s/REFP/$5/;" ../template.mdp > run.mdp
gmx grompp -f run.mdp -c $2 -p $1 -o run.tpr >>gro.log &>> gro.log
gmx mdrun -deffnm run -nt 10 -v >>gro.log &>> gro.log
