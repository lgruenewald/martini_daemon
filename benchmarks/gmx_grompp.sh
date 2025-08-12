#!/usr/bin/bash

# args: topol, gro, nsteps, pcoupl_type, pcoupl_pressure

sed "s/XSTEPS/$3/;" ../inputs/template.mdp > run.mdp
gmx grompp -f run.mdp -c $2 -p $1 -o run.tpr >>gro.log &>> gro.log
