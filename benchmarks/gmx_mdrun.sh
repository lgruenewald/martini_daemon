#!/usr/bin/bash

gmx mdrun -deffnm run -nt 10 -pin on -v >>gro.log &>> gro.log
