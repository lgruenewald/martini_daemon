#!/usr/bin/bash

gmx mdrun -deffnm run -nt 10 -v >>gro.log &>> gro.log
