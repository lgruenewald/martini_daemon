#!/usr/bin/bash

source /usr/local/gromacs-2024.1-double/bin/GMXRC
if [ -e respos.gro ]; then
  gmx_d grompp -f ../md.mdp -c system.gro -p system.top -o run.tpr -r respos.gro >>log.txt &>> log.txt
else
  gmx_d grompp -f ../md.mdp -c system.gro -p system.top -o run.tpr >>log.txt &>> log.txt
fi
gmx_d mdrun -deffnm run -rerun system.gro -nt 1 -v >>log.txt &>> log.txt
echo pot | gmx_d energy -f run.edr >>log.txt &>> log.txt
echo 0 | gmx_d traj -f run.trr -s run.tpr -of forces >>log.txt &>>log.txt
rm mdout.mdp
rm log.txt
rm run*
rm \#*
