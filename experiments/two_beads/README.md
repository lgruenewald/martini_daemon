_These examples were initially hosted at <https://git.sr.ht/~ma3ke/openmm_breaking_bonds>._

# Breaking the bond between two beads in OpenMM

In this example program, a two-bead system is set up with arbitrary parameters.
At the start, the beads start out with a harmonic bond potential between them,
until this bond is broken at the halfway point of the simulation.

There are two example programs here. One is written in C++ and the other in
Python. Both files are very much mirrored.

This example serves to demonstrate how this can be accomplished with the C++
interface to the [OpenMM][openmm] molecular mechanics library. By no means do I
intend to suggest that this is the only way. The same may be very feasible from
the Python interface as well, or perhaps a different route of changing the
bonds definitions is possible after all.

## Prerequisites

I wrote and ran this on Linux.

### Python version (`two_beads.py`)

Make sure that the `openmm` Python package is installed in your current
environment.

### C++ version (`two_beads.cpp`)

OpenMM must be installed somewhere on the system, for instance in
`/usr/local/openmm`. Here, I will assume it is installed at that location, but
it is fine to substitute it with your own installation location where relevant.
[Installing OpenMM from source][install] is pretty straightforward and will do
this by default. Make sure that _doxygen_ and _swig_ are installed so the
installation proceeds.

To produce the distance graph from the final trajectory, [Gromacs][gromacs] 
commands are used.
The graph can be displayed with _xmgrace_ (or [_phrace_][phrace]! :).

## Compile (C++)

Before we can execute this example, we need to compile it.

```console
g++ -I/usr/local/openmm/include two_beads.cpp -L/usr/local/openmm/lib -lOpenMM -o two_beads
```

(Or `clang++` for that matter. They should behave interchangibly in this case.)

## Run

```console
# C++ version
./two_beads > out.gro 

# Python version
python two_beads.py > out.gro
```

Since the program prints the trajectory as gro frames to _stdout_, we redirect
that stream into `out.gro`. This file contains the trajectory.

Quite possibly when you try to run the C++ version, a problem with loading the
OpenMM library will occur. In this case, you may need to append the path to 
the OpenMM library to the `LD_LIBRARY_PATH` environment variable. You can do
this as follows:

```console
# bash, etc.
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/openmm/lib # Or wherever.
# fish
set -Ux LD_LIBRARY_PATH /usr/local/openmm/lib # Or wherever.
```

## Visualize and analyze

What we expect to see is that the beads will move around a fixed distance in
the first half of the trajectory since they are bonded, and will move away from
each other when that bond is broken. Since the beads both have a positive
charge, they should repel each other to a distance greater than the bond
length, once the bond is broken.

Let's study our result.

First, we can look at the trajectory. For this, you can use, for example, _vmd_
or _pymol_.

Second, we can study the average distance between the beads over the
trajectory. What we expect to see is a rise in that distance at the midpoint of
the time series.

We can use `gmx distance` to calculate the distances.

```console
gmx distance -s out.gro -f out.gro -oav distave.xvg
```

If we inspect the graph, we will see something like the following.

```
                                Average distance                                
                                                                          ░▒██░ 
                                                                       ██░▒   ░▒
                                                                  ░▒▒█▒  ░      
D                                                       █░█▒▒▒ ░▒▒▒ ░           
i                                                   ░██▒ ░ ░ ░▒▒░               
s                                                   ░                           
t                                               █▒██                            
a                                             ░▒                                
n                                             ▒                                 
c                                           ▒█                                  
e                             ░▒░▒        ░▒                                    
       ░▒█▒ ░█▒█▒█▒ █▒▒█░░░█▒▒▒░░░░ ░░█▒  ░░                                    
(   ░░▒▒  ░▒▒     ░▒  ░ ░▒░ ░     ░█▒░  ██                                      
n   ░▒                                                                          
m  ░░                                                                           
)  ░                                                                            
  ░                                                                             
                                                                                
                                                                                
  ░                                                                             
                                  Time (ps)                                   
Summary:  202 items,  mean ± σ  2.2327778 ± 0.7983015,  min … max  0.2 … 3.654
```

(This graph was produced using _phrace_: `phrace distave.xvg -w 80 -h 25`.)

[openmm]: https://openmm.org/
[install]: http://docs.openmm.org/latest/userguide/library/02_compiling.html
[gromacs]: https://www.gromacs.org/
[phrace]: https://git.sr.ht/~ma3ke/phrace

---

Written by Marieke Westendorp in collaboration with Linus Grünewald, October 2024.
