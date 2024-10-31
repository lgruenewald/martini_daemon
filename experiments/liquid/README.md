In this example, a box with n particles is generated. Initially, there are no
bonds in the system, but if two particles come within a cutoff distance, a
harmonic bond is formed.

Run using
python liquid.py > traj.xyz

Visualize using Pymol
pymol traj.xyz vis.pml
