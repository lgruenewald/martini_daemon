In this example, a box with n particles is generated. Initially, there are no
bonds in the system, but if two particles come within a cutoff distance, a
harmonic bond is formed. Every particle can form bonds with 2 atoms (not just
one like in liquid.py). Bond angles are set to 180 degrees, with a weak force.

Run using
python chain.py > traj.gro

Visualize using
vmd traj.gro
