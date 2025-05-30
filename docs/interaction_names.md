List of filters useable for filtering interactions in the graph language.
Also doubles as a list of interactions implemented in martini_daemon.

## Bonds

- `bond` - any bond listed here
- `harmonic_bond` - bond type 1 and 6
- `g96_bond` - bond type 2
- `morse_bond` - bond type 3
- `cubic_bond` - bond type 4
- `connection` - bond type 5
- `fene_bond` - bond type 7
- `distance_restraint` - bond type 10
- `constraint` - constraint type 1 and 2

## Angles

- `angle` - any angle listed here
- `harmonic_angle` - angle type 1
- `g96_angle` - angle type 2
- `cross_bond_bond` - angle type 3
- `cross_bond_angle` - angle type 4
- `urey_bradley` - angle type 5
- `quartic_angle` - angle type 6
- `linear_angle` - angle type 9
- `restricted_angle` - angle type 10

## Dihedrals

- `dihedral` - any dihedral listed here
- `proper_dihedral` - dihedral type 1 and 9
- `improper_dihedral` - dihedral type 2
- `rb_torsion` - dihedral type 3 and 5
- `restricted_dihedral` - dihedral type 10
- `combined_bending_torsion` - dihedral type 11

## Virtual sites

- `virtual_site` - any virtual site, matches the virtual particle as well as constructing particles
- `vsite` - alias for virtual_site
- `vsite1` - virtual_sites1 type 1
- `vsite2` - virtual_sites2 type 1
- `2fd` - virtual_sites2 type 2
- `vsite3` - virtual_sites3 type 1
- `3fd` - virtual_sites3 type 2
- `3fad` - virtual_sites3 type 3
- `3out` - virtual_sites3 type 4
- `4fdn` - virtual_sites4 type 2
- `com` - virtual_sitesn type 2
- `weighed_average` - virtual_sitesn type 1 and 3

## Other

- `pair` - only matches pairs
