This document describes the additional directives that can be used for
creating reaction templates within martini_daemon.

# General recommendations

Martini_daemon specific things, such as [frag], [rx] directives, and the
product [moleculetype] directives (and all the directives within these
"parent" directives) should go to a separate file, currently recommended
with the ".rx" extension, which should be included in the main topology file,
after the .itp files.

If the topology files are being parsed by DaemonTopFile, the define DAEMON is
set, so putting the .rx files within an #ifdef DAEMON is possible, if the same
topology file is used in other software.


# Defining fragments - the [frag] directive

The first line in the [frag] directive should define the fragment name. The
fragment name should be unique (there should be no other [frag] or 
[moleculetype] with the same name). Fragment names can only contain ASCII 
letters, underscores and numbers. The [frag_*] subdirectives should
follow the [frag] directive directly for clarity.

## The [frag_atoms] directive

Optional. Must be after a [frag] directive was already defined. It applies
to the last [frag] directive. It is only valid to define this once per
[frag] directive.

Defines a list of atom types and name patterns to match against when
instantiating this fragment. If there is a mismatch, an error is raised.
Can be used for debugging templates, or to limit their misuse.

<!-- TODO: might be removed for simplicity if it is not used in practice -->

## The [frag_from] directive

Must be after a [frag] directive was already defined. If [frag_atoms] is present
for that [frag] directive, must come after [frag_atoms]. It applies to the
last [frag] directive.

Every line in this directive defines how this fragment should be instantiated.
The first column is the name of the parent fragment, and subsequent columns
contain the particle id's in the parent fragment that get mapped to this
fragment's particles. The number of particles in a fragment should be the same
for all lines, and same with the number of particles in [frag_atoms], if it
is defined.

All [frag_from] directives should come after the [frag] and [moleculetype]
directives they reference.

# Defining reactions - the [rx] directive

The first line in the [rx] directive should give the name of the reaction.
Reaction names can only contain ASCII letters, underscores and numbers.
The reaction name is only used for debugging at the moment, but in principle
they should be unique.

All [rx] directives should come after the [frag] and [moleculetype]
directives they reference, but before [system].

[rx_reactants], [rx_products] and [rx_conditions] are consdiered mandatory.
The directives within a [rx] directive should follow the [rx] directive
directly for clarity. Any of the [rx_*] subdirectives apply to the last
defined [rx] directive. It is only valid to define the [rx_*] subdirectives
once per [rx] directive, and they can come in any order.

## The [rx_reactants] directive

Contains the list of reactants for this reaction on a single line.
Can be [moleculetype] or [frag]. The numbering of the atoms for the product
and conditions will be the same as in the reactants. There can be 1 or 2 
reactants, other numbers of reactants are not supported.

## The [rx_products] directive

Contains the list of products for this reaction. Currently only 1 product
is supported, and it is recommended to make a new [moleculetype] entry
tailored to this specific [rx]. More flexibility might come later.

## The [rx_conditions] directive

List of reaction conditions to fulfill. See `docs/architecture.md`, section
Detection algorithm for details. One condition per line. Particles are indexed
in the following manner: first they index in reactant 1's atoms, then, if
there is a second reactant further indices index in reactant 2's atoms. So,
for example, if both reactants contain 3 particles, indices 1, 2 and 3 index
the first reactant and indices 4, 5 and 6 index the second reactant.

The following conditions are supported:

- `r_max`, arguments `i`, `j`, `r`, which rejects the reaction if the distance
  between particles `i` and `j` is above `r`.
- `r_min`, arguments `i`, `j`, `r`, which rejects the reaction if the distance
  between particles `i` and `j` is below `r`.
- `banned_angle`, arguments `i`, `j`, `k`, `min`, `max`, which rejects the
  reaction if the angle `ijk` (central atom `j`) falls between `min` and `max`.
  Both `min` and `max` are in degrees and the valid results for the angle
  calculation range from 0 to 180 degrees.
- `banned_dihedral`, arguments `i`, `j`, `k`, `l`, `min`, `max`, which rejects
  the reaction if the dihedral `ijkl` falls in the range `min` to `max`. Both
  `min` and `max` are in degrees. Valid results for the dihedral calculation
  range from 0 to 360 degrees (360 not included). If `min` is larger than `max`,
  it is assumed that it was meant to wrap around (for example, if `min` is 270
  and `max` is 90, it disallows it 0 to 90 and 270 to 360).
- `p`, argument `probability`, which rejects the reaction based on random
  chance. The `probability` is between 0.0 and 1.0, and represents the chance to
  accept the reaction (1.0 = always accept, 0.0 = always reject).
- `limit`, argument `count`, which rejects each reaction after `count` of this
  reaction have been accepted in this current reaction step.

## The [rx_break] and [rx_update] directives

Every line within these directives represents a group of particle indices to
be included in a break group. Each such group is used during a modification
algorithm to remove some interactions from reactants. Indices in these
directives work identically to how they work in [rx_conditions].

The difference between the two, is how interactions that are broken are
selected, see the section in `docs/architecture.md` for this for more
detail.

The main use case for [rx_break] is to specify a pair of particles, and to
break the molecule apart into two between those two particles.

The main use case for [rx_update] is to change how a part of the molecule
is constructed. It will break interactions within the group, but will not
break interactions that contain particles outside the group.

# Defining reaction products

Currently, it is recommended to write new [moleculetype] entries that define
a reaction product in the .rx file, specifically tailored with the reaction
template in mind. This is, because it is convenient to only define the
interactions that get formed during the reaction, while keeping most of
the interactions in the reactants that do not need updating. If an interaction
needs updating, they can be broken in the reaction template and readded by
the product.

## [atoms] directive in reaction products

An important difference between defining reaction products and just regular
[moleculetype] molecules is, that the [atoms] directive in reaction products 
does not represent new particles to be added to the system. The atom type,
charge and mass given here change the type, charge and mass of existing
particles. Other things, such as atom name, residue name, are currently
unchanged.

There are additional reaction product specific changes to the [atoms] directive
syntax in daemon:
- `*` is a valid atom type, which signifies that the atom type in the reactant
  should remain unchanged during the reaction. Note, that atom types are not
  patterns, only the `*` is recognized as a special value.
- Atom names are patterns that get checked against during the modification
  algorithm. If the pattern does not match an error is raised, this is
  present to help users catch and debug misbehaving reaction templates.

# Using non-daemon [moleculetype]s as reaction products

It is possible to use a predefined martini molecule as a reaction product.

To do this:
- add all reactant atoms to a single line within the [rx_update] directive
  for all reactants, so that all interactions that existed are broken first.
- add the name of the molecule to the [rx_product] directive

Defining new [moleculetype] directives is generally preferred, though, and
might be necessary, if the reactants and product contain constraints or
virtual sites, which cannot be broken or formed during reactions (also 
discussed in `docs/architecture.md`).

# Pattern matching

Pattern matching within the context of daemon is defined here. Since it is
mentioned in different later sections, it is described in a separate section.

A pattern matches a word if they are identical, but a few special characters
in patterns are also possible: `*` matches any string, `?` matches any single
character and curly braces (`{ab}`) match any of the characters enclosed
within the curly braces. A special character enclosed in curly braces matches
that special character literally.

Example use cases:
- `*` will match anything
- `TC*` will match any string starting with TC
- `C?` will match any string that contains C and one other character
- `C{123}` will match `C1`, `C2` and `C3`

Note: currently the unix fnmatch function is used to power this, which is also
used to match filenames, so it should be familiar.
