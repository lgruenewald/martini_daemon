#!/usr/bin/env python3

# ==== System layer ====
class ParticleList():
    """All particles in the system
    indices are called part_id
    values are (part_name: string, part_type: string, charge: float,
    mass: float)
    """
    _part_list: list[(str, str)]

    def __init__(self):
        self._part_list = []

    def add(self, part_name, part_type, charge, mass):
        """Adds a particle to the list, and returns its part_id"""
        self._part_list.append((part_name, part_type, charge, mass))
        return len(self._part_list) - 1

    def len(self):
        return len(self._part_list)


class BondList():
    """All bonds in the system
    indices are called bond_id
    values are (i: part_id, j: part_id, length: float, force: float)
    """
    _harmonic_bond_list: list[(int, int, float, float)]

    def __init__(self):
        self._harmonic_bond_list = []

    def add(self, part_id_i, part_id_j, length, force):
        """Adds a bond to the list, returns its bond_id"""
        self._harmonic_bond_list.append((part_id_i, part_id_j, length, force))
        return len(self._bond_list) - 1


class AngleList():
    """All the angles in the system
    indices are called angle_id
    values are (i: part_id, j: part_id, k: part_id, theta: float, force: float)
    """

    _harmonic_angle_list: list[(int, int, int, float, float)]

    def __init__(self):
        self._harmonic_angle_list = []

    def add(self, i, j, k, theta, force):
        """Adds an angle to the angle list, returns its angle_id"""
        self._harmonic_angle_list.append(
            (i, j, k, theta, force)
        )
        return len(self._harmonic_angle_list) - 1


class ProperDihedralList():
    """All the proper dihedrals in the system
    indices are called proper_dihedral_id
    values are (i: part_id, j: part_id, k: part_id, l: part_id,
     theta: float, force: float, multiplicity: int)
    """

    _dihedral_list: list[(int, int, int, int, float, float, int)]

    def __init__(self):
        self._dihedral_list = []

    def add(self, i, j, k, l, theta, force, multiplicity):
        """Adds a dihedral to the list, returns its proper_dihedral_id"""
        self._dihedral_list.append(
            (i, j, k, l, theta, force, multiplicity)
        )
        return len(self._dihedral_list) - 1


class ExclusionList():
    """All the exclusions in the system
    indices are caled excl_id
    values are (i: part_id, j: part_id)
    """

    _exclusion_list: list[(int, int)]

    def __init__(self):
        self._exclusion_list = []

    def add(self, i, j):
        self._exclusion_list.append((i, j))
        return len(self._exclusion_list) - 1


class ConstraintList():
    """All the constraints in the system
    indices are called constraint_id
    values are (i: part_id, j: part_id, length: float)
    """

    _constraint_list: list[(int, int, float)]

    def __init__(self):
        self._constraint_list = []

    def add(self, i, j, length):
        """Adds a constraint to the list, returns its constraint_id"""
        self._constraint_list.append((i, j, length))
        return len(self._constraint_list) - 1


class ImproperDihedralList():
    """All the improper dihedrals in the system
    indices are called improper_dihedral_id
    values are (i: part_id, j: part_id, k: part_id, l: part_id,
     theta: float, force: float)
    """

    _dihedral_list: list[(int, int, int, int, float, float)]

    def __init__(self):
        self._dihedral_list = []

    def add(self, i, j, k, l, theta, force):
        """Adds a dihedral to the angle list, returns its improper_dihedral_id
        """
        self._dihedral_list.append(
            (i, j, k, l, theta, force)
        )
        return len(self._dihedral_list) - 1


class SysStar():
    def __init__(self):
        self.particles = ParticleList()
        self.bonds = BondList()
        self.angles = AngleList()
        self.dihedrals = ProperDihedralList()
        self.impropers = ImproperDihedralList()
        self.exclusions = ExclusionList()
        self.constraints = ConstraintList()

    def dump(self):
        print("Particles:", self.particles._part_list)
        print("Bonds:", self.bonds._harmonic_bond_list)
        print("Angles:", self.angles._harmonic_angle_list)
        print("Dihedrals:", self.dihedrals._dihedral_list)
        print("Impropers:", self.impropers._dihedral_list)
        print("Exclusions:", self.exclusions._exclusion_list)
        print("Constraints:", self.constraints._constraint_list)
