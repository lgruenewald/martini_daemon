from .force import Force
import openmm as mm


class Constraint(Force):
    _members = 2

    def _add_to_force_obj(self):
        raise Exception("Can't add constraints after starting the run due"
                        " to periodic boundary conditions")

    def remove(self, i):
        raise Exception("Can't remove constraints")

    def build(self):
        if self._rebuild:
            for (i, j, length) in self._list:
                self._sysstar._system.addConstraint(i, j, length)
            self._rebuild = False

    _destroyable = False

    def destroy(self):
        raise Exception("Can't destroy constraints.")

    _filters = {"bond", "constraint"}

    def set_force_group(self, fg):
        return False

    def constraints_to_harmonic_bonds(self, harmonic=False):
        if harmonic:
            self.harmonic = mm.HarmonicBondForce()
            for i in range(len(self._list)):
                self._sysstar._system.removeConstraint(0)
                i, j, r = self._list[i]
                self.harmonic.addBond(i, j, r, 10000)
            # this sets reinitialize to True
            self._sysstar.add_force(self.harmonic)
        else:
            # we assume that only a small minimization has taken place
            # and thus no constraints were broken across pbc
            self._sysstar.remove_force(self.harmonic)
            self.harmonic = 0
            self._rebuild = True
            self.build()
            self._sysstar._reinitialize = True
