from .force import Force


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
