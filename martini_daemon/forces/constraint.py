from .force import Force, Interaction
import openmm as mm  # type: ignore[import-untyped]


class Constraint(Force):

    _list: list
    _rebuild: bool
    # or when removing elements from the bond

    def __init__(self, sysstar):
        self._list = []
        self._sysstar = sysstar
        self._rebuild = True

    def _build(self):
        for (i, j, length) in filter(None, self._list):
            self._sysstar._system.addConstraint(i, j, length)

    def add(self, members, params):
        if not self._rebuild:
            raise Exception("Can't add constraints after starting the run due"
                            " to periodic boundary conditions")
        i, j = members
        length = params[0]
        self._list.append((i, j, length))
        return self._interaction()

    def get_members(self, id):
        i, j, _ = self._list[id]
        return [i, j]

    def update_params(self, i, *params):
        raise NotImplementedError

    def remove(self, i):
        raise Exception("Can't remove constraints")

    def build(self):
        if self._rebuild:
            self._build()
            self._rebuild = False

    def destroy(self):
        # invalid op for constraints
        raise Exception("Can't destroy.")

    def _interaction(self):
        return Interaction(self, len(self._list) - 1)

    def is_instance(self, filter):
        return filter in {"bond", "constraint"}
