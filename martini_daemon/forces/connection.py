from .force import Force


class Connection(Force):
    # Dummy force that is only there for graphs and exclusions
    def _build(self):
        self._force_obj = None

    def add(self, members, params):
        i, j = members
        self._list.append((i, j))
        return self._interaction()

    def get_members(self, id):
        i, j = self._list[id]
        return [i, j]

    def update_params(self, id):
        raise NotImplementedError

    def build(self):
        pass

    def destroy(self):
        pass

    def is_instance(self, filter):
        return filter in {"bond", "connection"}
