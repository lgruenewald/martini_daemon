from .force import Force


class Connection(Force):
    _members = 2

    def _set_force_obj(self):
        pass

    def _add_to_force_obj(self, params):
        pass

    def build(self):
        pass

    def destroy(self):
        return False

    def is_instance(self, filter):
        return filter in {"bond", "connection"}
