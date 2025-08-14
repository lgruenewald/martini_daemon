from .force import Force


class Connection(Force):
    _members = 2

    def _set_force_obj(self):
        pass

    def _add_to_force_obj(self, params):
        pass

    def build(self):
        pass

    _destroyable = False
    
    def destroy(self):
        return False

    _filters = {"bond", "connection"}
