

def replay(self, reactions: list[tuple[int, str, list[list[int]]]], until_frame=-1):
    """
    Replays reactions
    = runs the modification algorithm based on reaction names and atom indices
    """
    for frame, rx_name, frags in reactions:
        if 0 <= until_frame <= frame:
            break
        # find rx
        rx = self.top.rx_by_name(rx_name)
        assert rx is not None, f"{rx} could not be found."
        # find frags
        frags = [
            self.top.frag_by_name_and_atoms(frag_name, frag)
            for frag_name, frag in zip(rx.__reactants, frags)
        ]
        assert all(frag is not None for frag in frags)
        # modifications
        self.top.modification(frame, [(frags, rx)])
    self.system.reinitialize()
    # note to self: we don't make constraints and __vsites whole again as reactions can't modify those
    # so if they are made whole when the sim is constructed that's enough