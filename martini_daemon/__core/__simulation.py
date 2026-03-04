
class Simulation:

    # user API
    def __init__(
        self,
        top_path: str,
        geom_path: str,
        sim_name: str,

    ):
        # TODO all default values should be in simulation only !!

        # the data based on which the simulation was constructed
        self.top_path = top_path
        self.geom_path = geom_path
        self.sim_name = sim_name

        # the "when" in the simulation
        self.simulation_step: int = 0
        self.simulation_time_ps: float = 0.

    def write_geometry(self, path: str) -> None:
        pass

    # programmer API

    # Custom logger - into the main log file / screen


    # Request trajectory, log file, text file, binary file, compressed binary file
    # - this enforce simulation name prefixes and backing up in a single place

    # event based architecture -- using not only reporters
    # TODO try to keep S* and T* free from reporter list and handle everything through simulation

    # Topology meta-info -- constraints, vsites, coupling, extra forces
