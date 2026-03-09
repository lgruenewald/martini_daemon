
class Interchange:
    """
    A class that handles some of the things "Simulation" and the base "Reporter" class used to.

    The goal of this class is to:
    * hold simulation metadata, such as current step, simulation name.
    * own all open file handles and loggers.
    * provide a friendly interface for requesting output files.
    * be passed around to all reporters to provide the required metadata for reporting.
    * own all reporters, and provide an API for calling them all whenever a certain event occurs.

    Simulation should inherit this class.
    """

    def write_geometry(self, path: str) -> None:
        pass

    # Custom logger - into the main log file / screen


    # Request trajectory, log file, text file, binary file, compressed binary file, ...
    # - this enforce simulation name prefixes and backing up instead of overwriting files in a single place
    # - the current __del__() thing in reporters is actually bad, make a single good implementation for these disk writes that ensures file closing
    # TODO all reporter custom formats should be versioned and stabilized

    # event based architecture -- using not only reporters
    # TODO try to keep S* and T* free from reporter list and handle everything through simulation

    # Topology meta-info -- constraints, vsites, coupling, extra forces
