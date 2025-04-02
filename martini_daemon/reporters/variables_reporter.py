from .reporter import Reporter
from ..utils import backup_try
import openmm as mm  # type: ignore[import-untyped]


class VariablesReporter(Reporter):

    def on_set_xtc_path(self, xtc_name: str) -> None:
        filename = xtc_name + "_ener.csv"
        backup_try(filename)

        with open(filename, "w") as file:
            file.write("Kinetic energy,Potential energy\n")

    def on_xtc_frame(self, pos, box, xtc_name: str) -> None:
        filename = xtc_name + "_ener.csv"
        state: mm.State = self._sysstar.get_state()
        ke = state.getKineticEnergy()
        pe = state.getPotentialEnergy()

        with open(filename, "a") as file:
            file.write(f"{ke},{pe}\n")
