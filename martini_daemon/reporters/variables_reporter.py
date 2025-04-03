from .reporter import Reporter
from ..utils import backup_try
import openmm as mm  # type: ignore[import-untyped]
from ..sysstar import SysStar


class VariablesReporter(Reporter):

    def __init__(self):
        pass

    def on_set_xtc_path(self, xtc_name: str) -> None:
        filename = xtc_name + ".ener"
        backup_try(filename)

        with open(filename, "w") as file:
            file.write("N,Kinetic energy (kJ/mol),Potential energy (kJ/mol),"
                       "Total energy (kJ/mol),"
                       "Temperature (Kelvin),"
                       "Box X (nm),Box Y (nm),Box Z (nm),Volume (nm^3)\n")

    def on_xtc_frame(self, pos, box, xtc_name: str) -> None:
        filename = xtc_name + ".ener"
        sysstar: SysStar = self._sysstar
        n = sysstar.len_particles()
        state: mm.State = sysstar.get_state()
        ke = state.getKineticEnergy().value_in_unit(mm.unit.kilojoule_per_mole)
        pe = state.getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole)
        te = ke + pe
        T = ke / n / 0.008314 / 3 * 2
        box = state.getPeriodicBoxVectors(asNumpy=True)
        box_x = box[0][0].value_in_unit(mm.unit.nanometer)
        box_y = box[1][1].value_in_unit(mm.unit.nanometer)
        box_z = box[2][2].value_in_unit(mm.unit.nanometer)
        V = state.getPeriodicBoxVolume().value_in_unit(mm.unit.nanometer ** 3)

        with open(filename, "a") as file:
            file.write(f"{n},{ke},{pe},{te},{T},{box_x},{box_y},{box_z},{V}\n")
