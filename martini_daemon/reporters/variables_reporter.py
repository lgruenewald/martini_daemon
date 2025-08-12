from .reporter import Reporter
import openmm as mm
from ..sysstar import SysStar


class VariablesReporter(Reporter):

    def __init__(self):
        pass

    def on_set_xtc_path(self, xtc_name: str) -> None:
        self._open(xtc_name + ".ener")
        self._print(
            "N,Kinetic energy (kJ/mol),Potential energy (kJ/mol),"
            "Total energy (kJ/mol),"
            "Temperature (Kelvin),"
            "Box X (nm),Box Y (nm),Box Z (nm),Volume (nm^3)"
        )

    def on_xtc_frame(self, i, pos, box, xtc_name: str) -> None:
        sysstar: SysStar = self._sysstar
        n = sysstar.len_atoms()
        state: mm.State = sysstar.get_state()
        ke = state.getKineticEnergy().value_in_unit(mm.unit.kilojoule_per_mole)
        pe = state.getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole)
        te = ke + pe
        n_constraints = len(sysstar.constraint)
        n_vsite = len(sysstar.vsites)
        degrees_of_freedom = n * 3 - n_constraints - n_vsite * 3
        if sysstar.remove_com:
            degrees_of_freedom -= 3
        assert degrees_of_freedom > 0
        T = ke / degrees_of_freedom / 0.008314 * 2
        box = state.getPeriodicBoxVectors(asNumpy=True)
        box_x = box[0][0].value_in_unit(mm.unit.nanometer)
        box_y = box[1][1].value_in_unit(mm.unit.nanometer)
        box_z = box[2][2].value_in_unit(mm.unit.nanometer)
        V = state.getPeriodicBoxVolume().value_in_unit(mm.unit.nanometer ** 3)

        self._print(f"{n},{ke},{pe},{te},{T},{box_x},{box_y},{box_z},{V}")
