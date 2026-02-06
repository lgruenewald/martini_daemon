from .reporter import Reporter
import openmm as mm
from ..sysstar import SysStar


class ReactionEnergyReporter(Reporter):

    def __init__(self, write_coords=False, force_groups=False):
        """
        Designed to work together with DaemonIntegrators (see martini_daemon.components subpackage)

        may not work without!

        write_coords, if set to True, it will print .gro files if daemon integrators are used
        force_groups, if set to True will print a breakdown of forces
        """
        self.write_coords = write_coords
        self.force_groups = force_groups

    def on_set_xtc_path(self, xtc_name: str) -> None:
        self._open(xtc_name + ".rxener")
        self._print(
            "Entry,N,Kinetic energy (kJ/mol),Potential energy (kJ/mol),"
            "Total energy (kJ/mol),"
            "Temperature (Kelvin),"
            "Box X (nm),Box Y (nm),Box Z (nm),Volume (nm^3)"
        )

    def write_energies(self, title):
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

        self._print(f"{title},{n},{ke},{pe},{te},{T},{box_x},{box_y},{box_z},{V}")

    def post_reinitialize(self, i):
        self.i = i
        self.write_energies(f"After reaction frame {i}")
        if self.write_coords:
            self._sysstar.write_gro(f"premin{i}.gro")
        self._write_force_groups(f"premin{self.i}")

    def post_sc_enable(self, i):
        assert i == self.i
        self.write_energies(f"After soft-core enable {i}")

    def post_di_minimize(self):
        self.write_energies(f"After post-reaction minimization {self.i}")
        if self.write_coords:
            self._sysstar.write_gro(f"postmin{self.i}.gro")
        self._write_force_groups(f"postmin{self.i}")

    def post_di_equilibrate(self):
        self.write_energies(f"After post-reaction equilibration {self.i}")
        if self.write_coords:
            self._sysstar.write_gro(f"posteq{self.i}.gro")
        self._write_force_groups(f"posteq{self.i}")

    def post_sc_disable(self, i):
        assert i == self.i
        self.write_energies(f"After soft-core disable {i}")

    def _write_force_groups(self, title):
        if self.force_groups:
            msg = self._sysstar.get_energies_and_forces_by_group()
            with open(title + ".txt", "w") as f:
                f.write(msg)
