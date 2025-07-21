
from .force import Force
import openmm as mm


class PositionRestraint(Force):
    # _list contains particle id, k_x, k_y, k_z, x0, y0, z0
    # 1 member per _list entry
    _members = 1
    # external forces don't auto handle pbc, we handle it ourselves here
    _pbc = False

    def _set_force_obj(self) -> None:
        # NOTE: this only works for 90 degree angles
        self._force_obj = mm.CustomExternalForce(
            "0.5*(fx+fy+fz);"
            "fx=kx*periodicdistance(x,0,0,x0,0,0)^2;"
            "fy=ky*periodicdistance(0,y,0,0,y0,0)^2;"
            "fz=kz*periodicdistance(0,0,z,0,0,z0)^2;"
        )
        self._force_obj.addPerParticleParameter("kx")
        self._force_obj.addPerParticleParameter("ky")
        self._force_obj.addPerParticleParameter("kz")
        self._force_obj.addPerParticleParameter("x0")
        self._force_obj.addPerParticleParameter("y0")
        self._force_obj.addPerParticleParameter("z0")

    def _add_to_force_obj(self, params) -> None:
        i, kx, ky, kz, x0, y0, z0 = params
        self._force_obj.addParticle(i, (kx, ky, kz, x0, y0, z0))
