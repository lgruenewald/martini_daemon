from .force import Force
import openmm as mm


class CustomDonorAcceptor(Force):
    _members = 3
    _pbc = False

    def _set_force_obj(self):
        cutoff = self._sysstar.nonbonded_force.cutoff_nm
        # TODO:
        # Bond - morse bond
        # bond - donor/acceptor parameter average
        # Angle - morse bond scale 1 to 0, 4 separate harmonic angles
        # Dihedral - morse bond scale 1 to 0,
        # 4 dihedrals of form:
        # sin(x/2 - phi/2)^2
        self._force_obj = mm.CustomHbondForce(
            "Fdist+Faa1+Faa2+Fad1+Fad2;"
            # acceptor angles
            "Faa1=a_k_angle1*(angle(a2, a1, d1)-a_theta1)^2;"
            "Faa2=a_k_angle2*(angle(a3, a1, d1)-a_theta2)^2;"
            # donor angles
            "Fad1=d_k_angle1*(angle(d2, d1, a1)-d_theta1)^2;"
            "Fad2=d_k_angle2*(angle(d3, d1, a1)-d_theta2)^2;"
            # distance force
            "Fdist=k*cdist;"
            "cdist=(distance(d1, a1)-d0)^2;"
            "d0=(d_r0+a_r0)/2;"
            "k=(d_kr0+a_kr0)/2;"
            
        )
        self._force_obj.setCutoffDistance(cutoff)
        self._force_obj.setNonbondedMethod(
            self._force_obj.CutoffPeriodic
        )
        params = [
            "r0", "kr0",
            "theta1", "k_angle1",
            "theta2", "k_angle2",
            "phi1", "k_dihedral1",
            "phi2", "k_dihedral2"
        ]

        for param in params:
            self._force_obj.addPerAcceptorParameter(
                f"a_{param}"
            )
            self._force_obj.addPerDonorParameter(
                f"d_{param}"
            )

    def _add_to_force_obj(self, params):
        i, j, k, type, *ps = params
        
        if type == "donor":
            self._force_obj.addDonor(i, j, k, ps)
        elif type == "acceptor":
            self._force_obj.addAcceptor(i, j, k, ps)

    def is_instance(self, filter):
        return filter == "custom_donor_acceptor"
