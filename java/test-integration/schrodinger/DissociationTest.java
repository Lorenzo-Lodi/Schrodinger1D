package schrodinger;

import org.junit.jupiter.api.Test;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.IntegratorFactory;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialLennardJones;
import schrodinger.potential.SchrodingerSystem;
import schrodinger.solver.RefinementStrategy;
import schrodinger.solver.ShootingSolver;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static schrodinger.PhysicalConstants.*;

public class DissociationTest {


    @Test
    void test001() {
        double wellDepthInverseCm = 3761;
        double wellDepthHartree = wellDepthInverseCm / HARTREE_TO_INVERSE_CM;
        double rMinAng = 1.;
        double rMinBohr = rMinAng / BOHR_TO_ANG;
        PhysicalPotential potential = new PhysicalPotentialLennardJones(rMinBohr, wellDepthHartree, 6);
        Integrator integrator = IntegratorFactory.getEfnFixedBeta();
        double mass = 16.85762920 * UMA_TO_ELECTRON_MASS;
        int nOfDesiredNodes = 14;
        int nOfPoints = 5000;
        Grid grid = GridFactory.generateUniformGrid(1.2d, 15., nOfPoints);

        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        for (nOfDesiredNodes = 0; nOfDesiredNodes <= 14; nOfDesiredNodes++) {
            QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes, RefinementStrategy.BISECTION_THEN_REGULA_FALSI);
            System.out.println("Eigenvalue: " + nOfDesiredNodes + " " + (ek.energy * HARTREE_TO_INVERSE_CM - 3761.));
        }
//        for (int i = 0; i < nOfPoints; i++) {
//            double r = grid.getRValue(i);
//            System.out.println(i + " " + r + " " + potential.value(r) * 219474.6313708 + " " + ek.psi[i]);
//        }

    }
}
