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
        OutputManager.initCommonOutputFile(this.getClass().getSimpleName() + ".log");
        double wellDepthInverseCm = 4050;
        double wellDepthHartree = wellDepthInverseCm / HARTREE_TO_INVERSE_CM;
        double rMinAng = 1.;
        double rMinBohr = rMinAng / BOHR_TO_ANG;
        PhysicalPotential potential = new PhysicalPotentialLennardJones(rMinBohr, wellDepthHartree, 6);
        Integrator integrator = IntegratorFactory.getExponentiallyFitted();
        double mass = 16.85762920 * UMA_TO_ELECTRON_MASS;
        int nOfPoints = 4000; // Adjusted so that the 14th state is weakly bound, needing an upper limit of around 50 or so
        Grid grid = GridFactory.generateUniformGrid(1.2d, 45., nOfPoints);

        // PROBLEMS:
        // 1. Often the answer is a NaN
        // 2. Often with BISECTION_THEN_REGULA_FALSI we get java.lang.IllegalArgumentException: The function values at the bounds must have opposite signs.

        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
//        for (int nOfDesiredNodes = 0; nOfDesiredNodes <= 14; nOfDesiredNodes++) {
        for (int nOfDesiredNodes = 0; nOfDesiredNodes <= 0; nOfDesiredNodes++) {
            QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes, RefinementStrategy.BISECTION_THEN_SECANT);
            System.out.println("Eigenvalue: " + nOfDesiredNodes + " " + (ek.energy * HARTREE_TO_INVERSE_CM - wellDepthInverseCm));

//            if (nOfDesiredNodes == 14) {
//                for (int i = 0; i < nOfPoints; i++) {
//                    double r = grid.getRValue(i);
////                    System.out.println(i + " " + r + " " + potential.value(r) * 219474.6313708 + " " + ek.psi[i]);
//                    System.out.println( ek.psi[i]);
//                }
//            }
        }


    }
}
