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
        Integrator integrator = IntegratorFactory.getNumerov();
        double mass = 16.85762920 * UMA_TO_ELECTRON_MASS;
        int nOfPoints = 1000; // Adjusted so that the 14th state is weakly bound, needing an upper limit of around 50 or so
//        Grid grid = GridFactory.generateSqrtGrid(1.2d, 45., nOfPoints,5);
        Grid grid = GridFactory.generateUniformGrid(1.2d, 45., nOfPoints);
//        OutputManager.write("Grid info");
//        OutputManager.write(String.format("%10s %22s %22s %22s", "i", "r", "y", "V"));
//        for (int i = 0; i < grid.getNumberOfPoints(); i++) {
//            double v = potential.value(grid.rAtGridPoint(i));
//            OutputManager.write(String.format("%10d %22.8f %22.8f %25.6f",
//                    i, grid.rAtGridPoint(i), grid.yAtGridPoint(i), toInverseCm(v)));
//        }

        // PROBLEMS:
        // * Often with BISECTION_THEN_REGULA_FALSI we get java.lang.IllegalArgumentException: The function values at the bounds must have opposite signs.

        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        for (int nOfDesiredNodes = 0; nOfDesiredNodes <= 0; nOfDesiredNodes++) {
            QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes, RefinementStrategy.BISECTION_ONLY);

            System.out.printf("Eigenvalue: %6d %25.14f (%25.14f)\n", nOfDesiredNodes, toInverseCm(ek.energy), (toInverseCm(ek.energy) - wellDepthInverseCm));

        }


    }
}
