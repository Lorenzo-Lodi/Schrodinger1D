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

public class DissociationTest {


    @Test
    void test001() {
        PhysicalPotential potential = new PhysicalPotentialLennardJones(1.8897259886, 3761. / 219474.6313708, 6);
        Integrator integrator = IntegratorFactory.getVignoli();
        double mass = 16.85762920 * 1822.89;
        int nOfDesiredNodes = 14;
        int nOfPoints = 5000;
        Grid grid = GridFactory.generateUniformGrid(1.2d, 60., nOfPoints);

        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        QuantumState ek = finder.findEigenvalue(nOfDesiredNodes, RefinementStrategy.BISECTION_THEN_REGULA_FALSI);
        System.out.println("Eigenvalue: " + (ek.energy * 219474.6313708 - 3761.));

//        for (int i = 0; i < nOfPoints; i++) {
//            double r = grid.getRValue(i);
//            System.out.println(i + " " + r + " " + potential.value(r) * 219474.6313708 + " " + ek.psi[i]);
//        }

    }
}
