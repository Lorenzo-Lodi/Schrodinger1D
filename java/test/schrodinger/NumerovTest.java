package schrodinger;

import org.junit.jupiter.api.Test;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.*;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.SchrodingerSystem;


public class NumerovTest {

    @Test
    public void propagate_test() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(1.2, 1.1);
        int nOfPoints = 10;
        Grid grid = GridFactory.generateUniformGrid(1.2, 1.6, nOfPoints);

        SchrodingerSystem system = new SchrodingerSystem(potential, 1.5, grid);
        QuantumState state = new QuantumState(system);
        double[] psi = new double[2];
        psi[0] = 1.0;
        psi[1] = 1.5;

        System.out.println(new TaylorThreePoints().propagate(psi, 1, state, Integrator.Direction.FORWARD));
        System.out.println(new Numerov().propagate(psi, 1, state, Integrator.Direction.FORWARD));
        System.out.println(new Vignoli().propagate(psi, 1, state, Integrator.Direction.FORWARD));
        System.out.println(new ExponentiallyFitted().propagate(psi, 1, state, Integrator.Direction.FORWARD));

    }

}
