package schrodinger;

import org.junit.jupiter.api.Test;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.*;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.SchrodingerSystem;


public class IntegratorsTest {

    @Test
    public void propagate_test() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(1.2, 1.1);
        int nOfPoints = 10;

        Grid grid = GridFactory.generateUniformGrid(1.2, 1.6, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, 1.5, grid);
        QuantumState state = new QuantumState(system);
        state.energy = -1000.;
        state.psi[0] = 1.0;
        state.psi[1] = 1.5;

        System.out.println(new TaylorThreePoints().propagate(state.psi, 1, state, Integrator.Direction.FORWARD));
        System.out.println(new Numerov().propagate(state.psi, 1, state, Integrator.Direction.FORWARD));
        System.out.println(new Vignoli().propagate(state.psi, 1, state, Integrator.Direction.FORWARD));
        System.out.println(new ExponentiallyFitted().propagate(state.psi, 1, state, Integrator.Direction.FORWARD));

    }


}
