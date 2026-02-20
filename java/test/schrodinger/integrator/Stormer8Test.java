package schrodinger.integrator;

import schrodinger.QuantumState;
import schrodinger.grid.Grid;

/**
 * Test class for the Stormer8 integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class Stormer8Test extends AbstractIntegratorTest {

    @Override
    protected void initializeState(QuantumState state, Grid grid) {
        // Stormer8 needs three further points
        super.initializeState(state, grid);
        state.psi[4] = exactSolution(grid.getRValue(4));
        state.psi[5] = exactSolution(grid.getRValue(5));
        state.psi[6] = exactSolution(grid.getRValue(6));
    }

    @Override
    protected Integrator getIntegrator() {
        return new Stormer8();
    }
}