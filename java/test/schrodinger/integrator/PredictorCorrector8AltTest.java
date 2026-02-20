package schrodinger.integrator;

import schrodinger.QuantumState;
import schrodinger.grid.Grid;

/**
 * Test class for the PredictorCorrector8Alt integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class PredictorCorrector8AltTest extends AbstractIntegratorTest {

    @Override
    protected void initializeState(QuantumState state, Grid grid) {
        // PredictorCorrector8Alt needs two further points
        super.initializeState(state, grid);
        state.psi[4] = exactSolution(grid.getRValue(4));
        state.psi[5] = exactSolution(grid.getRValue(5));
    }

    @Override
    protected Integrator getIntegrator() {
        return new PredictorCorrector8Alt();
    }
}