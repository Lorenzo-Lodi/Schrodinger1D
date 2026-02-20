package schrodinger.integrator;

import schrodinger.QuantumState;
import schrodinger.grid.Grid;

/**
 * Test class for the PredictorCorrector8 integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class PredictorCorrector8Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new PredictorCorrector8();
    }

    @Override
    protected double getGlobalConvergenceOrder() {
        return 5.8;  // TODO investigate why it's so poor
    }
}
