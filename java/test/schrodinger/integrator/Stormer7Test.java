package schrodinger.integrator;

import schrodinger.QuantumState;
import schrodinger.grid.Grid;

/**
 * Test class for the Stormer7 integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class Stormer7Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new Stormer7();
    }

    @Override
    protected double getGlobalConvergenceOrder() {
        return 4.9;  // TODO investigate why it's so poor
    }
}
