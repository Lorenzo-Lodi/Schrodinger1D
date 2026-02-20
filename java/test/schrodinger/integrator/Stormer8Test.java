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
    protected Integrator getIntegrator() {
        return new Stormer8();
    }

    @Override
    protected double getGlobalConvergenceOrder() {
        // TODO verify whether the poorer-than-expected rate is correct!
        return 5.5;
    }
}
