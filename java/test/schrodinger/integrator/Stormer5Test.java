package schrodinger.integrator;

/**
 * Test class for the Stormer5 integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class Stormer5Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new Stormer5();
    }

    @Override
    protected double getGlobalConvergenceOrder() {
        return 5.0;
    }
}
