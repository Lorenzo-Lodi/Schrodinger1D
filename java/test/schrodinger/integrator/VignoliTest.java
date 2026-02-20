package schrodinger.integrator;

/**
 * Test class for Vignoli integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class VignoliTest extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new Vignoli();
    }

    @Override
    protected double getGlobalConvergenceOrder() {
        return 4.0;
    }
}
