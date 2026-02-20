package schrodinger.integrator;

/**
 * Test class for the TaylorThreePoints integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class TaylorThreePointsTest extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new TaylorThreePoints();
    }

    @Override
    protected double getGlobalConvergenceOrder() {
        return 2.0;
    }
}
