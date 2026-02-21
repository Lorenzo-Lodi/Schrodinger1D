package schrodinger.integrator;

/**
 * Test class for the PredictorCorrector8Abstract integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class PredictorCorrector8NumerovIter1Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new PredictorCorrector8NumerovIter1();
    }

    @Override
    protected double getGlobalConvergenceOrder() {
        return 6.0;
    }
}
