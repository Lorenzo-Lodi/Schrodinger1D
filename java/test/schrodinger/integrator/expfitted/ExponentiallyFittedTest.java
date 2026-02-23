package schrodinger.integrator.expfitted;

import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.Integrator;

/**
 * Test class for the ExponentiallyFitted integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class ExponentiallyFittedTest extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new ExponentiallyFitted();
    }

}
