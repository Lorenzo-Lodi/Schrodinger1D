package schrodinger.integrator;

/**
 * Test class for the Numerov integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class NumerovTest extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new Numerov();
    }

}
