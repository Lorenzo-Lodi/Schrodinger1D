package schrodinger.integrator.predcorr;

import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.predcorr.PredictorCorrector8NumerovIter2;

/**
 * Test class for the PredictorCorrector8Abstract integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class PredictorCorrector8NumerovIter2Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new PredictorCorrector8NumerovIter2();
    }

}
