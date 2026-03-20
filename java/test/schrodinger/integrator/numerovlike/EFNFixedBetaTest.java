package schrodinger.integrator.numerovlike;

import org.junit.jupiter.api.Test;
import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.verification.ConvergenceParams;
import schrodinger.integrator.Integrator;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test class for "Exponentially-Fitted Numerov with fixed beta" integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class EFNFixedBetaTest extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new EFNFixedBeta();
    }

    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(2.081570296966459, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(1.6758023010055854, val, 1e-14);
    }

    @Test
    public void integrateForwardNegativeQTest() {
        double val = integrateOneStepNegativeQ(Integrator.Direction.FORWARD);
        assertEquals(-22.088386420352787, val, 1e-14);
    }

    @Test
    public void integrateBackwardNegativeQTest() {
        double val = integrateOneStepNegativeQ(Integrator.Direction.BACKWARD);
        assertEquals(-8.88050151774844, val, 1e-14);
    }

    @Test
    public void integrateForwardSmallZTest() {
        double val = integrateOneStepSmallZ(Integrator.Direction.FORWARD);
        assertEquals(-0.7668908395120492, val, 1e-14);
    }

    @Test
    public void integrateBackwardSmallZTest() {
        double val = integrateOneStepSmallZ(Integrator.Direction.BACKWARD);
        assertEquals(-1.6824522708379095, val, 1e-14);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x() {
        return new ConvergenceParams(3.54, 0.9);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_minus_x() {
        // It is actually exact for this function (to roundoff), hence the huge characteristicStepSize
        return new ConvergenceParams(1.00, 30000000000.);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x_minus_x() {
        return new ConvergenceParams(4.00, 0.65);
    }

    @Override
    protected ConvergenceParams convergence_params_tanh_x() {
        return new ConvergenceParams(4.00, 1.7);
    }

}
