package schrodinger.integrator.rungekutta;

import org.junit.jupiter.api.Test;
import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.verification.ConvergenceParams;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class CFMagnus6e4Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new CFMagnus6e4();
    }

    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(0.7494191043668578, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(0.7516639646717973, val, 1e-14);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x() {
        return new ConvergenceParams(6.00, 0.6);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_minus_x() {
        return new ConvergenceParams(1.00, 5e11);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x_minus_x() {
        return new ConvergenceParams(6.00, 0.6);
    }

    @Override
    protected ConvergenceParams convergence_params_tanh_x() {
        return new ConvergenceParams(6.00, 1.6);
    }

}
