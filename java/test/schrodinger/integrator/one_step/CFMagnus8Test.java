package schrodinger.integrator.one_step;

import org.junit.jupiter.api.Test;
import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.verification.ConvergenceParams;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class CFMagnus8Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new CFMagnus8();
    }

    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(0.7527195988991293, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(0.7469433084401753, val, 1e-14);
    }

    @Override
    protected Double error_threshold_harmonic_oscillator(int quantumNumber) {
        if(quantumNumber == 0) {
            return 1.5e-12;
        }
        return 5e-12;
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x() {
        return new ConvergenceParams(4.01, 0.53);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_minus_x() {
        return new ConvergenceParams(1.00, 5e11);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x_minus_x() {
        return new ConvergenceParams(4.00, 0.5);
    }

    @Override
    protected ConvergenceParams convergence_params_tanh_x() {
        return new ConvergenceParams(4.00, 0.9);
    }

}
