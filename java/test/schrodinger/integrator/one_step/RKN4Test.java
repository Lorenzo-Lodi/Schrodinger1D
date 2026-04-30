package schrodinger.integrator.one_step;

import org.junit.jupiter.api.Test;
import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.verification.ConvergenceParams;
import schrodinger.integrator.Integrator;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class RKN4Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new RKN4();
    }

    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(0.3996918280374979, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(0.336941647798199, val, 1e-14);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x() {
        return new ConvergenceParams(4.01, 0.53);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_minus_x() {
        return new ConvergenceParams(4.00, 1.3);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x_minus_x() {
        return new ConvergenceParams(4.00, 0.5);
    }

    @Override
    protected ConvergenceParams convergence_params_tanh_x() {
        return new ConvergenceParams(4.00, 1.9);
    }

}
