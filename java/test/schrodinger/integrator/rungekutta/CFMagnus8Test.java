package schrodinger.integrator.rungekutta;

import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.verification.ConvergenceParams;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class CFMagnus8Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new CFMagnus8();
    }

//    @Test
//    public void integrateForwardReferenceTest() {
//        double val = integrateOneStep(Integrator.Direction.FORWARD);
//        assertEquals(0.3996918280374979, val, 1e-14);
//    }
//
//    @Test
//    public void integrateBackwardReferenceTest() {
//        double val = integrateOneStep(Integrator.Direction.BACKWARD);
//        assertEquals(0.336941647798199, val, 1e-14);
//    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x() {
        return new ConvergenceParams(4.01, 0.00053);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_minus_x() {
        return new ConvergenceParams(1.00, 1e12);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x_minus_x() {
        return new ConvergenceParams(4.00, 0.0005);
    }

    @Override
    protected ConvergenceParams convergence_params_tanh_x() {
        return new ConvergenceParams(4.00, 0.0009);
    }

}
