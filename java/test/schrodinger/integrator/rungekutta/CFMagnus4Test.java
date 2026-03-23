package schrodinger.integrator.rungekutta;

import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.verification.ConvergenceParams;

public class CFMagnus4Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new CFMagnus4();
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
        return new ConvergenceParams(4.00, 1.00);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_minus_x() {
        return new ConvergenceParams(1.00, 1e12); // Exact for this case
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x_minus_x() {
        return new ConvergenceParams(4.00, 1.1);
    }

    @Override
    protected ConvergenceParams convergence_params_tanh_x() {
        return new ConvergenceParams(4.00, 3.0);
    }

}
