package schrodinger.integrator.verification;

import schrodinger.integrator.Integrator;

public class ExpMinusXVerifier extends ManufacturedSolutionVerifier {

    public ExpMinusXVerifier(Integrator integrator) {
        super(x -> Math.exp(-1.5 * x),
                x -> -1.5 * Math.exp(-1.5 * x),
                x -> -2.25,
                x -> 0.,
                x -> 0.,
                integrator);
    }

    public void propagate(Integrator.Direction direction, ConvergenceParams params) {
        super.propagate(0., 2., 10, 100, 10, direction, params);
    }

}
