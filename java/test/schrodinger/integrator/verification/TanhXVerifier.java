package schrodinger.integrator.verification;

import schrodinger.integrator.Integrator;

public class TanhXVerifier extends ManufacturedSolutionVerifier {

    public TanhXVerifier(Integrator integrator) {
        super(x -> Math.tanh(x),
                x -> Math.pow(Math.cosh(x), -2),
                x -> 2. * Math.pow(Math.cosh(x), -2),
                x -> -4. * Math.pow(Math.cosh(x), -2) * Math.tanh(x),
                x -> -4. * Math.pow(Math.cosh(x), -4) + 8. * Math.pow(Math.cosh(x), -2) * Math.pow(Math.tanh(x), 2),
                integrator);
    }

    public void propagate(Integrator.Direction direction, ConvergenceParams params) {
        super.propagate(-4., 5., 10, 100, 10, direction, params);
    }

}
