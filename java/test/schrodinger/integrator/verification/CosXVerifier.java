package schrodinger.integrator.verification;

import schrodinger.integrator.Integrator;

public class CosXVerifier extends ManufacturedSolutionVerifier {

    public CosXVerifier(Integrator integrator) {

        super(x -> Math.cos(x),
                x -> -Math.sin(x),
                x -> 1.,
                x -> 0.,
                x -> 0.,
                integrator);
    }

    public void propagate(Integrator.Direction direction, ConvergenceParams params) {
        super.propagate(0., 30, 40, 500, 10, direction, params);
    }

}
