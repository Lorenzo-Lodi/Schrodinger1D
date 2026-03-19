package schrodinger.integrator.verification;

import schrodinger.integrator.Integrator;

public class ExpCosXVerifier extends ManufacturedSolutionVerifier {

    public ExpCosXVerifier(Integrator integrator) {

        super(x -> Math.exp(Math.cos(x)),
                x -> -Math.sin(x) * Math.exp(Math.cos(x)),
                x -> Math.cos(x) - Math.pow(Math.sin(x), 2),
                x -> -Math.sin(x) - 2. * Math.cos(x) * Math.sin(x),
                x -> -Math.cos(x) - 2. * Math.pow(Math.cos(x), 2) + 2. * Math.pow(Math.sin(x), 2),
                integrator);
    }

}
