package schrodinger.integrator.verification;

import schrodinger.integrator.Integrator;

public class ExpCosXMinusXVerifier extends ManufacturedSolutionVerifier {

    public ExpCosXMinusXVerifier(Integrator integrator) {
        super(x -> Math.exp(Math.cos(x) - x / (4. * Math.PI)),
                x -> -Math.exp(Math.cos(x) - x / (4. * Math.PI)) * (Math.sin(x) + 1. / (4. * Math.PI)),
                x -> Math.cos(x) - Math.pow(1. + 4. * Math.PI * Math.sin(x), 2) / (16. * Math.PI * Math.PI),
                x -> -Math.sin(x) - (Math.cos(x) * (1. + 4. * Math.PI * Math.sin(x))) / (2. * Math.PI),
                x -> -Math.cos(x) - 2. * Math.pow(Math.cos(x), 2) + (Math.sin(x) * (1. + 4. * Math.PI * Math.sin(x))) / (2. * Math.PI),
                integrator);
    }

}
