package schrodinger.integrator.verification;

import schrodinger.integrator.Integrator;

public class MathieuSVerifier extends ManufacturedSolutionVerifier {
    // Solution of y''[x] == -(100 - 60 Cos[2 x]) y[x]
    // The solution is MathieuS[100, 30, x]
    public MathieuSVerifier(Integrator integrator) {

        super(MathieuSVerifier::MathieuS,
                MathieuSVerifier::MathieuSPrime,
                x -> 100. - 60. * Math.cos(2. * x),
                x -> 120. * Math.sin(2. * x),
                x -> 240. * Math.cos(2. * x),
                integrator);
    }

    public void propagate(Integrator.Direction direction, ConvergenceParams params) {
        super.propagate(0., 10, 40, 500, 10, direction, params);
    }

    private static double MathieuS(double x) {
        return 0.; // TODO: implement MathieuS[100, 30, x]
    }

    private static double MathieuSPrime(double x) {
        return 0.; // TODO: implement MathieuSPrime[100, 30, x]
    }


}
