package schrodinger.integrator.verification;

import schrodinger.integrator.Integrator;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public class MathieuSVerifier extends ManufacturedSolutionVerifier {
    // Solution of y''[x] == -(100 - 60 Cos[2 x]) y[x]
    // The solution is MathieuS[100, 30, x]

    private static final double XMIN = 0.;
    private static final double XMAX = 10.;

    // *****************************************************************************************************************
    private static final int N_OF_SAMPLED_GRID_POINTS = 3025;
    private static final Path file = Path.of("test", "resources", "MathieuS_reference_values_3025.csv");
    private static final int[] points = new int[]{
            43, 49, 55, 57, 64, 73, 85, 109, 113, 127, 145, 169, 190, 217, 253, 337, 379, 433, 505, 757, 1009, 1513, 3025
    };
    // *****************************************************************************************************************


    private static final double STEP = (XMAX - XMIN) / (N_OF_SAMPLED_GRID_POINTS - 1);

    private static final double[] sampledx = new double[N_OF_SAMPLED_GRID_POINTS];
    private static final double[] sampledMathieuS = new double[N_OF_SAMPLED_GRID_POINTS];
    private static final double[] sampledmathieuSPrime = new double[N_OF_SAMPLED_GRID_POINTS];

    static {
        Utils.loadReferenceValuesFromFile(file, sampledx, sampledMathieuS, sampledmathieuSPrime);
    }


    public MathieuSVerifier(Integrator integrator) {

        super(MathieuSVerifier::MathieuS,
                MathieuSVerifier::MathieuSPrime,
                x -> 100. - 60. * Math.cos(2. * x),
                x -> 120. * Math.sin(2. * x),
                x -> 240. * Math.cos(2. * x),
                integrator);
    }

    public void propagate(Integrator.Direction direction, ConvergenceParams params) {
        super.propagate(XMIN, XMAX, points, direction, params);
    }

    private static double MathieuS(double x) {
        int i = Utils.computeIndex(x, XMIN, STEP, N_OF_SAMPLED_GRID_POINTS, sampledx);
        return sampledMathieuS[i];
    }

    private static double MathieuSPrime(double x) {
        int i = Utils.computeIndex(x, XMIN, STEP, N_OF_SAMPLED_GRID_POINTS, sampledx);
        return sampledmathieuSPrime[i];
    }

}
