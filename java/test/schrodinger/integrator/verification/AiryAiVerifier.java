package schrodinger.integrator.verification;

import schrodinger.integrator.Integrator;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public class AiryAiVerifier extends ManufacturedSolutionVerifier {
    // Solution of y''[x] == x y[x]
    // The solution is AiryAi[x]

    private static final double XMIN = -20.;
    private static final double XMAX = 5.;

    // *****************************************************************************************************************
    private static final int N_OF_SAMPLED_GRID_POINTS = 3025;
    private static final Path file = Path.of("test", "resources", "AiryAi_reference_values_3025.csv");
    private static final int[] points = new int[]{
            43, 49, 55, 57, 64, 73, 85, 109, 113, 127, 145, 169, 190, 217, 253, 337, 379, 433, 505, 757, 1009, 1513, 3025
    };
    // *****************************************************************************************************************


    private static final double STEP = (XMAX - XMIN) / (N_OF_SAMPLED_GRID_POINTS - 1);

    private static final double[] sampledx = new double[N_OF_SAMPLED_GRID_POINTS];
    private static final double[] sampledAiryAi = new double[N_OF_SAMPLED_GRID_POINTS];
    private static final double[] sampledAiryAiPrime = new double[N_OF_SAMPLED_GRID_POINTS];

    static {
        Utils.loadReferenceValuesFromFile(file, sampledx, sampledAiryAi, sampledAiryAiPrime);
    }

    public AiryAiVerifier(Integrator integrator) {

        super(AiryAiVerifier::AiryAi,
                AiryAiVerifier::AiryAiPrime,
                x -> -x,
                x -> -1.,
                x -> 0.,
                integrator);
    }

    public void propagate(Integrator.Direction direction, ConvergenceParams params) {
        super.propagate(XMIN, XMAX, points, direction, params);
    }

    private static double AiryAi(double x) {
        int i = Utils.computeIndex(x, XMIN, STEP, N_OF_SAMPLED_GRID_POINTS, sampledx);
        return sampledAiryAi[i];
    }

    private static double AiryAiPrime(double x) {
        int i = Utils.computeIndex(x, XMIN, STEP, N_OF_SAMPLED_GRID_POINTS, sampledx);
        return sampledAiryAiPrime[i];
    }

}
