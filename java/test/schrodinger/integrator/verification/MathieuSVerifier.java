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
//    private static final int N_OF_SAMPLED_GRID_POINTS = 2561;
//    private static final Path file = Path.of("test", "resources", "MathieuS_reference_values_2561.csv");
//    private static final int[] points = new int[]{
//            81, 161, 321, 641, 1281, 2561
//    };
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
        try (BufferedReader br = Files.newBufferedReader(file)) {
            br.readLine(); // skip header row

            String line;
            int i = 0;
            while ((line = br.readLine()) != null) {
                String[] parts = line.split(",");
                double x = Double.parseDouble(parts[0].trim());
                double mathieuS = Double.parseDouble(parts[1].trim());
                double mathieuSPrime = Double.parseDouble(parts[2].trim());
                sampledx[i] = x;
                sampledMathieuS[i] = mathieuS;
                sampledmathieuSPrime[i] = mathieuSPrime;
                i++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

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
        int i = computeIndex(x);
        return sampledMathieuS[i];
    }

    private static double MathieuSPrime(double x) {
        int i = computeIndex(x);
        return sampledmathieuSPrime[i];
    }

    private static int computeIndex(double x) {
        double index = (x - XMIN) / STEP;
        int i = (int) Math.round(index);
        double residualIndex = Math.abs(i - index);
        double residualX = Math.abs(x - sampledx[i]);
        if (index < 0. || index > N_OF_SAMPLED_GRID_POINTS - 1 || residualIndex > 5e-13 || residualX > 2.e-15) {
            String msg = String.format("Illegal value x = %25.16f. This leads to i = %d and an expected x = %25.16f; "
                    + "residualIndex = %20.4e, residualX = %20.4e ", x, i, sampledx[i], residualIndex, residualX);
            throw new RuntimeException(msg);
        }
        return i;

    }


}
