package schrodinger.integrator.verification;

import schrodinger.integrator.Integrator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class MathieuSVerifier extends ManufacturedSolutionVerifier {
    // Solution of y''[x] == -(100 - 60 Cos[2 x]) y[x]
    // The solution is MathieuS[100, 30, x]

    private static final double XMIN = 0.;
    private static final double XMAX = 10.;
    private static final int N_OF_SAMPLED_GRID_POINTS = 2561;

    private static final double STEP = (XMAX - XMIN) / (N_OF_SAMPLED_GRID_POINTS - 1);

    private static final double[] sampledx = new double[N_OF_SAMPLED_GRID_POINTS];
    private static final double[] sampledMathieuS = new double[N_OF_SAMPLED_GRID_POINTS];
    private static final double[] sampledmathieuSPrime = new double[N_OF_SAMPLED_GRID_POINTS];

    static {
        Path file = Path.of("test", "resources", "MathieuS_reference_values.csv");

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
        int[] points = new int[]{
                81, 161, 321, 641, 1281, 2561
        };
        super.propagate(XMIN, XMAX, points, direction, params);
    }

    private static double MathieuS(double x) {
        double index = (x - XMIN) / STEP;
        int i = (int) Math.round(index);
        double residual = Math.abs(i - index);
        if (index < 0. || index > N_OF_SAMPLED_GRID_POINTS - 1 || residual > 1e-16 || Math.abs(x - sampledx[i])
                > 1.e-16) {
            String msg = String.format("Illegal value for x = %20.6f, i = %d, expected x = %20.6f", x, i, sampledx[i]);
            throw new RuntimeException(msg);
        }
//        System.out.println("called for x = " + x + "  index=" + index + " residual = " + residual);
        return 0.; // TODO: implement MathieuS[100, 30, x]
    }

    private static double MathieuSPrime(double x) {
        return 0.; // TODO: implement MathieuSPrime[100, 30, x]
    }


}
