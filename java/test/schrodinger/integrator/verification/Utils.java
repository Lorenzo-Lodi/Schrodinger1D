package schrodinger.integrator.verification;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public class Utils {

    private Utils() {
    }

    static int computeIndex(double x, double XMIN, double STEP, int N_OF_SAMPLED_GRID_POINTS, double[] sampledx) {
        double index = (x - XMIN) / STEP;
        int i = (int) Math.round(index);
        double residualIndex = Math.abs(i - index);
        double residualX = Math.abs(x - sampledx[i]);
        if (index < 0. || index > N_OF_SAMPLED_GRID_POINTS - 1 || residualIndex > 5e-13 || residualX > 4.e-15) {
            String msg = String.format("Illegal value x = %25.16f. This leads to i = %d and an expected x = %25.16f; "
                    + "residualIndex = %20.4e, residualX = %20.4e ", x, i, sampledx[i], residualIndex, residualX);
            throw new RuntimeException(msg);
        }
        return i;

    }

    static void loadReferenceValuesFromFile(Path file, double[] xArray, double[] fArray, double[] fPrimeArray) {
        try (BufferedReader br = Files.newBufferedReader(file)) {
            br.readLine(); // skip header row

            String line;
            int i = 0;
            while ((line = br.readLine()) != null) {
                String[] parts = line.split(",");
                double x = Double.parseDouble(parts[0].trim());
                double mathieuS = Double.parseDouble(parts[1].trim());
                double mathieuSPrime = Double.parseDouble(parts[2].trim());
                xArray[i] = x;
                fArray[i] = mathieuS;
                fPrimeArray[i] = mathieuSPrime;
                i++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
