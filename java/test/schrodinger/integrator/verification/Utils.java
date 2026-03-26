package schrodinger.integrator.verification;

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

}
