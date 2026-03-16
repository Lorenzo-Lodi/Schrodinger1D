package schrodinger.integrator.verification;

public class LinearFit {

    /**
     * Fits data to y = a + b*x using ordinary least squares.
     *
     * @param x independent variable
     * @param y dependent variable
     * @return double[]{a, b, rSquared}
     */
    public static double[] fit(double[] x, double[] y) {
        if (x == null || y == null || x.length != y.length || x.length < 2) {
            throw new IllegalArgumentException("Need at least 2 matching data points");
        }

        int n = x.length;
        double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0, sumYY = 0;

        // First pass: compute all necessary sums
        for (int i = 0; i < n; i++) {
            double xi = x[i], yi = y[i];
            sumX += xi;
            sumY += yi;
            sumXY += xi * yi;
            sumXX += xi * xi;
            sumYY += yi * yi;  // Required for R² calculation
        }

        double denom = n * sumXX - sumX * sumX;
        if (denom == 0) throw new ArithmeticException("Singular matrix (all x values identical)");

        // Calculate slope (b) and intercept (a)
        double b = (n * sumXY - sumX * sumY) / denom;
        double a = (sumY - b * sumX) / n;

        // Calculate R² = 1 - (SS_res / SS_tot)
        double meanY = sumY / n;
        double ssTot = sumYY - (sumY * sumY) / n;  // Total Sum of Squares

        // If all y values are identical, R² is 1.0 by convention if fit is perfect
        if (ssTot == 0.0) return new double[]{a, b, 1.0};

        double ssRes = 0;  // Residual Sum of Squares
        for (int i = 0; i < n; i++) {
            double yPred = a + b * x[i];
            ssRes += (y[i] - yPred) * (y[i] - yPred);
        }

        double rSquared = 1.0 - (ssRes / ssTot);
        return new double[]{a, b, rSquared};
    }
}