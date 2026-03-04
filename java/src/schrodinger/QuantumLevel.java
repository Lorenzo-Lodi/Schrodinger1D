package schrodinger;

import schrodinger.grid.Grid;
import schrodinger.potential.SchrodingerSystem;

import java.util.ArrayList;
import java.util.List;

public class QuantumLevel {
    private final SchrodingerSystem system;
    public int nodesUpper = -1; // Set to "sentinel" unphysical value
    public int nodesLower = -1; // Set to "sentinel" unphysical value
    public double energy;
    public double upperBound;
    public double lowerBound;
    public double[] psi;
    public double perturbativeCorrectionToEnergy;
    public List<ConvergenceInfo> convergenceInfo = new ArrayList<>();


    public QuantumLevel(SchrodingerSystem system) {
        this.system = system;
        this.psi = new double[system.getGrid().getNumberOfPoints()];
    }

    /**
     * Q-function Q(r) for the equation: ψ''(r) = -Q(r)ψ(r)
     */
    public double Q(double r) {
        return 2.d * system.getMass() * (energy - system.U(r));
    }

    /**
     * Transformed Q-function Q̃(y) for the mapped equation: ϕ''(y) = -Q̃(y)ϕ(y)
     * MSL Eq. (8): Q̃(y) = g²(y)·Q(r(y)) + F(y)
     */
    private double QTilde(double y) {
        return Math.pow(getGrid().g(y), 2) * Q(getGrid().r(y)) + getGrid().F(y);
    }

    public double QTildeAtGridPoint(double i) {
        double gy = getGrid().gAtGridPoint(i);
        return gy * gy * (2. * system.getMass() * (energy - system.UTildeAtGridPoint(i)));
    }

    // Note: because the wavefunctions are exponentially decreasing (or faster), the
    // trapezoidal rule is very quickly convergent (exponentially or so) until the
    // truncation error
    // at the borders takes over the global error.
    // The rectangle rule is practically the same.
    public double normalizePsi() {
        // Pre-calculate end points to avoid if-checks or method calls in loop
        int maxIndex = getGrid().getNumberOfPoints() - 1;

        // Trapezoidal rule boundaries
        double valStart = psi[0] * getGrid().gAtGridPoint(0);
        double valEnd = psi[maxIndex] * getGrid().gAtGridPoint(maxIndex);
        double sum = 0.5 * (valStart * valStart + valEnd * valEnd);

        for (int i = 1; i < maxIndex; i++) {
            double val = psi[i] * getGrid().gAtGridPoint(i);
            sum += val * val;
        }

        sum = sum * getGrid().getStepSizeYCoordinate();
        double normalizationFactor = 1.0 / Math.sqrt(sum);
        for (int i = 0; i <= maxIndex; i++) {
            psi[i] *= normalizationFactor;
        }

        return normalizationFactor;
    }

    public Grid getGrid() {
        return system.getGrid();
    }

    public double getMass() {
        return system.getMass();
    }

    public static class ConvergenceInfo {
        public String convergengeStage;
        public int iterations;
    }

}
