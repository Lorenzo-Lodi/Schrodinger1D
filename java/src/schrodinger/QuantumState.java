package schrodinger;

import schrodinger.grid.Grid;
import schrodinger.potential.SchrodingerSystem;

public class QuantumState {
    public int numberOfNodes;
    public double energy;
    public double upperBound;
    public double lowerBound;
    public int numberOfBisections;
    public double[] psi;
    public double perturbativeCorrectionToEnergy;
    private final Grid grid;
    private boolean isPsiNormalized = false;

    public QuantumState(SchrodingerSystem system) {
        this.grid = system.getGrid();
        this.psi = new double[grid.getNumberOfPoints()];
    }

    // Note: because the wavefunctions are exponentially decreasing (or faster), the
    // trapezoidal rule is very quickly convergent (exponentially or so) until the
    // truncation error
    // at the borders takes over the global error.
    // The rectangle rule is practically the same.
    public double normalizePsi() {
        // Pre-calculate end points to avoid if-checks or method calls in loop
        int maxIndex = grid.getNumberOfPoints() - 1;

        // Trapezoidal rule boundaries
        double valStart = psi[0] * grid.g(0);
        double valEnd = psi[maxIndex] * grid.g(maxIndex);
        double sum = 0.5 * (valStart * valStart + valEnd * valEnd);

        for (int i = 1; i < maxIndex; i++) {
            double val = psi[i] * grid.g(i);
            sum += val * val;
        }

        sum = sum * grid.getStepSizeYCoordinate();
        double normalizationFactor = 1.0 / Math.sqrt(sum);
        for (int i = 0; i <= maxIndex; i++) {
            psi[i] *= normalizationFactor;
        }

        isPsiNormalized = true;
        return normalizationFactor;
    }

    public boolean isPsiNormalized() {
        return isPsiNormalized;
    }

}
