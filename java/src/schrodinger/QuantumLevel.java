package schrodinger;

import schrodinger.grid.Grid;

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
    public double[] currentPsiPrime; // Only by one-step methods such as RKN and CFMagnus
    public double perturbativeCorrectionToEnergy;
    public List<ConvergenceInfo> convergenceInfo = new ArrayList<>();
    private boolean isCapPotential;

    public QuantumLevel(SchrodingerSystem system) {
        this.system = system;
        this.psi = new double[system.getGrid().getNumberOfPoints()];
        this.currentPsiPrime = new double[1];
    }

    /**
     * Q-function Q(r) for the equation: ψ''(r) = -Q(r)ψ(r)
     */
    private double Q(double r) {
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
        double qTilde = gy * gy * 2. * system.getMass() * (energy - system.UTildeAtGridPoint(i));
        if (isCapPotential) {
            qTilde = Math.max(qTilde, system.getQMin());
        }
        return qTilde;
    }

    public double QTildePrimeAtGridPoint(double i) {
        double eps = 1e-6;  // TODO
        double y = getGrid().yAtGridPoint(i);
        double der1 = (QTilde(y + eps) - QTilde(y - eps)) / (2. * eps);
        double der2 = (QTilde(y + 2. * eps) - QTilde(y - 2. * eps)) / (4. * eps);
        return (4. / 3.) * der1 - (1. / 3.) * der2;
    }

    public double QTildeDoublePrimeAtGridPoint(double i) {
        double eps = 1e-4;  // TODO
        double y = getGrid().yAtGridPoint(i);
        double der1 = (QTilde(y + eps) + QTilde(y - eps) - 2. * QTilde(y)) / (eps * eps);
        double der2 = (QTilde(y + 2. * eps) + QTilde(y - 2. * eps) - 2. * QTilde(y)) / (4. * eps * eps);
        return (4. / 3.) * der1 - (1. / 3.) * der2;
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

    public int countNodes() {
        int nodes = 0;
        for (int n = 0; n < this.psi.length - 1; n++) {
            // if psi is very small, we risk underflowing when we perfom the multiplication.
            // Use this trick to get the sign bit; works also for 0.0 vs -0.0
            if ((Double.doubleToRawLongBits(psi[n]) ^ Double.doubleToRawLongBits(psi[n + 1])) < 0L) {
                nodes++;
            }

        }

        return nodes;
    }

    public Grid getGrid() {
        return system.getGrid();
    }

    public double getMass() {
        return system.getMass();
    }

    public void setCapPotential(boolean capPotential) {
        isCapPotential = capPotential;
    }

    public static class ConvergenceInfo {
        public String convergengeStage;
        public int iterations;
    }

    public int countTotalScans() {
        int scans = 0;
        for (ConvergenceInfo c : this.convergenceInfo) {
            scans += c.iterations;
        }
        return scans;
    }

    public double maximumStepSize() {
        double Qmax = QTildeAtGridPoint(system.UTildeMinimumGridIndex);
        double smallestLambda = 2. * Math.PI / Math.sqrt(Qmax);
        return smallestLambda / 8.; // Divide the minimum lambda by a factor 2-10;
    }

    public void verifyStepSize() {
        double stepSize = this.getGrid().getStepSizeYCoordinate();
        double maxStepSize = this.maximumStepSize();
        if (maxStepSize < stepSize) {
            OutputManager.write(String.format("WARNING: current step size is %20.6f but the computed minimum step size is %20.6f",
                    stepSize, maxStepSize));
            OutputManager.write(String.format("Ratio currentStepSize/maxStepSize (should be less than 1.000): %20.4f",
                    stepSize / maxStepSize));
            // For now just warn, throw an exception in the future...
//            String msg = String.format("Current step size %20.6f is too large! Maximum step size = %20.6f", stepSize, maxStepSize);
//                    throw new RuntimeException(msg);
        }
    }

}
