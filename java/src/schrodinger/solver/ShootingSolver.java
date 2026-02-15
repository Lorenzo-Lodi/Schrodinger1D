package schrodinger.solver;

import schrodinger.QuantumState;
import schrodinger.grid.Grid;
import schrodinger.integrator.Integrator;
import schrodinger.potential.SchrodingerSystem;

import java.util.Arrays;

public class ShootingSolver {
    private static final double TARGET_ABSOLUTE_ERROR = 1e-13;
    private static final int MAXIMUM_NUMBER_OF_BISECTIONS = 60; // reduces error by 2**n
    private static final double PSI_MAX = 1e30; // stop integrating forward if wave function exceeds this value
    private final Integrator integrator;
    private final SchrodingerSystem system;

    public ShootingSolver(SchrodingerSystem system, Integrator integrator) {
        this.system = system;
        this.integrator = integrator;
    }

    /**
     * Locates the energy interval [lowerBound, upperBound] containing the
     * state with 'nOfDesiredNodes' nodes.
     */
    private QuantumState findInitialEnergyBracket(int nOfDesiredNodes) {
        QuantumState level = new QuantumState(system);

        double energyScale = estimateEnergyScaleAndLowerBound(level);

        // 3. Exponential Scan to find upper bound to the energy
        double currentEnergy = level.energy;

        for (int i = 0; i < MAXIMUM_NUMBER_OF_BISECTIONS; i++) {
            // Update the energy in the system
            level.energy = currentEnergy;
            int nodes = countNodes(level);

            if (nodes > nOfDesiredNodes) {
                level.upperBound = currentEnergy;
                level.nodesUpper = nodes;
                return level;
            } else {
                level.lowerBound = currentEnergy;
                level.nodesLower = nodes;
                energyScale *= 2.0;
                currentEnergy += energyScale;
            }
        }

        throw new RuntimeException("Failed to bracket energy level.");
    }

    private Double estimateEnergyScaleAndLowerBound(QuantumState level) {
        Grid grid = system.getGrid();

        // 1. Scan the *Effective Potential* U_tilde(y) for the minimum
        int minIndex = 0;
        double uMin = Double.MAX_VALUE;

        // We scan the uniform y-grid
        int nPoints = grid.getNumberOfPoints();
        for (int i = 1; i < nPoints - 1; i++) {
            // Use the effective potential that includes mapping corrections
            double val = system.UTildeValueAt(i);
            if (val < uMin) {
                uMin = val;
                minIndex = i;
            }
        }
        level.lowerBound = uMin;
        level.nodesLower = 0; // Should be always correct

        // 2. Estimate Step Size (Energy Scale)
        Double energyScale = null;

        // Attempt A: Harmonic Curvature of U_tilde on uniform y-grid
        if (minIndex > 0 && minIndex < nPoints - 1) {

            // U_tilde values at minimum and neighbors
            double u0 = system.UTildeValueAt(minIndex);
            double uL = system.UTildeValueAt(minIndex - 1);
            double uR = system.UTildeValueAt(minIndex + 1);

            // second derivative  d^2(U_tilde)/dy^2
            double hy = grid.getStepSizeYCoordinate();
            double der2 = (uR - 2.0 * u0 + uL) / (hy * hy);

            if (der2 > 1e-15) {
                // Harmonic oscillator: omega = sqrt(K / M)
                energyScale = Math.sqrt(der2 / system.getMass());
            }
        }

        // Attempt B: Particle-in-a-Box (Fallback)
        if (energyScale == null) {
            double L = grid.getLastYValue() - grid.getFirstYValue();
            energyScale = (Math.PI * Math.PI) / (2.0 * system.getMass() * L * L);
        }

        level.energy = uMin + energyScale;

        return energyScale;

    }

    public QuantumState findEigenvalueByBisection(int nOfDesiredNodes) {

        QuantumState level = this.findInitialEnergyBracket(nOfDesiredNodes);
        refineByBisection(level, nOfDesiredNodes, TARGET_ABSOLUTE_ERROR, 0);
        level.normalizePsi();
        integrator.computePerturbativeCorrection(level);

        return level;
    }

    private void refineByBisection(QuantumState level, int nOfDesiredNodes, double maxAbsError, int minBisections) {
        int nOfBisectionAfterStrictBracketing = 0;
        for (level.numberOfBisections = 1; level.numberOfBisections <= MAXIMUM_NUMBER_OF_BISECTIONS; level.numberOfBisections++) {
            level.energy = (level.lowerBound + level.upperBound) * 0.5;
            int nodes = countNodes(level);

            if (nodes > nOfDesiredNodes) {
                level.upperBound = level.energy;
                level.nodesUpper = nodes;
            } else {
                level.lowerBound = level.energy;
                level.nodesLower = nodes;
            }

            if (level.nodesLower == nOfDesiredNodes && level.nodesUpper == nOfDesiredNodes + 1) {
                nOfBisectionAfterStrictBracketing++;
            }

            if (Math.abs(level.upperBound - level.lowerBound) < maxAbsError &&
                    nOfBisectionAfterStrictBracketing >= minBisections) {
                break;
            }

        }

    }

    public QuantumState findEigenvalueByHybridMethod(int nOfDesiredNodes) {
        QuantumState level = this.findInitialEnergyBracket(nOfDesiredNodes);
        refineByBisection(level, nOfDesiredNodes, 1e-2, 3);
        refineByBidirectionalMatching(level);
        level.normalizePsi();
        integrator.computePerturbativeCorrection(level);
        return level;
    }

    private void refineByBidirectionalMatching(QuantumState level) {

        level.energy = level.upperBound;
        double diffUpper = computeDerivativeMismatch(level);

        level.energy = level.lowerBound;
        double diffLower = computeDerivativeMismatch(level);

        // Regula falsi (false position) iteration
        double x0 = level.lowerBound;
        double f0 = diffLower;
        double x1 = level.upperBound;
        double f1 = diffUpper;

        // Ensure the bracket is valid (f0 and f1 have opposite signs)
        if (f0 * f1 > 0) {
            // Handle error: the function does not bracket a root
            throw new IllegalArgumentException("The function values at the bounds must have opposite signs.");
        }

        double x2 = x0; // initialize
        double f2;
        int maxIter = 50;   // prevent infinite loops
        double tol = 1e-12;  // tolerance on function value

        for (int iter = 0; iter < maxIter; iter++) {
            // Compute the false position point (secant line crossing zero)
            // Avoid division by zero (should not happen if f0 and f1 have opposite signs)
            x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

            // Evaluate function at x2
            level.energy = x2;
            f2 = computeDerivativeMismatch(level);

            // Check for convergence
            if (Math.abs(f2) < tol) {
                break;
            }

            // Update the bracket while keeping the root inside
            if (f0 * f2 < 0) {
                // Root lies between x0 and x2
                x1 = x2;
                f1 = f2;
            } else {
                // Root lies between x2 and x1
                x0 = x2;
                f0 = f2;
            }
        }

        // Store the final approximation
        level.energy = x2;
    }

    private double computeDerivativeMismatch(QuantumState level) {
        int matchIndex = findMatchingIndex(level.energy);
        double hy = system.getGrid().getStepSizeYCoordinate();

        // --- Shoot Forward
        Arrays.fill(level.psi, 0.0d); // Let us zero the wave function for clarity (not necessary).
        level.psi[0] = 0.0;
        level.psi[1] = 1e-16;
        for (int n = 1; n < matchIndex + 1; n++) {
            level.psi[n + 1] = integrator.propagate(level.psi, n, level, Integrator.Direction.FORWARD);
        }

        double forwardDer = (level.psi[matchIndex + 1] - level.psi[matchIndex - 1]) / (level.psi[matchIndex] * 2. * hy);
        double forwardPsiAtMatchIndexMinusOne = level.psi[matchIndex - 1];
        double forwardPsiAtMatchIndex = level.psi[matchIndex];

        // --- Shoot Backward
        int np = system.getGrid().getNumberOfPoints();
        level.psi[np - 1] = 0.0;
        level.psi[np - 2] = 1.e-16;

        for (int n = np - 2; n > matchIndex - 1; n--) {
            level.psi[n - 1] = integrator.propagate(level.psi, n, level, Integrator.Direction.BACKWARD);
        }
        double backwardDer = (level.psi[matchIndex + 1] - level.psi[matchIndex - 1]) / (level.psi[matchIndex] * 2. * hy);

        // Let us rescale the correct psi (probably unnecessary doing this at each step).
        level.psi[matchIndex - 1] = forwardPsiAtMatchIndexMinusOne;
        for (int n = 0; n < matchIndex; n++) {
            level.psi[n] = level.psi[n] / forwardPsiAtMatchIndex;
        }
        for (int n = matchIndex + 1; n < np - 1; n++) {
            level.psi[n] = level.psi[n] / level.psi[matchIndex];
        }
        level.psi[matchIndex] = 1.;

        return forwardDer - backwardDer;
    }


    private int countNodes(QuantumState level) {
        int nPoints = system.getGrid().getNumberOfPoints();

        // --- Shoot Forward
        level.psi[0] = 0.0;
        level.psi[1] = 1e-16;

        int nodes = 0;
        for (int n = 1; n < nPoints - 1; n++) {
            level.psi[n + 1] = integrator.propagate(level.psi, n, level, Integrator.Direction.FORWARD);
            if (level.psi[n] * level.psi[n + 1] < 0.0) {
                nodes++;
            }

            // If we are deep in the classically-forbidden region and the wavefunction is blowing up, we can stop early
            if (level.QTildeValueAt(n) > 2.0 * level.energy && Math.abs(level.psi[n + 1]) > PSI_MAX) {
                break;
            }
        }
        return nodes;
    }

    private int findMatchingIndex(double energy) {
        Grid grid = system.getGrid();
        // Scan from right to left until we reach the classically-allowed region
        for (int i = grid.getNumberOfPoints() - 3; i >= 2; i--) {
            if (system.UTildeValueAt(i) <= energy) {
                return i;
            }
        }

        // Fallback to midpoint if always classically forbidden
        return grid.getNumberOfPoints() / 2;
    }

}
