package schrodinger.solver;

import schrodinger.OutputManager;
import schrodinger.QuantumLevel;
import schrodinger.pt_correction.PTCorrection;
import schrodinger.grid.Grid;
import schrodinger.integrator.Integrator;
import schrodinger.potential.SchrodingerSystem;

import java.util.Arrays;
import java.util.function.DoubleUnaryOperator;

public class ShootingSolver {
    private static final double TARGET_ABSOLUTE_ERROR = 1e-13;
    private static final int MAXIMUM_NUMBER_OF_BISECTIONS = 60; // reduces error by 2**n
    private static final double PSI_MAX = 1e30; // stop integrating forward if wave function exceeds this value
    private final Integrator integrator;
    private final PTCorrection correction;
    private final SchrodingerSystem system;
    private RefinementStrategy strategy;

    public ShootingSolver(SchrodingerSystem system, Integrator integrator, PTCorrection correction) {
        this.system = system;
        this.integrator = integrator;
        this.correction = correction;
    }

    // Backward-compatible constructor (correction set to null)
    public ShootingSolver(SchrodingerSystem system, Integrator integrator) {
        this(system, integrator, null);
    }

    /**
     * Locates the energy interval [lowerBound, upperBound] containing the
     * state with 'nOfDesiredNodes' nodes.
     */
    private QuantumLevel findInitialEnergyBracket(int nOfDesiredNodes) {
        QuantumLevel level = new QuantumLevel(system);

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
                level.energy = 0.5 * (level.lowerBound + level.upperBound);
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

    private Double estimateEnergyScaleAndLowerBound(QuantumLevel level) {
        Grid grid = system.getGrid();

        // 1. Scan the *Effective Potential* U_tilde(y) for the minimum
        int minIndex = 0;
        double uMin = Double.MAX_VALUE;

        // We scan the uniform y-grid
        int nPoints = grid.getNumberOfPoints();
        for (int i = 1; i < nPoints - 1; i++) {
            // Use the effective potential that includes mapping corrections
            double val = system.UTildeAtGridPoint(i);
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
            double u0 = system.UTildeAtGridPoint(minIndex);
            double uL = system.UTildeAtGridPoint(minIndex - 1);
            double uR = system.UTildeAtGridPoint(minIndex + 1);

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

    private QuantumLevel findEigenvalueByBisection(int nOfDesiredNodes) {

        QuantumLevel level = this.findInitialEnergyBracket(nOfDesiredNodes);
        refineByBisection(level, nOfDesiredNodes, TARGET_ABSOLUTE_ERROR, 0);
        level.normalizePsi();
        if (correction != null) {
            correction.compute(level);
        }

        return level;
    }

    private void refineByBisection(QuantumLevel level, int nOfDesiredNodes, double maxAbsError, int minBisections) {
        QuantumLevel.ConvergenceInfo info1 = new QuantumLevel.ConvergenceInfo();
        info1.convergengeStage = "Bisection (to strick bracketing)";
        level.convergenceInfo.add(info1);
        QuantumLevel.ConvergenceInfo info2 = new QuantumLevel.ConvergenceInfo();
        info2.convergengeStage = "Bisection (after strict bracketing)";
        info2.iterations = 0;
        level.convergenceInfo.add(info2);

        for (int i = 1; i <= MAXIMUM_NUMBER_OF_BISECTIONS; i++) {
            info1.iterations = i;
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
                info2.iterations++;
            }

            if (Math.abs(level.upperBound - level.lowerBound) < maxAbsError &&
                    info2.iterations >= minBisections) {
                break;
            }

        }
        level.energy = (level.lowerBound + level.upperBound) * 0.5;
    }

    public QuantumLevel findEigenvalue(int nOfDesiredNodes) {
        return findEigenvalue(nOfDesiredNodes, RefinementStrategy.BISECTION_THEN_REGULA_FALSI);
    }

    public QuantumLevel findEigenvalue(int nOfDesiredNodes, RefinementStrategy refinementStrategy) {
        OutputManager.write(String.format("Finding eigenvalue with %d nodes, strategy %s", nOfDesiredNodes, refinementStrategy.toString()));
        this.strategy = refinementStrategy;
        switch (strategy) {
            case BISECTION_ONLY:
                return findEigenvalueByBisection(nOfDesiredNodes);
            case BISECTION_THEN_SECANT:
            case BISECTION_THEN_REGULA_FALSI:
                return findEigenvalueByHybridMethod(nOfDesiredNodes);
        }
        return null;
    }

    private QuantumLevel findEigenvalueByHybridMethod(int nOfDesiredNodes) {
        QuantumLevel level = this.findInitialEnergyBracket(nOfDesiredNodes);
        OutputManager.write(String.format("Initial LOWER energy is: %23.14f ", level.lowerBound));
        OutputManager.write(String.format("Initial GUESS energy is: %23.14f ", level.energy));
        OutputManager.write(String.format("Initial UPPER energy is: %23.14f ", level.upperBound));
        refineByBisection(level, nOfDesiredNodes, 1e-2, 3);
        OutputManager.write("Initial bisection steps finished. The new brackets are:");
        OutputManager.write(String.format("Initial LOWER energy is: %23.14f ", level.lowerBound));
        OutputManager.write(String.format("Initial GUESS energy is: %23.14f ", level.energy));
        OutputManager.write(String.format("Initial UPPER energy is: %23.14f ", level.upperBound));
        refineByBidirectionalMatching(level);
        level.normalizePsi();
        if (correction != null) {
            correction.compute(level);
        }
        return level;
    }

    private void refineByBidirectionalMatching(QuantumLevel level) {
        QuantumLevel.ConvergenceInfo info = new QuantumLevel.ConvergenceInfo();
        info.convergengeStage = "Refinement by regula falsi or secant";
        info.iterations = 0;
        level.convergenceInfo.add(info);

        // It seems preferable to compute the matching index once and for all
        int matchIndex = findMatchingIndex(level.energy);

        level.energy = level.upperBound;
        double diffUpper = computeDerivativeMismatch(level, matchIndex);
        info.iterations++;

        level.energy = level.lowerBound;
        double diffLower = computeDerivativeMismatch(level, matchIndex);
        info.iterations++;

        OutputManager.write(String.format("Energy refinement stage. iterations = %d", info.iterations));
        OutputManager.write(String.format("Derivative mismatch for UPPER energy: %25.10f", diffUpper));
        OutputManager.write(String.format("Derivative mismatch for LOWER energy: %25.10f", diffLower));

        // Regula falsi (false position) iteration
        double x0 = level.lowerBound;
        double f0 = diffLower;
        double x1 = level.upperBound;
        double f1 = diffUpper;

        boolean isFalsePosition = (strategy == RefinementStrategy.BISECTION_THEN_REGULA_FALSI);

        // Ensure the bracket is valid (f0 and f1 have opposite signs)
        if (isFalsePosition && f0 * f1 > 0) {
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
            f2 = computeDerivativeMismatch(level, matchIndex);
            info.iterations++;

            // Check for convergence
            if (Math.abs(f2) < tol) {
                break;
            }

            if (isFalsePosition) {
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
            } else {
                // Secant update: always shift forward, discard oldest point
                x0 = x1;
                f0 = f1;
                x1 = x2;
                f1 = f2;
            }
        }

        // Store the final approximation
        level.energy = x2;
    }

    private double computeDerivativeMismatch(QuantumLevel level, int matchIndex) {
        double hy = system.getGrid().getStepSizeYCoordinate();
        DoubleUnaryOperator qTildeFunction = level::QTildeAtGridPoint;

        // --- Shoot Forward
        Arrays.fill(level.psi, 0.0d); // Let us zero the wave function for clarity (not necessary).
        level.psi[0] = 0.0;
        level.psi[1] = 1e-16;
        for (int n = 1; n < matchIndex + 1; n++) {
            level.psi[n + 1] = integrator.propagate(level.psi, n, hy, qTildeFunction, Integrator.Direction.FORWARD);

            // Check for potential overflow
            if (n % 16 == 0 && Math.abs(level.psi[n + 1]) > 1e100) {
                for (int i = 0; i <= n + 1; i++) { // Rescale computed points
                    level.psi[i] /= 1e200;
                }
            }

        }

        double forwardDer = (level.psi[matchIndex + 1] - level.psi[matchIndex - 1]) / (level.psi[matchIndex] * 2. * hy);
        double forwardPsiAtMatchIndexMinusOne = level.psi[matchIndex - 1];
        double forwardPsiAtMatchIndex = level.psi[matchIndex];

        // --- Shoot Backward
        int np = system.getGrid().getNumberOfPoints();
        level.psi[np - 1] = 0.0;
        level.psi[np - 2] = 1.e-16;

        for (int n = np - 2; n > matchIndex - 1; n--) {
            level.psi[n - 1] = integrator.propagate(level.psi, n, hy, qTildeFunction, Integrator.Direction.BACKWARD);

            // Check for potential overflow
            if (n % 16 == 0 && Math.abs(level.psi[n - 1]) > 1e100) {
                for (int i = np - 1; i >= n - 1; i--) { // Rescale computed points
                    level.psi[i] /= 1e200;
                }
            }

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


    private int countNodes(QuantumLevel level) {
        int nPoints = system.getGrid().getNumberOfPoints();
        double hy = system.getGrid().getStepSizeYCoordinate();
        DoubleUnaryOperator qTildeFunction = level::QTildeAtGridPoint;

        // --- Shoot Forward
        level.psi[0] = 0.0;
        level.psi[1] = 1e-16;

        int nodes = 0;
        for (int n = 1; n < nPoints - 1; n++) {
            level.psi[n + 1] = integrator.propagate(level.psi, n, hy, qTildeFunction, Integrator.Direction.FORWARD);
            if (level.psi[n] * level.psi[n + 1] < 0.0) {
                nodes++;
            }

            // If we are deep in the classically-forbidden region and the wavefunction is blowing up, we can stop early
            if (level.QTildeAtGridPoint(n) > 2.0 * level.energy && Math.abs(level.psi[n + 1]) > PSI_MAX) {
                break;
            }
        }
        return nodes;
    }

    private int findMatchingIndex(double energy) {
        Grid grid = system.getGrid();
        // Scan from right to left until we reach the classically-allowed region
        for (int i = grid.getNumberOfPoints() - 3; i >= 2; i--) {
            if (system.UTildeAtGridPoint(i) <= energy) {
                return i;
            }
        }

        // Fallback to midpoint if always classically forbidden
        return grid.getNumberOfPoints() / 2;
    }

}
