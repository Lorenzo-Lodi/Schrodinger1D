package schrodinger.solver;

import schrodinger.OutputManager;
import schrodinger.QuantumLevel;
import schrodinger.integrator.IntegratorFactory;
import schrodinger.integrator.pt_correction.PTCorrector;
import schrodinger.grid.Grid;
import schrodinger.integrator.Integrator;
import schrodinger.potential.SchrodingerSystem;

import java.util.Arrays;
import java.util.function.DoubleUnaryOperator;

import static schrodinger.PhysicalConstants.toInverseCm;

public class ShootingSolver {
    private static final double TARGET_ABSOLUTE_ERROR = 1e-13;
    private static final int MAXIMUM_NUMBER_OF_BISECTIONS = 60; // reduces error by 2**n
    private static final double PSI_MAX = 1e140; // rescale if psi exceeds this value
    private final Integrator integrator;
    private final SchrodingerSystem system;
    private final PTCorrector corrector;
    private RefinementStrategy strategy;

    private final Integrator integratorBest = IntegratorFactory.getBestOneStepIntegrator();

    public ShootingSolver(SchrodingerSystem system, Integrator integrator) {
        this.system = system;
        this.integrator = integrator;
        this.corrector = integrator.getPertubativeCorrector();
    }

    /**
     * Locates the energy interval [lowerBound, upperBound] containing the
     * state with 'nOfDesiredNodes' nodes.
     */
    private QuantumLevel findInitialEnergyBracket(int nOfDesiredNodes) {
        OutputManager.writeBlankLine();
        OutputManager.write(String.format("Trying to find initial energy bracketing for state with n = %d", nOfDesiredNodes));
        QuantumLevel level = new QuantumLevel(system);

        if (integrator.needsPotentialCapping()) {
            double Umax = -system.getQMin() / (2. * system.getMass());
            OutputManager.write(String.format("Integrator is set to %s, and this integrator needs potential capping.", integrator.getClass().getSimpleName()));
            OutputManager.write(String.format("The potential Q(x) = 2m [E-U(r)] will be capped from below to %15.4e.", system.getQMin()));
            OutputManager.write(String.format("This means U(r) < E + %15.4e Eh (%15.4e cm-1)", Umax, toInverseCm(Umax)));
            OutputManager.writeBlankLine();
            level.setCapPotential(true);
        } else {
            OutputManager.write(String.format("Integrator is set to %s, and this integrator does NOT needs potential capping.",
                    integrator.getClass().getSimpleName()));
        }


        double energyScale = estimateEnergyScaleAndLowerBound(level);

        // 3. Exponential Scan to find upper bound to the energy
        double currentEnergy = level.energy;

        QuantumLevel.ConvergenceInfo info = new QuantumLevel.ConvergenceInfo();
        info.convergengeStage = "Initial energy bracketing (find upper bound)";
        info.iterations = 0;

        OutputManager.writeBlankLine();
        for (int i = 0; i < MAXIMUM_NUMBER_OF_BISECTIONS; i++) {
            // Update the energy in the system
            level.energy = currentEnergy;
            int nodes = countNodes(level);
            info.iterations++;
            ;

            for (int k = 0; k < level.psi.length; k++) {
                OutputManager.writeData(String.format("%8d %20.6e", k, level.psi[k]));
            }

            if (nodes > nOfDesiredNodes) {
                level.upperBound = currentEnergy;
                level.nodesUpper = nodes;
                level.energy = 0.5 * (level.lowerBound + level.upperBound);
                OutputManager.write(String.format("Upper bound found, current energy set to %s", fmtEnergy(level.energy)));
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

        OutputManager.write(String.format("I scanned the potential and found a minimum value %23.14f (%25.6f cm-1) for i = %d",
                uMin, toInverseCm(uMin), minIndex));

        // 2. Estimate Step Size (Energy Scale)
        Double energyScale = null;

        // Attempt A: Harmonic Curvature of U_tilde on uniform y-grid
        String energyScaleMethod = "";
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
                energyScaleMethod = "harmonic constant at equilibrium";
            }
        }

        // Attempt B: Particle-in-a-Box (Fallback)
        if (energyScale == null) {
            double L = grid.getLastYValue() - grid.getFirstYValue();
            energyScale = (Math.PI * Math.PI) / (2.0 * system.getMass() * L * L);
            energyScaleMethod = "Particle-in-a-Box";
        }

        OutputManager.write(String.format("The energy scale was set to %23.14f (%25.6f cm-1) using as method: %s",
                energyScale, toInverseCm(energyScale), energyScaleMethod));

        level.lowerBound = uMin - energyScale * 0.05; // Set minimum a bit lower than minimum of the potential on the grid.
        level.nodesLower = 0; // Should be always correct

        level.energy = uMin + energyScale;
        return energyScale;
    }

    private QuantumLevel findEigenvalueByBisection(int nOfDesiredNodes) {

        QuantumLevel level = this.findInitialEnergyBracket(nOfDesiredNodes);
        refineByBisection(level, nOfDesiredNodes, TARGET_ABSOLUTE_ERROR, 0);
        level.normalizePsi();
        if (corrector != null) {
            corrector.computeAndSet(level);
        }

        return level;
    }

    private void refineByBisection(QuantumLevel level, int nOfDesiredNodes, double maxAbsError, int minBisections) {
        QuantumLevel.ConvergenceInfo info1 = new QuantumLevel.ConvergenceInfo();
        info1.convergengeStage = "Bisection (to strict bracketing)";
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
        OutputManager.write("************************************************************************************");
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
        QuantumLevel l = this.findInitialEnergyBracket(nOfDesiredNodes);
        OutputManager.write(String.format("Initial LOWER energy is: %23.14f (%25.6f cm-1)", l.lowerBound, toInverseCm(l.lowerBound)));
        OutputManager.write(String.format("Initial GUESS energy is: %23.14f (%25.6f cm-1)", l.energy, toInverseCm(l.energy)));
        OutputManager.write(String.format("Initial UPPER energy is: %23.14f (%25.6f cm-1)", l.upperBound, toInverseCm(l.upperBound)));
        refineByBisection(l, nOfDesiredNodes, 1e-2, 3);
        OutputManager.write("Initial bisection steps finished. The new brackets are:");
        OutputManager.write(String.format("Initial LOWER energy is: %23.14f (%25.6f cm-1)", l.lowerBound, toInverseCm(l.lowerBound)));
        OutputManager.write(String.format("Initial GUESS energy is: %23.14f (%25.6f cm-1)", l.energy, toInverseCm(l.energy)));
        OutputManager.write(String.format("Initial UPPER energy is: %23.14f (%25.6f cm-1)", l.upperBound, toInverseCm(l.upperBound)));
        refineByBidirectionalMatching(l);
        l.normalizePsi();
        if (corrector != null) {
            corrector.computeAndSet(l);
        }
        return l;
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

        OutputManager.writeBlankLine();
        OutputManager.write(String.format("Energy refinement stage. iterations = %d", info.iterations));
        OutputManager.write(String.format("Derivative mismatch for UPPER energy: %25.12f (%25.8f cm-1/a0)", diffUpper, toInverseCm(diffUpper)));
        OutputManager.write(String.format("Derivative mismatch for LOWER energy: %25.12f (%25.8f cm-1/a0)", diffLower, toInverseCm(diffLower)));

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

            OutputManager.writeBlankLine();
            OutputManager.write(String.format("Energy refinement stage. iterations = %d", info.iterations));
            OutputManager.write(String.format("Current energy is: %25.12f (%25.8f cm-1)", x2, toInverseCm(x2)));
            OutputManager.write(String.format("Derivative mismatch for  current energy: %25.12f (%25.8f cm-1/a0)", f2, toInverseCm(f2)));

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
        int startIndex = initialize(level, Integrator.Direction.FORWARD);

        for (int n = startIndex; n < matchIndex + 1; n++) {
            level.psi[n + 1] = integrator.propagate(level.psi, level.currentPsiPrime, n, hy, qTildeFunction,
                    level::QTildePrimeAtGridPoint, level::QTildeDoublePrimeAtGridPoint,
                    Integrator.Direction.FORWARD);

            // Check for potential overflow
            if (n % 16 == 0 && Math.abs(level.psi[n + 1]) > PSI_MAX) {
                double factor = level.psi[n + 1];
                for (int i = 0; i <= n + 1; i++) { // Rescale computed points
                    level.psi[i] /= factor;
                }
                level.currentPsiPrime[0] /= factor;
            }

        }

        double forwardDer = (level.psi[matchIndex + 1] - level.psi[matchIndex - 1]) / (level.psi[matchIndex] * 2. * hy);
        double forwardPsiAtMatchIndexMinusOne = level.psi[matchIndex - 1];
        double forwardPsiAtMatchIndex = level.psi[matchIndex];

        // --- Shoot Backward
        int np = system.getGrid().getNumberOfPoints();
        startIndex = initialize(level, Integrator.Direction.BACKWARD);

        for (int n = startIndex; n > matchIndex - 1; n--) {
            level.psi[n - 1] = integrator.propagate(level.psi, level.currentPsiPrime, n, hy,
                    qTildeFunction, level::QTildePrimeAtGridPoint, level::QTildeDoublePrimeAtGridPoint,
                    Integrator.Direction.BACKWARD);

            // Check for potential overflow
            if (n % 16 == 0 && Math.abs(level.psi[n - 1]) > PSI_MAX) {
                double factor = level.psi[n + 1];
                for (int i = level.psi.length - 1; i >= n - 1; i--) { // Rescale computed points
                    level.psi[i] /= factor;
                }
                level.currentPsiPrime[0] /= factor;
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

        // --- Shoot Forward
        int startIndex = initialize(level, Integrator.Direction.FORWARD);

        for (int n = startIndex; n < nPoints - 1; n++) {
            level.psi[n + 1] = integrator.propagate(level.psi, level.currentPsiPrime, n, hy, level::QTildeAtGridPoint,
                    level::QTildePrimeAtGridPoint, level::QTildeDoublePrimeAtGridPoint,
                    Integrator.Direction.FORWARD);

            // Check for potential overflow
            if (n % 16 == 0 && Math.abs(level.psi[n + 1]) > PSI_MAX) {
                double factor = level.psi[n + 1];
                for (int i = 0; i <= n + 1; i++) { // Rescale computed points
                    level.psi[i] /= factor;
                }
                level.currentPsiPrime[0] /= factor;
            }

        }
        return level.countNodes();
    }

    private int findMatchingIndex(double energy) {
        Grid grid = system.getGrid();
        // Scan from right to left until we reach the classically-allowed region
        int minIndex = Math.max(2, integrator.minHistoryLength());
        for (int i = grid.getNumberOfPoints() - (minIndex + 1); i >= minIndex; i--) {
            if (system.UTildeAtGridPoint(i) <= energy) {
                return i;
            }
        }

        // Fallback to midpoint if always classically forbidden
        return grid.getNumberOfPoints() / 2;
    }

    private int initialize(QuantumLevel level, Integrator.Direction direction) {
        int d = direction.getValue();
        int i = (d == 1) ? 0 : level.psi.length - 1;
        level.psi[i] = 0.;
        level.currentPsiPrime[0] = d * 1.e-16;

        int minHistory = integrator.minHistoryLength();
        if (minHistory == 1) {
            return i;
        }

        i += d;
        level.psi[i] = 1e-16;
        if (minHistory == 2) {
            return i;
        }

        // Initialize additional points for multi-step integrators
        double hy = system.getGrid().getStepSizeYCoordinate();
        for (int k = 3; k <= minHistory; k++) {
            double nextValue = integratorBest.propagate(level.psi, level.currentPsiPrime, i, hy, level::QTildeAtGridPoint,
                    level::QTildePrimeAtGridPoint, level::QTildeDoublePrimeAtGridPoint, direction);
            i += d;
            level.psi[i] = nextValue;
        }

        return i;

    }

    private String fmtEnergy(double energy) {
        return String.format("%25.14f (%25.6f cm-1)", energy, toInverseCm(energy));
    }

}
