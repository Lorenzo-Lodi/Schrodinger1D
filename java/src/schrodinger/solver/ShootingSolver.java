package schrodinger.solver;

import schrodinger.OutputManager;
import schrodinger.QuantumLevel;
import schrodinger.integrator.IntegratorFactory;
import schrodinger.integrator.pt_correction.PTCorrector;
import schrodinger.grid.Grid;
import schrodinger.integrator.Integrator;
import schrodinger.SchrodingerSystem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.DoubleUnaryOperator;

import static schrodinger.PhysicalConstants.toHartree;
import static schrodinger.PhysicalConstants.toInverseCm;

public class ShootingSolver {
    private static final double TARGET_ABSOLUTE_ERROR = 1e-13;
    private static final int MAXIMUM_NUMBER_OF_BISECTIONS = 60; // reduces error by 2**n
    private static final double PSI_MAX = 1e140; // rescale if psi exceeds this value
    private final Integrator integrator;
    private final SchrodingerSystem system;
    private final Integrator integratorBest = IntegratorFactory.getBestOneStepIntegrator();
    private RefinementStrategy strategy;
    public Bounds bounds;

    public ShootingSolver(SchrodingerSystem system, Integrator integrator) {
        this.system = system;
        this.integrator = integrator;
        this.system.initializeCache(integrator.getFractionalOffsets());
    }

    public List<QuantumLevel> findEigenvaluesUpTo(int maximumQuantumNumber, RefinementStrategy refinementStrategy) {
        List<QuantumLevel> quantumLevels = new ArrayList<>();
        bounds = new Bounds(maximumQuantumNumber);


        for (int v = maximumQuantumNumber; v >= 0; v--) {
            QuantumLevel level = findEigenvalue(v, refinementStrategy);
            quantumLevels.add(level);
        }

        return quantumLevels;
    }

    public QuantumLevel findEigenvalue(int nOfDesiredNodes, RefinementStrategy refinementStrategy) {
        OutputManager.write("************************************************************************************");
        OutputManager.write(String.format("Finding eigenvalue with %d nodes, strategy %s", nOfDesiredNodes, refinementStrategy.toString()));
        this.strategy = refinementStrategy;
        QuantumLevel level = this.findInitialEnergyBracket(nOfDesiredNodes);

        switch (strategy) {
            case BISECTION_ONLY:
                findEigenvalueByBisection(level, nOfDesiredNodes);
                break;
            case BISECTION_THEN_BIDIRECTIONAL:
                findEigenvalueByHybridMethod(level, nOfDesiredNodes);
                break;
        }

        level.normalizePsi();
        PTCorrector corrector = integrator.getPertubativeCorrector();
        if (corrector != null) {
            corrector.computeAndSet(level);
        }

        return level;
    }

    public QuantumLevel findEigenvalue(int nOfDesiredNodes) {
        return findEigenvalue(nOfDesiredNodes, RefinementStrategy.BISECTION_THEN_BIDIRECTIONAL);
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
            OutputManager.write(String.format("Integrator is set to %s, and this integrator needs potential capping.", integrator.getClass().getSimpleName()));
            OutputManager.write(String.format("The potential Q(x) = 2m [E-U(r)] will be capped from below to %15.4e.", system.getQMin()));
            OutputManager.write(String.format("This means U(r) < E + %15.4e Eh (%15.4e cm-1)", system.Ucapped(), toInverseCm(system.Ucapped())));
            OutputManager.writeBlankLine();
            level.setCapPotential(true);
        } else {
            OutputManager.write(String.format("Integrator is set to %s, and this integrator does NOT needs potential capping.",
                    integrator.getClass().getSimpleName()));
        }
        double energyScale = system.estimateEnergyScale();

        // Set minimum a bit lower than minimum of the potential on the grid
        level.lowerBound = system.UTildeMinimumGridValue - energyScale * 0.05;
        level.nodesLower = 0; // Should be always correct
        bounds.updateBounds(level.lowerBound, level.nodesLower);
        bounds.updateBounds(system.uMaxRight, nOfDesiredNodes);

        // 3. Exponential Scan to find upper bound to the energy
        double currentEnergy = level.lowerBound + energyScale * (nOfDesiredNodes + 1);

        QuantumLevel.ConvergenceInfo info = new QuantumLevel.ConvergenceInfo();
        info.convergengeStage = "Initial energy bracketing (find upper bound)";
        info.iterations = 0;
        level.convergenceInfo.add(info);

        OutputManager.writeBlankLine();
        for (int i = 0; i < MAXIMUM_NUMBER_OF_BISECTIONS; i++) {
            // Update the energy in the system
            level.energy = currentEnergy;
            int nodes = countNodes(level);
            bounds.updateBounds(level.energy, nodes);
            {
                String msg = ((nodes - nOfDesiredNodes) < 0) ? " too few" : " OK";
                OutputManager.write(String.format("Trying to find an upper bound. Current energy = %20.6f cm-1, nodes = %10d, " +
                                "(nodes - nOfDesiredNodes) = %10d, so nodes is %s",
                        toInverseCm(currentEnergy), nodes, (nodes - nOfDesiredNodes), msg));
            }
            info.iterations++;

            if (nodes > nOfDesiredNodes) {
                level.upperBound = currentEnergy;
                level.nodesUpper = nodes;
                level.energy = 0.5 * (level.lowerBound + level.upperBound);
                OutputManager.write(String.format("Upper bound found, current energy set to %s", fmtEnergy(level.energy)));
                level.verifyStepSize();
                return level;
            } else {
                level.lowerBound = currentEnergy;
                level.nodesLower = nodes;
                currentEnergy += energyScale;
                energyScale *= 1.5;
            }
        }

        throw new RuntimeException("Failed to bracket energy level.");
    }

    private void findEigenvalueByBisection(QuantumLevel level, int nOfDesiredNodes) {
        refineByBisection(level, nOfDesiredNodes, TARGET_ABSOLUTE_ERROR, 0);
    }

    private void refineByBisection(QuantumLevel level, int nOfDesiredNodes, double maxAbsError, int extraBisections) {
        QuantumLevel.ConvergenceInfo info1 = new QuantumLevel.ConvergenceInfo();
        info1.convergengeStage = "Bisections (to strict bracketing)";
        info1.iterations = 0;
        level.convergenceInfo.add(info1);
        QuantumLevel.ConvergenceInfo info2 = new QuantumLevel.ConvergenceInfo();
        info2.convergengeStage = "Bisections (after strict bracketing)";
        info2.iterations = 0;
        level.convergenceInfo.add(info2);


        // Tentative code for linear interpolation-search (should be better)
//        {
//            double span = level.nodesUpper - level.nodesLower;
//            double fraction = (nodes - level.nodesLower) / span;
//            double epsilon = 0.5 / span;
//            fraction = Math.clamp(fraction, epsilon, 1.0 - epsilon);
//            level.energy = level.lowerBound + fraction * (level.upperBound - level.lowerBound);
//        }

        for (int i = 1; i <= MAXIMUM_NUMBER_OF_BISECTIONS; i++) {

            level.energy = (level.lowerBound + level.upperBound) * 0.5;
            int nodes = countNodes(level);

            OutputManager.write(String.format("Bisection refinement. i = %5d; nodes= %5d %22.10f %22.10f %22.10f", i, nodes,
                    toInverseCm(level.lowerBound), toInverseCm(level.energy), toInverseCm(level.upperBound)));

            if (nodes > nOfDesiredNodes) {
                level.upperBound = level.energy;
                level.nodesUpper = nodes;
            } else {
                level.lowerBound = level.energy;
                level.nodesLower = nodes;
            }

            if (level.nodesLower == nOfDesiredNodes && level.nodesUpper == nOfDesiredNodes + 1) {
                info2.iterations++; // We have only one state (the desired one) in the energy bracket
            } else {
                info1.iterations++;
            }

            if (Math.abs(level.upperBound - level.lowerBound) < maxAbsError &&
                    info2.iterations >= extraBisections) {
                break;
            }

        }
        level.energy = (level.lowerBound + level.upperBound) * 0.5;
    }

    private void findEigenvalueByHybridMethod(QuantumLevel l, int nOfDesiredNodes) {
        int maxMacroIterations = 3;
        int matchIndex = -1;
        double diffLower = 0;
        double diffUpper = 0;

        // Loop:
        // 1. Do some bisection
        // 2. Check if sign of derivative differences for upper/lower energy have different signs
        // 3. different sign => break loop and start refinement with improved regula falsi
        // 4. loop again (more bisection and then test again the signs)

        for (int i = 1; i <= maxMacroIterations; i++) {
            OutputManager.write(String.format("Initial LOWER energy is: %23.14f (%25.6f cm-1)", l.lowerBound, toInverseCm(l.lowerBound)));
            OutputManager.write(String.format("Initial GUESS energy is: %23.14f (%25.6f cm-1)", l.energy, toInverseCm(l.energy)));
            OutputManager.write(String.format("Initial UPPER energy is: %23.14f (%25.6f cm-1)", l.upperBound, toInverseCm(l.upperBound)));

            refineByBisection(l, nOfDesiredNodes, toHartree(1.0), 4);

            OutputManager.write("Bisection finished. The new brackets are:");
            OutputManager.write(String.format("LOWER energy is: %23.14f (%25.6f cm-1)", l.lowerBound, toInverseCm(l.lowerBound)));
            OutputManager.write(String.format("GUESS energy is: %23.14f (%25.6f cm-1)", l.energy, toInverseCm(l.energy)));
            OutputManager.write(String.format("UPPER energy is: %23.14f (%25.6f cm-1)", l.upperBound, toInverseCm(l.upperBound)));

            // It seems preferable to compute the matching index once and for all
            matchIndex = findMatchingIndex(l.energy);
            OutputManager.write(String.format("Energy refinement stage (initialization). Macro iteration = %10d", i));
            OutputManager.write(String.format("Matching index set to %10d", matchIndex));

            // Here we compute the current derivative difference for the upper/lower energies
            QuantumLevel.ConvergenceInfo info = new QuantumLevel.ConvergenceInfo();
            info.convergengeStage = "Compute derivative differences to see if we can switch away from bisection...";
            info.iterations = 0;
            l.convergenceInfo.add(info);
            l.energy = l.upperBound;
            diffUpper = computeDerivativeMismatch(l, matchIndex);
            info.iterations++;

            l.energy = l.lowerBound;
            diffLower = computeDerivativeMismatch(l, matchIndex);
            info.iterations++;

            OutputManager.writeBlankLine();
            OutputManager.write(String.format("Derivative mismatch for UPPER energy: %25.12f (%25.8f cm-1/a0)", diffUpper, toInverseCm(diffUpper)));
            OutputManager.write(String.format("Derivative mismatch for LOWER energy: %25.12f (%25.8f cm-1/a0)", diffLower, toInverseCm(diffLower)));

            if (diffUpper * diffLower < 0) {
                OutputManager.write(String.format("Macroiteration i = %10d, the sign of the derivative mismatch for the " +
                        "upper and lower energy are different, I can start iterating the fast root-finding algorithm", i));
                break;
            } else {
                if (i < maxMacroIterations) {
                    OutputManager.write(String.format("Macroiteration i = %10d, the sign of the derivative mismatch for the " +
                            "upper and lower energy are THE SAME, I will try some more bisection and try again!", i));
                } else {
                    OutputManager.write(String.format("Macroiteration i = %10d, the sign of the derivative mismatch for the " +
                            "upper and lower energy are THE SAME, and I've reached maxMacroIterations! I WILL RETURN NOW", i));
                    return;
                }
            }
        }

        refineByBidirectionalMatching(l, matchIndex, diffLower, diffUpper);

    }

    private void refineByBidirectionalMatching(QuantumLevel level, int matchIndex, double diffLower, double diffUpper) {

        QuantumLevel.ConvergenceInfo info = new QuantumLevel.ConvergenceInfo();
        info.convergengeStage = "Refinement by " + strategy.toString().toLowerCase().replace("bisection_then_", "");
        info.iterations = 0;
        level.convergenceInfo.add(info);

        int maxIter = 50;   // prevent infinite loops
        double tol = 1e-11; // tolerance on derivative difference (value appropriate when log-der is used)
//        double tol = 1e-18; // tolerance on derivative difference (value appropriate when wroksian is used)

        double diffEnergy;
        for (int iter = 0; iter < maxIter; iter++) {
            // Compute the false position point directly into level.energy
            level.energy = level.upperBound - diffUpper * (level.upperBound - level.lowerBound) / (diffUpper - diffLower);

            // Evaluate function
            diffEnergy = computeDerivativeMismatch(level, matchIndex);
            info.iterations++;

            OutputManager.writeBlankLine();
            OutputManager.write(String.format("Energy refinement stage. iterations = %d", info.iterations));
            OutputManager.write(String.format("Current energy is: %25.12f (%25.8f cm-1)", level.energy, toInverseCm(level.energy)));
            OutputManager.write(String.format("Derivative mismatch for current energy: %25.12f (%25.8f cm-1/a0)", diffEnergy, toInverseCm(diffEnergy)));

            // Check for convergence on the derivative difference
            if (Math.abs(diffEnergy) < tol) {
                break;
            }

            // Modified regula falsi update logic (Illinois, Pegasus, Anderson-Björck)
            if (diffLower * diffEnergy < 0) {
                // Root lies between lowerBound and energy. Endpoint lowerBound is retained.
                level.upperBound = level.energy;
                diffLower = diffLower * computeM(diffEnergy, diffUpper);
                diffUpper = diffEnergy;
            } else {
                // Root lies between energy and upperBound. Endpoint upperBound is retained.
                level.lowerBound = level.energy;
                diffUpper = diffUpper * computeM(diffEnergy, diffUpper);
                diffLower = diffEnergy;
            }
        }

        // TODO rescale wavefunction from zero to matchIndex

    }

    private double computeM(double fEnergy, double fReplaced) {
//  Illinois or Pegasus may be more stable. Consider switching to either of them if , eg, iter reaches 20 or so.
//         return 0.5;  // Illinois method.
//        return fReplaced / (fReplaced + fEnergy); // Pegasus method.
        double m = 1.0 - fEnergy / fReplaced;  // Anderson-Björck method.
        return (m <= 0.0) ? 0.5 : m;
    }

    private double computeDerivativeMismatch(QuantumLevel level, int matchIndex) {
        double hy = system.getGrid().getStepSizeYCoordinate();
        DoubleUnaryOperator qTildeFunction = level::QTildeAtGridPoint;

        // --- Shoot Forward
        Arrays.fill(level.psi, 0.0); // Let us zero the wave function for clarity (not necessary).
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

        double forwardDer = (level.psi[matchIndex + 1] - level.psi[matchIndex - 1]) / (2. * hy);
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
                double factor = level.psi[n - 1];
                for (int i = level.psi.length - 1; i >= n - 1; i--) { // Rescale computed points
                    level.psi[i] /= factor;
                }
                level.currentPsiPrime[0] /= factor;
            }

        }
        double backwardDer = (level.psi[matchIndex + 1] - level.psi[matchIndex - 1]) / (2. * hy);
        double backwardPsiAtMatchIndex = level.psi[matchIndex];

        // Let us rescale the correct psi (probably unnecessary doing this at each step). It should be done only once
        // after convergence is reached.
        level.psi[matchIndex - 1] = forwardPsiAtMatchIndexMinusOne;
        for (int n = 0; n < matchIndex; n++) {
            level.psi[n] = level.psi[n] / forwardPsiAtMatchIndex;
        }
        for (int n = matchIndex + 1; n < np - 1; n++) {
            level.psi[n] = level.psi[n] / level.psi[matchIndex];
        }
        level.psi[matchIndex] = 1.;

        double result = forwardDer / forwardPsiAtMatchIndex - backwardDer / backwardPsiAtMatchIndex; // Return logarithmic derivative
//        double result = forwardDer * backwardPsiAtMatchIndex - backwardDer * forwardPsiAtMatchIndex; // Return Wroksian; seems to give problems!!!!
//        result = result / (Math.max(Math.abs(forwardPsiAtMatchIndex), Math.abs(backwardPsiAtMatchIndex)));
        return result;
    }

    private int countNodes(QuantumLevel level) {
        int nPoints = system.getGrid().getNumberOfPoints();
        double hy = system.getGrid().getStepSizeYCoordinate();

        // --- Shoot Forward
        int startIndex = initialize(level, Integrator.Direction.FORWARD);
        int rightmostInversionPoint = findMatchingIndex(level.energy);
        if (rightmostInversionPoint == level.getGrid().getNumberOfPoints() / 2) {
            rightmostInversionPoint = 0;
        }

        double gamma = 0.0; // To track how deep into the forbidden region we are
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

            // When we start integrating into the final forbidden region, keep track of gamma
            if (n > rightmostInversionPoint) {
                gamma += hy * Math.sqrt(-level.QTildeAtGridPoint(n));
            }

            if (gamma > 16.0) { // Safety margin before the 18.4 wall
                // Stop outward integration! The wavefunction has decayed by a factor of e^-16 (~1e-7)
                // and numerical noise is about to take over.
                // Let us fill the wavefunction with extrapolated values (maybe not necessary, but it could be useful
                // for correct node-counting and for later computing matrix elements with this wavefunction.
                double decayFactor = level.psi[n] / level.psi[n + 1];
                OutputManager.write("Stopping integrating into the forbidden region because the gamma threshold was reached" +
                        " at n = " + n + "; extrapolating the remaining points using decayFactor = " + decayFactor);
                for (int k = n + 1; k < nPoints - 1; k++) {
                    level.psi[k + 1] = level.psi[k] * decayFactor;
                }
                break;
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
