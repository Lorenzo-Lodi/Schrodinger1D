package schrodinger.integrator.verification;

import schrodinger.QuantumLevel;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.Integrator;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.SchrodingerSystem;

import java.util.HashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class HarmonicOscillatorConvergenceVerifier {

    private final static int INITIALIZATION_N_MAX = 7;
    private static final double DEFAULT_R0 = 0.0;
    private static final double DEFAULT_ALPHA = 1.0;
    private HarmonicOscillatorExactSolution exactSolutionInstance;
    private final Integrator integrator;
    private final Double maxErrorThreshold;

    public HarmonicOscillatorConvergenceVerifier(Integrator integrator, Double maxErrorThreshold) {
        this.integrator = integrator;
        this.maxErrorThreshold = maxErrorThreshold;
    }

    /**
     * Returns the grid boundaries for testing, which depend on the quantum state.
     *
     * @param quantumNumber the quantum number
     * @return array {rMin, rMax}
     */
    private double[] getGridBoundaries(int quantumNumber) {
        double energy = 0.5 + quantumNumber;
        double turningPoint = Math.sqrt(energy);
        return new double[]{DEFAULT_R0 - turningPoint - 3.0, DEFAULT_R0 + turningPoint + 3.0};
    }

    /**
     * Represents the convergence data at different points in the grid.
     * Uses a map to store error values for different fractions of the grid.
     */
    private static class ConvergenceData {
        private final Map<String, Double> errors = new HashMap<>();
        private long elapsedNanos;

        public void addError(String pointFraction, double error) {
            errors.put(pointFraction, error);
        }

        public double getError(String pointFraction) {
            return errors.getOrDefault(pointFraction, 0.0);
        }

        public void setElapsedNanos(long ns) {
            this.elapsedNanos = ns;
        }

        public long getElapsedNanos() {
            return elapsedNanos;
        }
    }

    /**
     * Tests the convergence of the integrator against the exact analytical solution
     * for a specified quantum state.
     *
     * @param quantumNumber the quantum number n (0 = ground state, 1 = first excited, etc.)
     */
    public void verify_harmonic_oscillator(int quantumNumber) {
        if (integrator.minHistoryLength() > INITIALIZATION_N_MAX + 1) {
            String className = this.getClass().getSimpleName();
            String msg = String.format("Integrator %s requires at least %d previously-computed points, but only %d are available!",
                    className, integrator.minHistoryLength(), INITIALIZATION_N_MAX);
            throw new RuntimeException(msg);
        }

        PhysicalPotential potential = new PhysicalPotentialHarmonic(DEFAULT_R0, DEFAULT_ALPHA);
        exactSolutionInstance = new HarmonicOscillatorExactSolution(DEFAULT_R0, DEFAULT_ALPHA, quantumNumber);

        Map<Integer, ConvergenceData> convergenceResults = new HashMap<>();

        // Test with different numbers of grid points
        // NOTE: nOfPoints should be divisible by 20, otherwise when computing nOfPoints/4 (etc.) we get fractions etc.
        for (int nOfPoints = 200; nOfPoints <= 800; nOfPoints += 20) {
            ConvergenceData data = testWithPoints(potential, integrator, nOfPoints, quantumNumber);
            convergenceResults.put(nOfPoints, data);
        }

        // Print results
        String stateName = (quantumNumber == 0) ? "Ground State" : quantumNumber + "th Excited State";
        System.out.println("\n=== Testing " + stateName + " (n=" + quantumNumber + ") ===\n");
        System.out.println("Results for: " + integrator.getClass().getSimpleName() + "\n");
        double[] bounds = getGridBoundaries(quantumNumber);
        System.out.printf("xmin = %25.14f, xmax = %25.14f, grid size = %25.14f\n%n", bounds[0], bounds[1], bounds[1] - bounds[0]);

        printResults(convergenceResults);

        // Calculate and print convergence rates
        printLinearFitResultsAndAssert(convergenceResults, quantumNumber, integrator);
    }

    /**
     * Tests the integrator with a specific number of grid points for a given quantum state.
     *
     * @param potential     The potential to use
     * @param integrator    The integrator to test
     * @param nOfPoints     The number of grid points
     * @param quantumNumber The quantum number for the state being tested
     * @return Convergence data containing errors at different points
     */
    private ConvergenceData testWithPoints(PhysicalPotential potential, Integrator integrator, int nOfPoints, int quantumNumber) {
        double[] gridBounds = getGridBoundaries(quantumNumber);
        double rMin = gridBounds[0];
        double rMax = gridBounds[1];

        Grid grid = GridFactory.generateUniformGrid(rMin, rMax, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, 2., grid);
        QuantumLevel level = new QuantumLevel(system);

        // Use the exact energy for the specified quantum level
        HarmonicOscillatorExactSolution exactSol = new HarmonicOscillatorExactSolution(DEFAULT_R0, DEFAULT_ALPHA, quantumNumber);
        level.energy = exactSol.getEnergy();

        // Initialize with exact solution
        for (int n = 0; n <= INITIALIZATION_N_MAX && n < nOfPoints; n++) {
            level.psi[n] = exactSol.evaluate(grid.rAtGridPoint(n));
        }

        {
            level.currentPsiPrime = new double[1]; // Initialize first derivative for RKN methods
            double eps = 5e-4;
            double r = grid.rAtGridPoint(INITIALIZATION_N_MAX);
            double psip1 = exactSol.evaluate(r + eps);
            double psim1 = exactSol.evaluate(r - eps);
            double der1 = (psip1 - psim1) / (2. * eps);

            double psipp1 = exactSol.evaluate(r + 2. * eps);
            double psimm1 = exactSol.evaluate(r - 2. * eps);
            double der2 = (psipp1 - psimm1) / (4. * eps);

            double der = (4. / 3.) * der1 - (1. / 3.) * der2;

            level.currentPsiPrime[0] = der;
        }

        // Propagate the wavefunction
        double step = level.getGrid().getStepSizeYCoordinate();
        long t0 = System.nanoTime();
        for (int n = INITIALIZATION_N_MAX; n < nOfPoints - 1; n++) {
            level.psi[n + 1] = integrator.propagate(level.psi, level.currentPsiPrime, n, step,
                    level::QTildeAtGridPoint,
                    level::QTildePrimeAtGridPoint, level::QTildeDoublePrimeAtGridPoint,
                    Integrator.Direction.FORWARD);
        }
        long elapsedNanos = System.nanoTime() - t0;

        // Calculate errors at different points
        ConvergenceData data = new ConvergenceData();
        data.setElapsedNanos(elapsedNanos);
        data.addError("10", exactSolution(grid.rAtGridPoint(nOfPoints / 10.)) - level.psi[nOfPoints / 10]);
        data.addError("4", exactSolution(grid.rAtGridPoint(nOfPoints / 4.)) - level.psi[nOfPoints / 4]);
        data.addError("2", exactSolution(grid.rAtGridPoint(nOfPoints / 2.)) - level.psi[nOfPoints / 2]);

        return data;
    }

    /**
     * Prints the convergence results including point-to-point convergence rates.
     */
    private void printResults(Map<Integer, ConvergenceData> convergenceResults) {
        System.out.println("Grid Points       Error at 10%            Error at 25%            Error at 50%            Rate 10%   Rate 25%   Rate 50%   Time (ms)");
        System.out.println("-----------------------------------------------------------------------------------------------------------------------------------------");

        // Convert to sorted list for easier access to previous entry
        java.util.List<Map.Entry<Integer, ConvergenceData>> sortedEntries = convergenceResults.entrySet().stream()
                .sorted(Map.Entry.comparingByKey())
                .collect(java.util.stream.Collectors.toList());

        for (int i = 0; i < sortedEntries.size(); i++) {
            Map.Entry<Integer, ConvergenceData> entry = sortedEntries.get(i);
            int nOfPoints = entry.getKey();
            ConvergenceData data = entry.getValue();

            // Print grid points and errors (using % flag to align positive/negative values)
            System.out.printf("%10d\t\t% 20.19f\t% 20.19f\t% 20.19f",
                    nOfPoints,
                    data.getError("10"),
                    data.getError("4"),
                    data.getError("2"));

            // Calculate and print convergence rates if not the first entry
            if (i > 0) {
                Map.Entry<Integer, ConvergenceData> prevEntry = sortedEntries.get(i - 1);
                int prevPoints = prevEntry.getKey();
                ConvergenceData prevData = prevEntry.getValue();

                double rate10 = calculateConvergenceRate(prevData.getError("10"), data.getError("10"), prevPoints, nOfPoints);
                double rate25 = calculateConvergenceRate(prevData.getError("4"), data.getError("4"), prevPoints, nOfPoints);
                double rate50 = calculateConvergenceRate(prevData.getError("2"), data.getError("2"), prevPoints, nOfPoints);

                System.out.printf("\t%9.4f  %9.4f  %9.4f", rate10, rate25, rate50);
            } else {
                System.out.printf("\t%9s  %9s  %9s", "-", "-", "-");
            }

            double ms = data.getElapsedNanos() / 1_000_000.0;
            System.out.printf("  %9.3f%n", ms);
        }

        System.out.println();

        long totalNanos = convergenceResults.values().stream()
                .mapToLong(ConvergenceData::getElapsedNanos).sum();
        System.out.printf("Total propagation time: %.3f ms%n%n", totalNanos / 1_000_000.0);
    }

    /**
     * Computes and prints linear fit parameters for convergence analysis.
     * Fits ln(|error|) = a + b * ln(nPoints) using least squares regression.
     *
     * @param convergenceResults the convergence data for different grid sizes
     * @param quantumNumber      the quantum number being tested
     */
    private void printLinearFitResultsAndAssert(Map<Integer, ConvergenceData> convergenceResults, int quantumNumber, Integrator integrator) {
        String className = this.getClass().getSimpleName();
        String stateName = (quantumNumber == 0) ? "Ground State" : "n=" + quantumNumber + " Excited State";
        System.out.println("Linear Fit Results for " + className + " (" + stateName + "):");
        System.out.println("Fitting ln(|error|) = a + b * ln(nPoints)");
        System.out.println("--------------------------------------------------------");

        String[] pointFractions = {"10", "4", "2"};
        String[] pointNames = {"10%", "25%", "50%"};
        double expectedBCoefficient = -integrator.globalConvergenceOrder();

        for (int idx = 0; idx < pointFractions.length; idx++) {
            String fraction = pointFractions[idx];
            String pointName = pointNames[idx];

            // Get point-specific and state-dependent tolerances
            double minRSquared = getMinRSquared(quantumNumber, pointName);
            double convergenceTol = getConvergenceOrderTolerance(quantumNumber, pointName);

            // Create arrays for ln(points) and ln(|error|)
            int n = convergenceResults.size();
            double[] lnPoints = new double[n];
            double[] lnErrors = new double[n];

            int i = 0;
            for (Map.Entry<Integer, ConvergenceData> entry : convergenceResults.entrySet()) {
                int nOfPoints = entry.getKey();
                double error = entry.getValue().getError(fraction);

                if (error == 0.) { // small hack to avoid exact zero
                    error = 1e-19;
                }

                lnPoints[i] = Math.log(nOfPoints);
                lnErrors[i] = Math.log(Math.abs(error));
                i++;
            }

            // Compute linear fit
            double[] fitParams = LinearFit.fit(lnPoints, lnErrors);
            double a = fitParams[0];  // intercept
            double b = fitParams[1];  // slope
            double rSquared = fitParams[2];  // R-squared

            System.out.printf("Point %s:\t\ta = %.6f\tb = %.6f\tR^2 = %.6f%n",
                    pointNames[idx], a, b, rSquared);

            if (maxErrorThreshold == null) {
                // Assert R² is greater than minimum threshold (state-dependent)
                assertTrue(rSquared > minRSquared,
                        String.format("%s (%s): R² (%.6f) should be > %.6f for point %s",
                                className, stateName, rSquared, minRSquared, pointNames[idx]));

                // Assert b coefficient is within tolerance of expected value (state-dependent)
                double bDeviation = Math.abs(b - expectedBCoefficient);
                assertTrue(bDeviation <= convergenceTol,
                        String.format("%s (%s): b coefficient (%.6f) deviates %.6f from expected %.6f (max allowed: %.6f) for point %s",
                                className, stateName, b, bDeviation, expectedBCoefficient, convergenceTol, pointNames[idx]));
            } else {
                // if maxErrorThreshold specified, just check the errors
                for (int j = 0; j < lnErrors.length - 1; j++) {
                    double error = Math.exp(lnErrors[j]);
                    assertTrue(error < maxErrorThreshold,
                            String.format("%s (%s): error %20.6e should be <= %20.6e for point %s",
                                    className, stateName, error, maxErrorThreshold, pointNames[idx]));

                }
            }
        }

        System.out.println();
    }

    /**
     * Calculates the convergence rate between two points.
     */
    private double calculateConvergenceRate(double error1, double error2, int n1, int n2) {
        if (error1 == 0 || error2 == 0) return 0.0;

        double denom = Math.log10((double) n2 / n1);
        return Math.log10(Math.abs(error2 / error1)) / denom;
    }

    /**
     * Returns the exact analytical solution at position x using the modular solution.
     *
     * @param x the position
     * @return the wavefunction value ψ(x)
     */
    private double exactSolution(double x) {
        return exactSolutionInstance.evaluate(x);
    }

    /**
     * Returns the minimum R² threshold for linear fit, which depends on quantum state,
     * integrator order, and evaluation point. Ground state can achieve very high R²,
     * but excited states with nodes show more variability, especially at points
     * near nodes (typically the 10% evaluation point).
     *
     * @param quantumNumber the quantum number
     * @param pointName     the evaluation point name (e.g., "10%", "25%", "50%")
     * @return minimum acceptable R²
     */
    private double getMinRSquared(int quantumNumber, String pointName) {
        int baseOrder = integrator.globalConvergenceOrder();

        if (quantumNumber == 0) {
            // Ground state: very high R² achievable for lower order methods
            if (baseOrder < 6) {
                return 0.995;  // Original strict tolerance
            } else {
                return 0.80;   // higher-order methods reach numerical noise for dense grids
            }
        }
        // Excited states: more variation due to nodes
        // The 10% point is often near a node and shows erratic behavior
        // The 50% point (center) can also show variation due to being at peak amplitude
        boolean isNearNode = pointName.equals("10%");
        boolean isCenter = pointName.equals("50%");

        if (baseOrder < 2) {
            if (isNearNode) return 0.70;
            if (isCenter) return 0.60;  // Also relaxed for 2nd order at center
            return 0.94;
        } else if (baseOrder <= 4.0) {
            if (isNearNode) return 0.85;
            if (isCenter) return 0.90;
            return 0.94;
        } else {
            if (isNearNode) return 0.69;
            if (isCenter) return 0.65;  // Very relaxed for high-order at center
            return 0.87;
        }
    }

    /**
     * Returns the convergence order tolerance, which depends on quantum state,
     * integrator order, and evaluation point. Higher-order methods show more variation,
     * and excited states with nodes show additional variation, especially near nodes.
     *
     * @param quantumNumber the quantum number
     * @param pointName     the evaluation point name (e.g., "10%", "25%", "50%")
     * @return maximum allowed deviation from expected convergence order
     */
    private double getConvergenceOrderTolerance(int quantumNumber, String pointName) {
        int baseOrder = integrator.globalConvergenceOrder();

        // Base tolerance depends on integrator order
        double baseTolerance;
        if (baseOrder < 8) {
            baseTolerance = 0.75;
        } else {
            baseTolerance = 2.2; // high order methods are affected by roundoff
        }

        // Add extra tolerance for excited states
        if (quantumNumber > 0 && baseOrder >= 5) {
            baseTolerance = 1.5;
        }

        return baseTolerance;
    }


}
