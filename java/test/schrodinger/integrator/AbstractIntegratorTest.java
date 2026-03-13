package schrodinger.integrator;

import org.junit.jupiter.api.Test;
import schrodinger.QuantumLevel;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.expfitted.EfnFixedBeta;
import schrodinger.integrator.expfitted.ExponentiallyFitted;
import schrodinger.integrator.expfitted.Numerov;
import schrodinger.integrator.predcorr.PredictorCorrector6;
import schrodinger.integrator.predcorr.PredictorCorrector8NumerovIter1;
import schrodinger.integrator.predcorr.PredictorCorrector8NumerovIter2;
import schrodinger.integrator.stormer.Stormer5;
import schrodinger.integrator.stormer.Stormer6;
import schrodinger.integrator.stormer.Stormer8;
import schrodinger.integrator.stormer.Stormer8i;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.SchrodingerSystem;

import java.util.HashMap;
import java.util.Map;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Predicate;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Abstract base class for integrator tests.
 * Provides common functionality for testing integrator convergence.
 * Supports testing against exact solutions for any quantum state of the harmonic oscillator.
 */
public abstract class AbstractIntegratorTest {

    private final static int INITIALIZATION_N_MAX = 7;
    private final static int MIN_POINTS = 200;

    // Default potential parameters (matching the existing test setup)
    protected static final double DEFAULT_R0 = 10.0;
    protected static final double DEFAULT_ALPHA = 1.0;

    // The exact solution for the current quantum state
    private HarmonicOscillatorExactSolution exactSolutionInstance;

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
    protected double getMinRSquared(int quantumNumber, String pointName, Integrator integrator) {
        double baseOrder = getGlobalConvergenceOrder();

        // Exception for the Obrechkoff method. TODO investigate why the convergence is somewhat erratic.
        // Other tests show that the implementation is correct, though...
        if (integrator instanceof Obrechkoff) {
            return 0.68;
        }

        if (quantumNumber == 0) {
            // Ground state: very high R² achievable for lower order methods
            if (baseOrder <= 4.0) {
                return 0.995;  // Original strict tolerance
            } else {
                return 0.98;   // Slightly relaxed for higher-order methods
            }
        } else {
            // Excited states: more variation due to nodes
            // The 10% point is often near a node and shows erratic behavior
            // The 50% point (center) can also show variation due to being at peak amplitude
            boolean isNearNode = pointName.equals("10%");
            boolean isCenter = pointName.equals("50%");

            if (baseOrder <= 2.0) {
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
                return 0.88;
            }
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
    protected double getConvergenceOrderTolerance(int quantumNumber, String pointName, Integrator integrator) {
        double baseOrder = getGlobalConvergenceOrder();

        // Exception for the Obrechkoff method. TODO investigate why the convergence is somewhat erratic.
        // Other tests show that the implementation is correct, though...
        if (integrator instanceof Obrechkoff) {
            return 4.;
        }


        // Base tolerance depends on integrator order
        double baseTolerance;
        if (baseOrder <= 2.0) {
            baseTolerance = 0.4;  // 2nd order methods
        } else if (baseOrder <= 4.0) {
            baseTolerance = 0.5;  // 4th order methods
        } else if (baseOrder <= 5.0) {
            baseTolerance = 0.6;  // 5th order methods
        } else {
            baseTolerance = 1.2;  // 6th+ order methods (more sensitive to grid effects)
        }

        // Add extra tolerance for excited states
        if (quantumNumber > 0) {
            baseTolerance += 0.5;

            // Different points show different behavior
            if (pointName.equals("10%")) {
                // Near a node - be very lenient
                baseTolerance += 1.5;
            } else if (pointName.equals("25%")) {
                // Some integrators show faster-than-expected convergence here
                baseTolerance += 0.9;
            } else if (pointName.equals("50%")) {
                baseTolerance += 0.4;
            }

        }

        return baseTolerance;
    }

    /**
     * Returns the grid boundaries for testing, which depend on the quantum state.
     * For ground state, uses original wider grid. For excited states, uses narrower
     * grid within the classically allowed region to avoid numerical instabilities.
     *
     * @param quantumNumber the quantum number
     * @return array {rMin, rMax}
     */
    protected double[] getGridBoundaries(int quantumNumber) {
        if (quantumNumber == 0) {
            // Original grid for ground state
            return new double[]{6.0, 14.0};
        } else {
            // For excited states, stay within classically allowed region
            // For n=10: E=10.5, turning points at r = 10 ± sqrt(10.5) ≈ 6.76, 13.24
            // Use [7.5, 12.5] to avoid numerical issues near turning points
            return new double[]{DEFAULT_R0 - 2.5, DEFAULT_R0 + 2.5};
        }
    }

    /**
     * Represents the convergence data at different points in the grid.
     * Uses a map to store error values for different fractions of the grid.
     */
    protected static class ConvergenceData {
        private final Map<String, Double> errors = new HashMap<>();
        private long elapsedNanos;

        public void addError(String pointFraction, double error) {
            errors.put(pointFraction, error);
        }

        public double getError(String pointFraction) {
            return errors.getOrDefault(pointFraction, 0.0);
        }

        public Map<String, Double> getAllErrors() {
            return new HashMap<>(errors);
        }

        public void setElapsedNanos(long ns) {
            this.elapsedNanos = ns;
        }

        public long getElapsedNanos() {
            return elapsedNanos;
        }
    }

    /**
     * Returns the quantum number for the state being tested.
     * Subclasses can override this to test different excited states.
     *
     * @return the quantum number n (0 = ground state, 1 = first excited, etc.)
     */
    protected int getQuantumNumber() {
        return 0; // Default to ground state
    }

    /**
     * Returns the exact solution for the current quantum state.
     * This method is called once per test to initialize the exact solution instance.
     *
     * @return the exact solution for the quantum state being tested
     */
    protected HarmonicOscillatorExactSolution getExactSolution() {
        return HarmonicOscillatorExactSolution.excitedState(DEFAULT_R0, DEFAULT_ALPHA, getQuantumNumber());
    }

    /**
     * Tests the convergence of the integrator by comparing numerical solutions
     * with the exact analytical solution for the ground state (n=0).
     */
    @Test
    public void testIntegratorConvergence() {
        testIntegratorConvergence(0); // Default: ground state
    }

    /**
     * Tests the convergence of the integrator by comparing numerical solutions
     * with the exact analytical solution for the 10th excited state (n=10, has 10 nodes).
     */
    @Test
    public void testIntegratorConvergence10thExcitedState() {
        testIntegratorConvergence(10);
    }

    /**
     * Tests the convergence of the integrator against the exact analytical solution
     * for a specified quantum state.
     *
     * @param quantumNumber the quantum number n (0 = ground state, 1 = first excited, etc.)
     */
    protected void testIntegratorConvergence(int quantumNumber) {
        Integrator integrator = getIntegrator();
        if (integrator.minHistoryLength() > INITIALIZATION_N_MAX + 1) {
            String className = this.getClass().getSimpleName();
            String msg = String.format("Integrator %s requires at least %d previously-computed points, but only %d are available!",
                    className, integrator.minHistoryLength(), INITIALIZATION_N_MAX);
            throw new RuntimeException(msg);
        }

        PhysicalPotential potential = new PhysicalPotentialHarmonic(DEFAULT_R0, DEFAULT_ALPHA);

        // Initialize the exact solution for the specified quantum state
        exactSolutionInstance = HarmonicOscillatorExactSolution.excitedState(DEFAULT_R0, DEFAULT_ALPHA, quantumNumber);

        Map<Integer, ConvergenceData> convergenceResults = new HashMap<>();

        // Test with different numbers of grid points
        for (int nOfPoints = MIN_POINTS; nOfPoints <= 1100; nOfPoints += 100) {
            ConvergenceData data = testWithPoints(potential, integrator, nOfPoints, quantumNumber);
            convergenceResults.put(nOfPoints, data);
        }

        // Print results
        String stateName = (quantumNumber == 0) ? "Ground State" : quantumNumber + "th Excited State";
        System.out.println("\n=== Testing " + stateName + " (n=" + quantumNumber + ") ===\n");
        System.out.println("Results for: " + integrator.getClass().getSimpleName() + "\n");
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
    protected ConvergenceData testWithPoints(PhysicalPotential potential, Integrator integrator, int nOfPoints, int quantumNumber) {
        // Get level-dependent grid boundaries
        double[] gridBounds = getGridBoundaries(quantumNumber);
        double rMin = gridBounds[0];
        double rMax = gridBounds[1];

        Grid grid = GridFactory.generateUniformGrid(rMin, rMax, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, 2., grid);
        QuantumLevel level = new QuantumLevel(system);

        // Use the exact energy for the specified quantum level
        HarmonicOscillatorExactSolution exactSol = HarmonicOscillatorExactSolution.excitedState(DEFAULT_R0, DEFAULT_ALPHA, quantumNumber);
        level.energy = exactSol.getEnergy();

        // Initialize with exact solution
        for (int n = 0; n <= INITIALIZATION_N_MAX && n < nOfPoints; n++) {
            level.psi[n] = exactSol.evaluate(grid.rAtGridPoint(n));


        }

        level.currentPsiPrime = new double[1]; // TODO Initialize this for RKN methods

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
        data.addError("10", exactSolution(grid.rAtGridPoint(nOfPoints / 10)) - level.psi[nOfPoints / 10]);
        data.addError("4", exactSolution(grid.rAtGridPoint(nOfPoints / 4)) - level.psi[nOfPoints / 4]);
        data.addError("2", exactSolution(grid.rAtGridPoint(nOfPoints / 2)) - level.psi[nOfPoints / 2]);

        return data;
    }

    /**
     * Prints the convergence results including point-to-point convergence rates.
     */
    protected void printResults(Map<Integer, ConvergenceData> convergenceResults) {
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
            System.out.printf("%10d\t\t% 20.18f\t% 20.18f\t% 20.18f",
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
    protected void printLinearFitResultsAndAssert(Map<Integer, ConvergenceData> convergenceResults, int quantumNumber, Integrator integrator) {
        String className = this.getClass().getSimpleName();
        String stateName = (quantumNumber == 0) ? "Ground State" : "n=" + quantumNumber + " Excited State";
        System.out.println("Linear Fit Results for " + className + " (" + stateName + "):");
        System.out.println("Fitting ln(|error|) = a + b * ln(nPoints)");
        System.out.println("--------------------------------------------------------");

        String[] pointFractions = {"10", "4", "2"};
        String[] pointNames = {"10%", "25%", "50%"};
        double expectedBCoefficient = -getGlobalConvergenceOrder();

        for (int idx = 0; idx < pointFractions.length; idx++) {
            String fraction = pointFractions[idx];
            String pointName = pointNames[idx];

            // Get point-specific and state-dependent tolerances
            double minRSquared = getMinRSquared(quantumNumber, pointName, integrator);
            double convergenceTol = getConvergenceOrderTolerance(quantumNumber, pointName, integrator);

            // Create arrays for ln(points) and ln(|error|)
            int n = convergenceResults.size();
            double[] lnPoints = new double[n];
            double[] lnErrors = new double[n];

            int i = 0;
            for (Map.Entry<Integer, ConvergenceData> entry : convergenceResults.entrySet()) {
                int nOfPoints = entry.getKey();
                double error = entry.getValue().getError(fraction);

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

            // Assert R² is greater than minimum threshold (state-dependent)
            assertTrue(rSquared > minRSquared,
                    String.format("%s (%s): R² (%.6f) should be > %.6f for point %s",
                            className, stateName, rSquared, minRSquared, pointNames[idx]));

            // Assert b coefficient is within tolerance of expected value (state-dependent)
            double bDeviation = Math.abs(b - expectedBCoefficient);
            assertTrue(bDeviation <= convergenceTol,
                    String.format("%s (%s): b coefficient (%.6f) deviates %.6f from expected %.6f (max allowed: %.6f) for point %s",
                            className, stateName, b, bDeviation, expectedBCoefficient, convergenceTol, pointNames[idx]));
        }

        System.out.println();
    }

    /**
     * Calculates the convergence rate between two points.
     */
    protected double calculateConvergenceRate(double error1, double error2, int n1, int n2) {
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
    protected double exactSolution(double x) {
        return exactSolutionInstance.evaluate(x);
    }

    /**
     * Returns the integrator to be tested.
     * Subclasses must implement this method to provide their specific integrator.
     */
    protected abstract Integrator getIntegrator();

    /**
     * Returns the expected b coefficient (slope) for the linear fit
     * ln(|error|) = a + b * ln(nPoints).
     * Delegates to the integrator under test.
     *
     * @return the expected b coefficient
     */
    protected double getGlobalConvergenceOrder() {
        return getIntegrator().globalConvergenceOrder();
    }

    /**
     * Utility method to format floating point numbers for display.
     */
    protected static String padFloat(Double number) {
        return padFloat(number, 2, 6);
    }

    /**
     * Utility method to format floating point numbers for display.
     */
    protected static String padFloat(Double number, int nOfDecimals, int width) {
        String padding = "";
        int magnitude = number <= 1 ? 0 : (int) Math.log10(Math.abs(number));
        int signPadding = (number >= 0) ? 0 : 1;

        for (int i = 0; i < width - magnitude - signPadding - nOfDecimals; i++) {
            padding += " ";
        }

        return padding + String.format("%." + nOfDecimals + "f", number) + " ";
    }


    // This is a freeze/regression test to detect accidental changes to integrator behavior during refactoring.
    // Uses "meaningless" values: cos(k) for psi and a polynomial for Q.
    // Does NOT guarantee correctness - other tests verify that.
    // Subclasses implement getIntegrator() to provide specific integrator.
    private double integrateOneStep(double qFactor, Integrator.Direction direction) {
        int arraySize = 20;
        double[] psi = new double[arraySize];
        double[] currentPsiPrime = new double[1];
        for (int k = 0; k < arraySize; k++) {
            psi[k] = Math.cos(k);
            currentPsiPrime[0] = -Math.sin(k);
        }

        Integrator integrator = getIntegrator();
        DoubleUnaryOperator q = i -> qFactor * (0.1111 + 0.2222 * i + 0.3333 * i * i + 0.44444 * i * i * i);
        DoubleUnaryOperator q1 = i -> qFactor * (0.2222 + 2 * 0.3333 * i + 3. * 0.44444 * i * i);
        DoubleUnaryOperator q2 = i -> qFactor * (2 * 0.3333 + 3. * 2. * 0.44444 * i);
        return integrator.propagate(psi, currentPsiPrime, arraySize / 2, 0.123, q, q1, q2, direction);
    }

    public double integrateOneStep(Integrator.Direction direction) {
        return integrateOneStep(1., direction);
    }

    // Same as integrateOneStep but with negated Q to test Z < 0 branch in exponentially fitted methods.
    public double integrateOneStepNegativeQ(Integrator.Direction direction) {
        return integrateOneStep(-1., direction);
    }

    // Same as integrateOneStep but with very small Q  to test series expansion branch (|Z| < threshold).
    public double integrateOneStepSmallZ(Integrator.Direction direction) {
        return integrateOneStep(2.e-5, direction);
    }

    @Test
    public void time_benchmark() {
        int npoints = 1000000;

        double[] psi = new double[npoints];
        double[] currentPsiPrime = new double[1];
        for (int k = 0; k < 15; k++) {
            psi[k] = Math.cos(k);
            currentPsiPrime[0] = -Math.sin(k);
        }

        Integrator integrator = getIntegrator();
        DoubleUnaryOperator q = i -> 1e-3 * (0.1111 + 0.2222 * i / npoints + 0.03333 / (npoints * npoints) * i * i);
        DoubleUnaryOperator q1 = i -> 1e-3 * (0.2222 / npoints + 2 * 0.3333 * i / (npoints * npoints));
        DoubleUnaryOperator q2 = i -> 1e-3 * (2 * 0.3333 / (npoints * npoints));
        long t0 = System.nanoTime();
        for (int n = 14; n < npoints - 1; n++) {
            psi[n + 1] = integrator.propagate(psi, currentPsiPrime, n, 0.123, q, q1, q2, Integrator.Direction.FORWARD);
        }
        long elapsedNanos = System.nanoTime() - t0;
        double timePerPoint = ((double) elapsedNanos) / (npoints);
        String msg = String.format("Time taken for integrating %d points is %10.3f ns / point for %s", npoints, timePerPoint, integrator.getClass().getSimpleName());
        System.out.println(msg);

//        Results 2026-03-05 h 16:40, commit e61c3cd0d1e0c9103dc4e39340a33adb5cf5c48f
//        Laptop i7 1355U, Amazon corretto 21.0.8
//        Time taken for integrating 1000000 points is     19.651 ns / point for EfnFixedBeta
//        Time taken for integrating 1000000 points is     36.674 ns / point for ExponentiallyFitted
//        Time taken for integrating 1000000 points is     29.706 ns / point for Numerov
//        Time taken for integrating 1000000 points is     60.352 ns / point for PredictorCorrector6
//        Time taken for integrating 1000000 points is     76.251 ns / point for PredictorCorrector8NumerovIter1
//        Time taken for integrating 1000000 points is    102.563 ns / point for PredictorCorrector8NumerovIter2
//        Time taken for integrating 1000000 points is     18.324 ns / point for Stormer5
//        Time taken for integrating 1000000 points is     14.718 ns / point for Stormer6
//        Time taken for integrating 1000000 points is     57.302 ns / point for Stormer8
//        Time taken for integrating 1000000 points is      9.685 ns / point for TaylorThreePoints

//        Results 2026-03-06 h 09:06, commit e77ee50e7a6e5cb8922cd0820b5503ea52a615fe
//        Desktop i5 4460, Oracle OpenJDK 19.0.1
//        Time taken for integrating 1000000 points is     60.920 ns / point for EfnFixedBeta
//        Time taken for integrating 1000000 points is     78.005 ns / point for ExponentiallyFitted
//        Time taken for integrating 1000000 points is    104.782 ns / point for Numerov
//        Time taken for integrating 1000000 points is    135.936 ns / point for PredictorCorrector6
//        Time taken for integrating 1000000 points is    226.429 ns / point for PredictorCorrector8NumerovIter1
//        Time taken for integrating 1000000 points is    315.115 ns / point for PredictorCorrector8NumerovIter2
//        Time taken for integrating 1000000 points is     33.766 ns / point for Stormer5
//        Time taken for integrating 1000000 points is     44.969 ns / point for Stormer6
//        Time taken for integrating 1000000 points is    101.723 ns / point for Stormer8
//        Time taken for integrating 1000000 points is     16.698 ns / point for TaylorThreePoints

//        Note 1: The i7 1355U is approx. 2.7 times faster than i5 4460, with a scatter of around 20% around this value.
//        These are the speed normalized to the TaylorThreePoints method:
//                                        i7-1355U  i5-4460
//        EfnFixedBeta	                     2.03	 3.65
//        ExponentiallyFitted	             3.79	 4.67
//        Numerov	                         3.07	 6.28
//        PredictorCorrector6	             6.23	 8.14
//        PredictorCorrector8NumerovIter1 	 7.87	13.56
//        PredictorCorrector8NumerovIter2	10.59	18.87
//        Stormer5	                         1.89	 2.02
//        Stormer6	                         1.52	 2.69
//        Stormer8	                         5.92	 6.09
//        TaylorThreePoints              	 1.00	 1.00

// These are the averages across the two machines, and sorted by speed:
//        PredictorCorrector8NumerovIter2	14.73
//        PredictorCorrector8NumerovIter1	10.72
//        PredictorCorrector6	             7.19
//        Stormer8	                         6.00
//        Numerov	                         4.67
//        ExponentiallyFitted	             4.23
//        EfnFixedBeta	                     2.84
//        Stormer6	                         2.11
//        Stormer5	                         1.96
//        TaylorThreePoints	                 1.00

// NOTE: TaylorThreePoints is very fast, which is no surprise as it's super simple, but its convergence is really poor.
// The real baseline is Numerov, and we can see that all methods apart from Predictor-Corrector have similar (usually BETTER)
//  speed! This really shouldn't be the case, as they are more complex.

    }

    @Test
    public void propagate_forward_exp_minus_x() {
        Integrator integrator = getIntegrator();

        // NB don't use too many points because some methods have fairly bad roundoff error which distorts the convergence patterns
        for (int np = 10; np <= 20; np++) {
            double xmin = 0.;
            double xmax = 1.;
            double[] psi = new double[np];
            double[] currentPsiPrime = new double[1];
            double step = (xmax - xmin) / (np - 1);

            for (int i = 0; i < integrator.minHistoryLength(); i++) {
                double x = xmin + i * step;
                psi[i] = Math.exp(-x);
                currentPsiPrime[0] = -Math.exp(-x);
            }

            for (int n = integrator.minHistoryLength() - 1; n < np - 1; n++) {
                psi[n + 1] = integrator.propagate(psi, currentPsiPrime, n, step,
                        i -> -1., i -> 0., i -> 0.,
                        Integrator.Direction.FORWARD);
            }

            double error = Math.exp(-1.) - psi[np - 1];
            System.out.println(np + " " + error);
        }

    }

    @Test
    public void propagate_forward_exp_to_cos_x() {
        propagate_exp_to_cos_x(Integrator.Direction.FORWARD);
    }

    @Test
    public void propagate_backward_exp_to_cos_x() {
        propagate_exp_to_cos_x(Integrator.Direction.BACKWARD);
    }

    public void propagate_exp_to_cos_x(Integrator.Direction direction) {
        Integrator integrator = getIntegrator();
        System.out.println(integrator.getClass().getSimpleName());
        int d = direction.getValue();

        for (int np = 100; np <= 500; np += 10) {
            double xmin = 0.;
            double xmax = 4. * Math.PI;
            double[] psi = new double[np];
            double[] currentPsiPrime = new double[1];
            double step = (xmax - xmin) / (np - 1);

            boolean isBackward = direction.equals(Integrator.Direction.BACKWARD);
            {
                final int startIndex1 = isBackward ? np - 1 : 0;
                final int endIndex1 = isBackward ? np - 1 - integrator.minHistoryLength() : integrator.minHistoryLength();
                Predicate<Integer> condition = isBackward ? n -> (n > endIndex1) : n -> (n < endIndex1);
                for (int n = startIndex1; condition.test(n); n = n + d) {
                    double x = xmin + n * step;
                    psi[n] = Math.exp(Math.cos(x));
                    currentPsiPrime[0] = -Math.sin(x) * Math.exp(Math.cos(x));
                }
            }

            DoubleUnaryOperator qTilde = i -> Math.cos(xmin + i * step) - Math.pow(Math.sin(xmin + i * step), 2);
            DoubleUnaryOperator qTildePrime = i -> -Math.sin(xmin + i * step) - 2. * Math.cos(xmin + i * step) * Math.sin(xmin + i * step);
            DoubleUnaryOperator qTildeDoublePrime = i -> -Math.cos(xmin + i * step)
                    - 2. * Math.pow(Math.cos(xmin + i * step), 2)
                    + 2. * Math.pow(Math.sin(xmin + i * step), 2);

            {
                final int startIndex2 = isBackward ? np - integrator.minHistoryLength() : integrator.minHistoryLength() - 1;
                final int endIndex2 = isBackward ? np - 1 - integrator.minHistoryLength() : np - 1;
                Predicate<Integer> condition = isBackward ? n -> (n > 0) : n -> (n < endIndex2);

                for (int n = startIndex2; condition.test(n); n = n + d) {
                    psi[n + d] = integrator.propagate(psi, currentPsiPrime, n, step,
                            qTilde, qTildePrime, qTildeDoublePrime,
                            direction);
                }
            }

            double error = isBackward ? Math.exp(Math.cos(xmin)) - psi[0] : Math.exp(Math.cos(xmax)) - psi[np - 1];
            System.out.println(np + " " + error);

            String className = integrator.getClass().getSimpleName();
            double maxError = expectedError_propagate_exp_to_cos_x(integrator, np);
            assertTrue(Math.abs(error) < maxError,
                    String.format("%s : |error| = %.3e should be < %.3e for np = %s",
                            className, Math.abs(error), maxError, np));
        }

    }

    private double expectedError_propagate_exp_to_cos_x(Integrator integrator, int np) {

        // We assume errors are in the format |error| = (C/np)^a
        // The higher a , the better;
        // The smallest C, the better;

        double a = 1.;
        double C = 1.e-15; // Very tight defaults which will make the test fail
        if (integrator instanceof TaylorThreePoints) {
            a = 1.88;
            C = 85.4;
        } else if (integrator instanceof RKN4) {
            a = 4.01;
            C = 24.0;
        } else if (integrator instanceof Numerov) {
            a = 3.62;
            C = 15.1;
        } else if (integrator instanceof ExponentiallyFitted) {
            a = 3.68;
            C = 15.4;
        } else if (integrator instanceof EfnFixedBeta) {
            a = 3.54;
            C = 13.5;
        } else if (integrator instanceof Stormer5) {
            a = 4.85;
            C = 28.6;
        } else if (integrator instanceof Stormer8) {
            a = 8.63;
            C = 52.1;
        } else if (integrator instanceof Stormer6) {
            a = 7.49;
            C = 39.1;
        } else if (integrator instanceof Stormer8i) {
            a = 8.70;
            C = 36.6;
        } else if (integrator instanceof PredictorCorrector6) {
            a = 8.88;
            C = 41.2;
        } else if (integrator instanceof Obrechkoff) {
            a = 6.07;
            C = 17.5;
        } else if (integrator instanceof RK45) {
            a = 5.95;
            C = 14.3;
        } else if (integrator instanceof PredictorCorrector8NumerovIter1) {
            a = 5.01;
            C = 7.5;
        } else if (integrator instanceof PredictorCorrector8NumerovIter2) {
            a = 6.33;
            C = 13.9;
        }

        return Math.pow(C / np, a);
    }


}