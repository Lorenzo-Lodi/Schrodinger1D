package schrodinger.integrator;

import org.junit.jupiter.api.Test;
import schrodinger.QuantumState;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.SchrodingerSystem;

import java.util.HashMap;
import java.util.Map;

/**
 * Abstract base class for integrator tests.
 * Provides common functionality for testing integrator convergence.
 * Supports testing against exact solutions for any quantum state of the harmonic oscillator.
 */
public abstract class AbstractIntegratorTest {

    private final static int INITIALIZATION_N_MAX = 4;
    private final static double MIN_R_SQUARED = 0.995;
    private final static double CONVERGENCE_ORDER_TOL = 0.3;
    private final static int MIN_POINTS = 200;
    
    // Default potential parameters (matching the existing test setup)
    protected static final double DEFAULT_R0 = 10.0;
    protected static final double DEFAULT_ALPHA = 1.0;
    
    // The exact solution for the current quantum state
    private HarmonicOscillatorExactSolution exactSolutionInstance;

    /**
     * Represents the convergence data at different points in the grid.
     * Uses a map to store error values for different fractions of the grid.
     */
    protected static class ConvergenceData {
        private final Map<String, Double> errors = new HashMap<>();

        public void addError(String pointFraction, double error) {
            errors.put(pointFraction, error);
        }

        public double getError(String pointFraction) {
            return errors.getOrDefault(pointFraction, 0.0);
        }

        public Map<String, Double> getAllErrors() {
            return new HashMap<>(errors);
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
        // Initialize the exact solution for the specified quantum state
        exactSolutionInstance = HarmonicOscillatorExactSolution.excitedState(DEFAULT_R0, DEFAULT_ALPHA, quantumNumber);
        
        PhysicalPotential potential = new PhysicalPotentialHarmonic(DEFAULT_R0, DEFAULT_ALPHA);
        Integrator integrator = getIntegrator();

        Map<Integer, ConvergenceData> convergenceResults = new HashMap<>();

        // Test with different numbers of grid points
        for (int nOfPoints = MIN_POINTS; nOfPoints <= 1100; nOfPoints += 100) {
            ConvergenceData data = testWithPoints(potential, integrator, nOfPoints, quantumNumber);
            convergenceResults.put(nOfPoints, data);
        }

        // Print results
        String stateName = (quantumNumber == 0) ? "Ground State" : quantumNumber + "th Excited State";
        System.out.println("\n=== Testing " + stateName + " (n=" + quantumNumber + ") ===\n");
        printResults(convergenceResults);

        // Calculate and print convergence rates
        printConvergenceRates(convergenceResults);
    }

    /**
     * Returns the exact energy for the current quantum state.
     * 
     * @return the exact energy eigenvalue
     */
    protected double getExactEnergy() {
        return exactSolutionInstance.getEnergy();
    }

    /**
     * Tests the integrator with a specific number of grid points.
     *
     * @param potential  The potential to use
     * @param integrator The integrator to test
     * @param nOfPoints  The number of grid points
     * @return Convergence data containing errors at different points
     */
    protected ConvergenceData testWithPoints(PhysicalPotential potential, Integrator integrator, int nOfPoints) {
        return testWithPoints(potential, integrator, nOfPoints, 0);
    }

    /**
     * Tests the integrator with a specific number of grid points for a given quantum state.
     *
     * @param potential     The potential to use
     * @param integrator    The integrator to test
     * @param nOfPoints    The number of grid points
     * @param quantumNumber The quantum number for the state being tested
     * @return Convergence data containing errors at different points
     */
    protected ConvergenceData testWithPoints(PhysicalPotential potential, Integrator integrator, int nOfPoints, int quantumNumber) {
        Grid grid = GridFactory.generateUniformGrid(6.0, 14.0, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, 2., grid);
        QuantumState state = new QuantumState(system);
        
        // Use the exact energy for the specified quantum state
        HarmonicOscillatorExactSolution exactSol = HarmonicOscillatorExactSolution.excitedState(DEFAULT_R0, DEFAULT_ALPHA, quantumNumber);
        state.energy = exactSol.getEnergy();

        // Initialize with exact solution - use more points for higher quantum numbers
        int initMax = Math.max(INITIALIZATION_N_MAX, quantumNumber + 1);
        for (int n = 0; n <= initMax && n < nOfPoints; n++) {
            state.psi[n] = exactSol.evaluate(grid.getRValue(n));
        }

        // Propagate the wavefunction
        double step = state.getGrid().getStepSizeYCoordinate();
        for (int n = initMax; n < nOfPoints - 1; n++) {
            state.psi[n + 1] = integrator.propagate(state.psi, n, step, state::QTildeValueAt, Integrator.Direction.FORWARD);
        }

        // Calculate errors at different points
        ConvergenceData data = new ConvergenceData();
        data.addError("10", exactSolution(grid.getRValue(nOfPoints / 10)) - state.psi[nOfPoints / 10]);
        data.addError("4", exactSolution(grid.getRValue(nOfPoints / 4)) - state.psi[nOfPoints / 4]);
        data.addError("2", exactSolution(grid.getRValue(nOfPoints / 2)) - state.psi[nOfPoints / 2]);

        return data;
    }

    /**
     * Initializes the quantum state with the exact solution.
     * Subclasses can override this method to handle different initialization requirements.
     */
    protected void initializeState(QuantumState state, Grid grid) {
        for (int n = 0; n <= INITIALIZATION_N_MAX; n++) {
            state.psi[n] = exactSolution(grid.getRValue(n));
        }
    }

    /**
     * Prints the convergence results.
     */
    protected void printResults(Map<Integer, ConvergenceData> convergenceResults) {
        System.out.println("Grid Points\t\tError at 10%\t\t\tError at 25%\t\t\tError at 50%");
        System.out.println("------------------------------------------------------------");

        // Sort the entries by number of points
        convergenceResults.entrySet().stream()
                .sorted(Map.Entry.comparingByKey())
                .forEach(entry -> {
                    int nOfPoints = entry.getKey();
                    ConvergenceData data = entry.getValue();

                    System.out.printf("%10d\t\t%20.16f\t%20.16f\t%20.16f%n",
                            nOfPoints,
                            data.getError("10"),
                            data.getError("4"),
                            data.getError("2"));
                });

        System.out.println();
    }

    /**
     * Calculates and prints the convergence rates.
     */
    protected void printConvergenceRates(Map<Integer, ConvergenceData> convergenceResults) {
        System.out.println("Convergence Rates:");
        System.out.println("Point 10%\tPoint 25%\tPoint 50%");
        System.out.println("------------------------------------");

        for (int i = 1; i < convergenceResults.size(); i++) {
            int n1 = MIN_POINTS + (i - 1) * 100;
            int n2 = n1 + 100;

            if (convergenceResults.containsKey(n1) && convergenceResults.containsKey(n2)) {
                ConvergenceData data1 = convergenceResults.get(n1);
                ConvergenceData data2 = convergenceResults.get(n2);

                double rate1 = calculateConvergenceRate(data1.getError("10"), data2.getError("10"), n1, n2);
                double rate2 = calculateConvergenceRate(data1.getError("4"), data2.getError("4"), n1, n2);
                double rate3 = calculateConvergenceRate(data1.getError("2"), data2.getError("2"), n1, n2);

                System.out.printf("%.4f\t\t%.4f\t\t%.4f%n", rate1, rate2, rate3);
            }
        }

        System.out.println();

        // Print linear fit results
        printLinearFitResults(convergenceResults);
    }

    /**
     * Computes and prints linear fit parameters for convergence analysis.
     * Fits ln(|error|) = a + b * ln(nPoints) using least squares regression.
     */
    protected void printLinearFitResults(Map<Integer, ConvergenceData> convergenceResults) {
        String className = this.getClass().getSimpleName();
        System.out.println("Linear Fit Results for " + className + ":");
        System.out.println("Fitting ln(|error|) = a + b * ln(nPoints)");
        System.out.println("--------------------------------------------------------");

        String[] pointFractions = {"10", "4", "2"};
        String[] pointNames = {"10%", "25%", "50%"};
        double expectedBCoefficient = -getGlobalConvergenceOrder();

        for (int idx = 0; idx < pointFractions.length; idx++) {
            String fraction = pointFractions[idx];

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

            System.out.printf("Point %s:\t\ta = %.6f\tb = %.6f\tR² = %.6f%n",
                    pointNames[idx], a, b, rSquared);

            // Assert R² is greater than 0.95
            org.junit.jupiter.api.Assertions.assertTrue(rSquared > MIN_R_SQUARED,
                    String.format("%s: R² (%.6f) should be > (%.6f) for point %s", className, rSquared, MIN_R_SQUARED, pointNames[idx]));

            // Assert b coefficient is within ±0.3 of expected value
            double bDeviation = Math.abs(b - expectedBCoefficient);
            org.junit.jupiter.api.Assertions.assertTrue(bDeviation <= CONVERGENCE_ORDER_TOL,
                    String.format("%s: b coefficient (%.6f) deviates %.6f from expected %.6f (max allowed: %.6f) for point %s",
                            className, b, bDeviation, expectedBCoefficient, CONVERGENCE_ORDER_TOL, pointNames[idx]));
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
     * Subclasses must implement this method to provide the theoretical convergence rate.
     *
     * @return the expected b coefficient
     */
    protected abstract double getGlobalConvergenceOrder();

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
}