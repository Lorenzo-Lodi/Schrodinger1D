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
 */
public abstract class AbstractIntegratorTest {

    private final static int INITIALIZATION_N_MAX = 4;

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
     * Tests the convergence of the integrator by comparing numerical solutions
     * with the exact analytical solution.
     */
    @Test
    public void testIntegratorConvergence() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(10, 1.0);
        Integrator integrator = getIntegrator();

        Map<Integer, ConvergenceData> convergenceResults = new HashMap<>();

        // Test with different numbers of grid points
        for (int nOfPoints = 100; nOfPoints <= 1100; nOfPoints += 100) {
            ConvergenceData data = testWithPoints(potential, integrator, nOfPoints);
            convergenceResults.put(nOfPoints, data);
        }

        // Print results
        printResults(convergenceResults);

        // Calculate and print convergence rates
        printConvergenceRates(convergenceResults);
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
        Grid grid = GridFactory.generateUniformGrid(6.0, 14.0, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, 2., grid);
        QuantumState state = new QuantumState(system);
        state.energy = 0.5;

        // Initialize with exact solution
        initializeState(state, grid);

        // Propagate the wavefunction
        double step = state.getGrid().getStepSizeYCoordinate();
        for (int n = INITIALIZATION_N_MAX; n < nOfPoints - 1; n++) {
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
            int n1 = 100 + (i - 1) * 100;
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
     * Returns the exact analytical solution at position x.
     */
    protected double exactSolution(double x) {
        return Math.exp(-(x - 10.) * (x - 10.));
    }

    /**
     * Returns the integrator to be tested.
     * Subclasses must implement this method to provide their specific integrator.
     */
    protected abstract Integrator getIntegrator();

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