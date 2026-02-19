package schrodinger.integrator;

import org.junit.jupiter.api.Test;
import schrodinger.QuantumState;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.SchrodingerSystem;


public class IntegratorsTest {

    @Test
    public void propagate_test() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(10, 1.0);
        int nOfPoints = 100;
        Grid grid = GridFactory.generateUniformGrid(5.0, 15.0, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, 2., grid);
        QuantumState state = new QuantumState(system);
        state.energy = 0.5;

        int n = nOfPoints / 10;
        for (int k = 0; k <= n; k++) {
            state.psi[k] = exactSolution(grid.getRValue(k));
        }

        double exact = exactSolution(grid.getRValue(n + 1));
        double f1 = new TaylorThreePoints().propagate(state.psi, n, state, Integrator.Direction.FORWARD);
        double f2 = new Numerov().propagate(state.psi, n, state, Integrator.Direction.FORWARD);
        double f3 = new Vignoli().propagate(state.psi, n, state, Integrator.Direction.FORWARD);
        double f4 = new RaptisAllison().propagate(state.psi, n, state, Integrator.Direction.FORWARD);
        double f5 = new RaptisAllison2().propagate(state.psi, n, state, Integrator.Direction.FORWARD);
        double f6 = new ExponentiallyFitted().propagate(state.psi, n, state, Integrator.Direction.FORWARD);
        double f7 = new Stormer7().propagate(state.psi, n, state, Integrator.Direction.FORWARD);
        double f8 = new Stormer8().propagate(state.psi, n, state, Integrator.Direction.FORWARD);
        double f9 = new PredictorCorrector6().propagate(state.psi, n, state, Integrator.Direction.FORWARD);
        double f10 = new PredictorCorrector8().propagate(state.psi, n, state, Integrator.Direction.FORWARD);

        myPrint("TaylorThreePoints", f1, exact);
        myPrint("Numerov", f2, exact);
        myPrint("Vignoli", f3, exact);
        myPrint("RaptisAllison", f4, exact);
        myPrint("RaptisAllison2", f5, exact);
        myPrint("ExponentiallyFitted", f6, exact);
        myPrint("Stormer7", f7, exact);
        myPrint("Stormer8", f8, exact);
        myPrint("PredictorCorrector6", f9, exact);
        myPrint("PredictorCorrector8", f10, exact);

    }

    private static void myPrint(String s, double f1, double f2) {
        System.out.println(padFloat(f1) + "  " + padFloat(f2 - f1) + " <- " + s);
    }

    private static double exactSolution(double x) {
        return Math.exp(-(x - 10.) * (x - 10.));
    }

    private static String padFloat(double number) {
        return padFloat(number, 15, 18);
    }

    private static String padInt(int number) {
        return padInt(number, 6);
    }

    private static String padding(double number, int width) {
        int magnitude = number <= 1 ? 0 : (int) Math.log10(Math.abs(number));
        String padding = "";
        int signPadding = (number >= 0) ? 0 : 1;
        for (int i = 0; i < width - magnitude - signPadding; i++) {
            padding += " ";
        }
        return padding;
    }

    private static String padFloat(Double number, int nOfDecimals, int width) {
        return padding(number, width - nOfDecimals - 2) + String.format("%." + nOfDecimals + "f", number);
    }

    private static String padInt(Integer number, int width) {
        return padding(number, width) + number;
    }

}
