package schrodinger;

import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.TaylorThreePoints;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.SchrodingerSystem;
import schrodinger.solver.ShootingSolver;

public class Main {

    public static void main(String[] args) {

        PhysicalPotential potential = new PhysicalPotentialHarmonic(0, 1);
        Integrator integrator = new TaylorThreePoints();
        double mass = 2.0d;
        int nOfDesiredNodes = 1;

        System.out.println("   i nPoints       1/Ystep          Lower               Upper" + "               Energy"
                + "          Energy + Pert. " + "     n bisec");
        for (int i = 0; i < 20; i++) {
            int nOfPoints = 20 + 50 * i;
//            int nOfPoints = 81;
            Grid grid = GridFactory.generateUniformGrid(-4.0d, 4.0d, nOfPoints);
            SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
            ShootingSolver finder = new ShootingSolver(system, integrator);
            EnergyLevel ek = finder.findEigenvalueByBisection(nOfDesiredNodes);

            int floatDecimals = 16;
            System.out.println(padInt(i, 4) + padInt(grid.getNumberOfPoints(), 6)
                    + padFloat(1.0 / grid.getStepSizeYCoordinate(), 3, 14)
                    + padFloat(ek.lowerBound, floatDecimals, floatDecimals + 4)
                    + padFloat(ek.upperBound, floatDecimals, floatDecimals + 4)
                    + padFloat(ek.energy, floatDecimals, floatDecimals + 4)
                    + padFloat(ek.energy + ek.perturbativeCorrectionToEnergy, floatDecimals, floatDecimals + 6)
                    + padInt(ek.numberOfBisections, 6));
        }

    }

    private static String padFloat(Double number, int nOfDecimals, int width) {
        return padding(number, width - nOfDecimals - 2) + String.format("%." + nOfDecimals + "f", number);
    }

    private static String padInt(Integer number, int width) {
        return padding(number, width) + number.toString();
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
}
