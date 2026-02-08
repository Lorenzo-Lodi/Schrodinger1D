package schrodinger;

import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.TaylorThreePoints;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.PhysicalPotential;

public class Main {

	public static void main(String[] args) {

		PhysicalPotential potential = new PhysicalPotentialHarmonic(15, 1);
		Integrator integrator = new TaylorThreePoints();

		System.out.println("   i nPoints       1/Ystep          Lower               Upper" + "               Energy"
				+ "          Energy + Pert. " + "     n bisec");
		for (int i = 0; i < 20; i++) {
//			int nOfPoints = 20 + (int) Math.pow(i, 3);
			double rmax = 16 + i;
			double rmin = 10 - i;
//			double rRef = 5;
//			double alpha = 0.5;
//			Grid grid = GridFactory.generateSurkusGrid(0.5d, 5.5d, nOfPoints, rRef, alpha);
//			Grid grid = GridFactory.generateSqrtGrid(10d, 20.0d, nOfPoints, rRef);
//			Grid grid = GridFactory.generateUniformGrid(10.0d, 20.0d, nOfPoints);
			int nOfPoints = (int) (rmax - rmin) * 100;
			Grid grid = GridFactory.generateUniformGrid(10.0d, rmax, nOfPoints);
//			System.out.println(grid.toString());
			double mass = 2.0d;
			EigenvalueFinder finder = new EigenvalueFinder(grid, potential, mass, integrator);

			int nOfDesiredNodes = 10;
			EnergyLevel ek = finder.findEigenvalueByBisection(nOfDesiredNodes);

			int floatDecimals = 16;
			System.out.println(padInt(i, 4) + padInt(grid.getNumberOfPoints(), 6)
					+ padFloat(1.0 / grid.getStepSizeYCoordinate(), 3, 14)
					+ padFloat(ek.lowerBound, floatDecimals, floatDecimals + 4)
					+ padFloat(ek.upperBound, floatDecimals, floatDecimals + 4)
					+ padFloat(ek.energy, floatDecimals, floatDecimals + 4)
					+ padFloat(ek.energy + ek.perturbativeCorrectionToEnergy, floatDecimals, floatDecimals + 6)
					+ padInt(ek.numberOfBisections, 6));
//			for (int j = 0; j < grid.getNumberOfPoints(); j++) {
//				 System.out.println(grid.getYValue(j) + " " + grid.getRValue(j) + " " +
//				 ek.psi[j]);
//			}
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
