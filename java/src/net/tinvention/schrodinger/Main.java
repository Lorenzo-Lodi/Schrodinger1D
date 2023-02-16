package net.tinvention.schrodinger;

import net.tinvention.schrodinger.grid.Grid;
import net.tinvention.schrodinger.integrator.Integrator;
import net.tinvention.schrodinger.integrator.TaylorThreePoints;
import net.tinvention.schrodinger.potential.HarmonicPotential;
import net.tinvention.schrodinger.potential.BarePotential;

public class Main {

	public static void main(String[] args) {

		BarePotential potential = new HarmonicPotential(15, 1);
		Integrator integrator = new TaylorThreePoints();

		for (int i = 1; i < 17; i++) {
			int nOfPoints = 20 + (int) Math.pow(i, 3);
//			double rRef = 5;
//			double alpha = 0.5;
//			Grid grid = Grid.generateSurkusGrid(0.5d, 5.5d, nOfPoints, rRef, alpha);
//			Grid grid = Grid.generateSqrtGrid(10d, 20.0d, nOfPoints, rRef);
			Grid grid = Grid.generateUniformGrid(5.0d, 25.0d, nOfPoints);
//			System.out.println(grid.toString());
			double mass = 2.0d;
			EigenvalueFinder finder = new EigenvalueFinder(grid, potential, mass, integrator);

			int nOfDesiredNodes = 0;
			EnergyLevel ek = finder.findEigenvalueByBisection(nOfDesiredNodes);

			int floatDecimals = 16;
			System.out.println(padInt(i, 4) + padInt(grid.getNumberOfPoints(), 6)
					+ padFloat(ek.lowerBound, floatDecimals, floatDecimals + 4)
					+ padFloat(ek.upperBound, floatDecimals, floatDecimals + 4)
					+ padFloat(ek.energy, floatDecimals, floatDecimals + 4)
					+ padFloat(ek.energy +ek.perturbativeCorrectionToEnergy, floatDecimals, floatDecimals + 7)
					+ padInt(ek.numberOfBisections, 6));
			for (int j = 0; j < grid.getNumberOfPoints(); j++) {
				// System.out.println(grid.getYValue(j) + " " + grid.getRValue(j) + " " +
				// ek.psi[j]);
			}
		}

	}

	private static String padFloat(Double number, int nOfDecimals, int width) {
		return padding(number, width - nOfDecimals - 2) + String.format("%." + nOfDecimals + "f", number);
	}

	private static String padInt(Integer number, int width) {
		return padding(number, width) + number.toString();
	}

	private static String padding(double number, int width) {
		int magnitude = number == 0 ? 0 : (int) Math.log10(Math.abs(number));
		String padding = "";
		int signPadding = (number >= 0) ? 0 : 1;
		for (int i = 0; i < width - magnitude - signPadding; i++) {
			padding += " ";
		}
		return padding;
	}
}
