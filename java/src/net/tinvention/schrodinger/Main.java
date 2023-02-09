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
		System.out.println();
		for (int i = 6; i < 7; i++) {
			int nOfPoints = 20 + (int) Math.pow(i, 3);
//			double rRef = 4;
//			double alpha = 0.5;
//			Grid grid = Grid.generateSurkusGrid(0.5d, 5.5d, nOfPoints, rRef, alpha);
//			Grid grid = Grid.generateSqrtGrid(0.5d, 5.5d, nOfPoints, rRef);
			Grid grid = Grid.generateUniformGrid(10.0d, 20.0d, nOfPoints);
//			System.out.println(grid.toString());
			double mass = 2.0d;
			EigenvalueFinder finder = new EigenvalueFinder(grid, potential, mass, integrator);

			int nOfDesiredNodes = 0;
			EnergyLevel ek = finder.findEigenvalue(nOfDesiredNodes);
			System.out.println(i + " " + grid.getNumberOfPoints() + " " + (ek.energy - 0.5) + " " + " "
					+ ek.perturbativeCorrectionToEnergy + " " + ek.numberOfBisections);
//			for (int j = 0; j < grid.getNumberOfPoints(); j++) {
//				System.out.println(grid.getRValue(j) + " " + grid.getRValue(j) + " " + ek.psi[j]);
//			}
		}

	}
}
