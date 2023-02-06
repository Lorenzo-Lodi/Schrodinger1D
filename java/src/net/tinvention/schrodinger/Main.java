package net.tinvention.schrodinger;

import net.tinvention.schrodinger.grid.Grid;
import net.tinvention.schrodinger.integrator.Integrator;
import net.tinvention.schrodinger.integrator.TaylorThreePoints;
import net.tinvention.schrodinger.potential.HarmonicPotential;
import net.tinvention.schrodinger.potential.ClampedNucleiPotential;

public class Main {

	public static void main(String[] args) {

		ClampedNucleiPotential potential = new HarmonicPotential(3, 1);
		Integrator integrator = new TaylorThreePoints();
		System.out.println();
		for (int i = 0; i < 11; i++) {
//			for (int i = 0; i < 20; i++) {
			int nOfPoints = 20 + (int) Math.pow(i, 3);
			double rRef = 30;
			double alpha = 0.5;
			Grid grid = Grid.generateSurkusGrid(0.5d, 5.5d, nOfPoints, rRef, alpha);
//			System.out.println(grid.toString());
			double mass = 2.0d;
			EigenvalueFinder finder = new EigenvalueFinder(grid, potential, mass, integrator);

			int nOfDesiredNodes = 0;
			EnergyLevel ek = finder.findEigenvalue(nOfDesiredNodes);
			System.out.println(i + " " + grid.getNumberOfPoints() + " " + ek.energy);
		}

	}
}
