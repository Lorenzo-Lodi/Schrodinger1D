package net.tinvention.schrodinger;

import net.tinvention.schrodinger.grid.Grid;
import net.tinvention.schrodinger.grid.SurkusGrid;
import net.tinvention.schrodinger.grid.UniformGrid;
import net.tinvention.schrodinger.integrator.Integrator;
import net.tinvention.schrodinger.integrator.TaylorThreePoints;
import net.tinvention.schrodinger.potential.HarmonicPotential;
import net.tinvention.schrodinger.potential.ClampedNucleiPotential;

public class Main {

	public static void main(String[] args) {

		Grid grid = new SurkusGrid(0.1, 5.5d, 15, 1, 1);
		System.out.println(grid.printGrid());

		ClampedNucleiPotential v = new HarmonicPotential();
		Integrator integrator = new TaylorThreePoints();

		for (int i = 0; i < 16; i++) {
			int nOfPoints = 20 + (int) Math.pow(i, 3);
			grid = new UniformGrid(-5.5d, 5.5d, nOfPoints);
			double mass = 2.0d;
			EigenvalueFinder finder = new EigenvalueFinder(grid, v, mass, integrator);

			int nOfDesiredNodes = 1;
			EnergyLevel ek = finder.findEigenvalue(nOfDesiredNodes);
			// System.out.println(i + " " + grid.getStepSize() + " " +
			// grid.getNumberOfPoints() + " " + ek.energy);
		}

	}
}
