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

		ClampedNucleiPotential potential = new HarmonicPotential(3, 1);
		Integrator integrator = new TaylorThreePoints();
//		{
//			Grid grid = new SurkusGrid(0.5d, 5.5d, 20, 3, 1);
//			for (int i = 0; i < 100; i++) {
//				double r = 0.01 + ((double) i) / 2.;
//				System.out.println(r + " " + grid.mappingFunctionYofR(r) + " "
//						+ grid.mappingFunctionRofY(grid.mappingFunctionYofR(r)));
//			}
//		}
//		{
//			Grid grid = new SurkusGrid(0.5d, 5.5d, 30, 3, 1);
//			double mass = 2.0d;
//		}
		System.out.println();
		for (int i = 0; i < 15; i++) {
			int nOfPoints = 20 + (int) Math.pow(i, 3);
//			Grid grid = new UniformGrid(0.5d, 5.5d, nOfPoints);
			Grid grid = new SurkusGrid(0.5d, 5.5d, nOfPoints, 3, 4);
			double mass = 2.0d;
			EigenvalueFinder finder = new EigenvalueFinder(grid, potential, mass, integrator);

			int nOfDesiredNodes = 1;
			EnergyLevel ek = finder.findEigenvalue(nOfDesiredNodes);
			System.out.println(
					i + " " + grid.getStepSizeYCoordinate() + " " + grid.getNumberOfPoints() + " " + ek.energy);
		}

	}
}
