package net.tinvention.schrodinger;

import net.tinvention.schrodinger.grid.Grid;
import net.tinvention.schrodinger.grid.QuadraticGrid;
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
		System.out.println();
		for (int i = 5; i < 6; i++) {
//			for (int i = 0; i < 20; i++) {
			int nOfPoints = 20 + (int) Math.pow(i, 3);

			double rRef = 3;
			double alpha = 1;
//			Grid grid = new SurkusGrid(0.5d, 5.5d, nOfPoints, rRef, alpha);
//			Grid grid = new UniformGrid(0.5d, 5.5d, nOfPoints);
			Grid grid = new QuadraticGrid(0.5d, 10.5d, nOfPoints, 3);
			System.out.println(grid.toString());
			double mass = 2.0d;
			EigenvalueFinder finder = new EigenvalueFinder(grid, potential, mass, integrator);

			int nOfDesiredNodes = 1;
			EnergyLevel ek = finder.findEigenvalue(nOfDesiredNodes);
			System.out.println(i + " " + grid.getNumberOfPoints() + " " + ek.energy);
		}

	}
}
