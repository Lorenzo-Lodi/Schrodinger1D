package net.tinvention.schrodinger;

import net.tinvention.schrodinger.grid.UniformGrid;
import net.tinvention.schrodinger.potential.HarmonicPotential;
import net.tinvention.schrodinger.potential.Potential;

public class Main {

	public static void main(String[] args) {

		Potential v = new HarmonicPotential();

		for (int i = 0; i < 16; i++) {
			double newh = Math.pow(10.d, -1.d - ((double) i) * 0.25d);
			UniformGrid grid = new UniformGrid(-5.5d, 5.5d, newh);
			EigenvalueFinder finder = new EigenvalueFinder(grid, v);

			for (int nOfDesiredNodes = 0; nOfDesiredNodes < 1; nOfDesiredNodes++) {
				double ek = finder.findEigenvalue(nOfDesiredNodes, grid, v);
				System.out.println(i + " " + grid.h + " " + ek);
			}
		}
	}

}
