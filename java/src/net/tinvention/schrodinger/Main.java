package net.tinvention.schrodinger;

import net.tinvention.schrodinger.grid.UniformGrid;
import net.tinvention.schrodinger.potential.HarmonicPotential;
import net.tinvention.schrodinger.potential.Potential;

public class Main {
	private static double mass = 2.0d;

	public static void main(String[] args) {

		Potential v = new HarmonicPotential();

		for (int i = 0; i < 16; i++) {
			double newh = Math.pow(10.d, -1.d - ((double) i) * 0.25d);
			UniformGrid grid = new UniformGrid(-5.5d, 5.5d, newh);
			EigenvalueFinder finder = new EigenvalueFinder(grid, v);

			for (int nOfDesiredNodes = 0; nOfDesiredNodes < 1; nOfDesiredNodes++) {
				EnergyLevel trial = finder.computeApproximateEnergyLevel(nOfDesiredNodes);
				double ek = findEigenvalue(nOfDesiredNodes, trial.lowerBound, trial.upperBound, trial.energy, grid, v);
				System.out.println(i + " " + grid.h + " " + ek);
			}
		}
	}

	private static double findEigenvalue(int nOfDesiredNodes, double eLowInput, double eHighInput, double eTrial,
			UniformGrid grid, Potential v) {
		double target_relative_error = Math.ulp(1.d); // change for single-precision float
		double e = eTrial;
		double eLow = eLowInput;
		double eHigh = eHighInput;

		// now we can bisect the energy
		int imax = 100; // maximum of 100 bisection, reduces by 2**100
		for (int i = 1; i <= imax; i++) {
			e = (eHigh + eLow) * 0.5d;
			if ((eHigh - eLow) / Math.abs(e) < target_relative_error) {
				break;
			}
			if (countNodes(grid, e, v) > nOfDesiredNodes) {
				eHigh = e;
			} else {
				eLow = e;
			}
		}
		return e;
	}

	private static int countNodes(UniformGrid grid, double e, Potential v) {

		// first (leftmost) point
		double x = grid.xmin;
		double f0 = 0.d;

		// second point
		x += grid.h;
		double f1 = 0.0000001d; // arbitrary initial value
		int nOfNodes = 0;
		while (true) {
			x = x + grid.h;
			double f2 = 2.0d * f1 * (1.0d - grid.h * grid.h * mass * (e - v.value(x))) - f0;
			if (f1 * f2 <= 0.d) {
				nOfNodes++;
			}
			f0 = f1;
			f1 = f2;
			if (x >= grid.xmax) {
				break;
			}
		}

		return nOfNodes;
	}

}
