package net.tinvention.schrodinger;

import net.tinvention.schrodinger.grid.UniformGrid;
import net.tinvention.schrodinger.potential.Potential;

public class EigenvalueFinder {
	private static final double mass = 2.0d;
	private static final double targetRelativeError = Math.ulp(1.d); // change for single-precision float
	UniformGrid grid;
	Potential v;

	public EigenvalueFinder(UniformGrid grid, Potential v) {
		this.grid = grid;
		this.v = v;
	}

	public EnergyLevel computeApproximateEnergyLevel(int nOfDesiredNodes) {
		EnergyLevel result = new EnergyLevel();

		result.upperBound = -Double.MAX_VALUE / 1000.d;
		result.lowerBound = Double.MAX_VALUE / 1000.d;
		double x = grid.xmin;

		while (true) {
			if (x >= grid.xmax) {
				break;
			}
			if (v.value(x) > result.upperBound) {
				result.upperBound = v.value(x);
			}
			if (v.value(x) < result.lowerBound) {
				result.lowerBound = v.value(x);
			}
			x += grid.h;
		}

		result.lowerBound = result.lowerBound - 0.05d * (result.upperBound - result.lowerBound);
		result.upperBound = result.upperBound + 0.05d * (result.upperBound - result.lowerBound);
		result.energy = (result.upperBound + result.lowerBound) / 2.d;

		return result;
	}

	public double findEigenvalue(int nOfDesiredNodes, UniformGrid grid, Potential v) {

		EnergyLevel trial = this.computeApproximateEnergyLevel(nOfDesiredNodes);

		double e = trial.energy;
		double eLow = trial.lowerBound;
		double eHigh = trial.upperBound;

		// now we can bisect the energy
		int imax = 100; // maximum of 100 bisection, reduces by 2**100
		for (int i = 1; i <= imax; i++) {
			e = (eHigh + eLow) * 0.5d;
			if ((eHigh - eLow) / Math.abs(e) < targetRelativeError) {
				break;
			}
			if (countNodes(e) > nOfDesiredNodes) {
				eHigh = e;
			} else {
				eLow = e;
			}
		}
		return e;
	}

	private int countNodes(double e) {

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
