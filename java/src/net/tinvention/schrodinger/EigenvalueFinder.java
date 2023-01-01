package net.tinvention.schrodinger;

import net.tinvention.schrodinger.grid.UniformGrid;
import net.tinvention.schrodinger.potential.Potential;

public class EigenvalueFinder {
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
	
}
