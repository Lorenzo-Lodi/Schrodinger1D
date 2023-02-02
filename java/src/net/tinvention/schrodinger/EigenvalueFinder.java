package net.tinvention.schrodinger;

import net.tinvention.schrodinger.grid.Grid;
import net.tinvention.schrodinger.integrator.Integrator;
import net.tinvention.schrodinger.potential.DressedPotential;
import net.tinvention.schrodinger.potential.ClampedNucleiPotential;

public class EigenvalueFinder {
	private static final double TARGET_RELATIVE_ERROR = 10000.d * Math.ulp(1.d); // change for single-precision float
	private Grid grid;
	private ClampedNucleiPotential v;
	private double mass;
	private Integrator integrator;

	public EigenvalueFinder(Grid grid, ClampedNucleiPotential v, double mass, Integrator integrator) {
		this.grid = grid;
		this.v = v;
		this.mass = mass;
		this.integrator = integrator;
	}

	public EnergyLevel findEigenvalue(int nOfDesiredNodes) {

		EnergyLevel level = this.computeApproximateEnergyLevel(nOfDesiredNodes);

		// now we can bisect the energy
		int imax = 100; // maximum of 100 bisection, reduces by 2**100
		for (int i = 1; i <= imax; i++) {
			level.energy = (level.upperBound + level.lowerBound) * 0.5d;
			if ((level.upperBound - level.lowerBound) / Math.abs(level.energy) < TARGET_RELATIVE_ERROR) {
				break;
			}
			if (countNodes(level.energy) > nOfDesiredNodes) {
				level.upperBound = level.energy;
			} else {
				level.lowerBound = level.energy;
			}
		}
		return level;
	}

// y'' == -2m*(e-v) y
// (y_{x+h}+y_{x-h}-2y_x)/(h**2) ==  -2m*(e-v) y 	
// y_{x+h}+y_{x-h}-2y_x ==  -2*(h**2 * m*(e-v) y_x   	
// y_{x+h} == 2(1 -(h**2 * m*(e-v)) -y_{x-h}   	
	private int countNodes(double energy) {

		// first (leftmost) point
		double x = grid.getFirstGridPointValue();
		double y0 = 0.d;

		// second point
		x += grid.getStepSize();
		double y1 = 0.0000001d; // arbitrary initial value
		int nOfNodes = 0;

		DressedPotential dressedPotential = new DressedPotential(v, mass, energy, grid);
		while (true) {
			x = x + grid.getStepSize();
			double y2 = integrator.propagate(x, y0, y1, grid.getStepSize(), dressedPotential);
			if (y1 * y2 <= 0.d) {
				nOfNodes++;
			}
			y0 = y1;
			y1 = y2;
			if (x >= grid.getLastGridPointValue()) {
				break;
			}
		}

		return nOfNodes;
	}

	private EnergyLevel computeApproximateEnergyLevel(int nOfDesiredNodes) {
		EnergyLevel result = new EnergyLevel();

		result.upperBound = -Double.MAX_VALUE / 10.d;
		result.lowerBound = Double.MAX_VALUE / 10.d;
		double x = grid.getFirstGridPointValue();

		while (true) {
			if (x >= grid.getLastGridPointValue()) {
				break;
			}
			if (v.value(x) > result.upperBound) {
				result.upperBound = v.value(x);
			}
			if (v.value(x) < result.lowerBound) {
				result.lowerBound = v.value(x);
			}
			x += grid.getStepSize();
		}

		result.lowerBound = result.lowerBound - 0.05d * (result.upperBound - result.lowerBound);
		result.upperBound = result.upperBound + 0.05d * (result.upperBound - result.lowerBound);
		result.energy = (result.upperBound + result.lowerBound) / 2.d;

		return result;
	}

}
