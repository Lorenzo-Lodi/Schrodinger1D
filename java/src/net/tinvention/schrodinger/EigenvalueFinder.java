package net.tinvention.schrodinger;

import net.tinvention.schrodinger.grid.Grid;
import net.tinvention.schrodinger.integrator.Integrator;
import net.tinvention.schrodinger.potential.DressedPotential;
import net.tinvention.schrodinger.potential.BarePotential;

public class EigenvalueFinder {
	private static final double TARGET_RELATIVE_ERROR = 2.d * Math.ulp(1.d); // change for single-precision float
	private static final int MAXIMUM_NUMBER_OF_BISECTIONS = 60; // reduces error by 2**n
	private Grid grid;
	private BarePotential barePotential;
	private double mass;
	private Integrator integrator;

	public EigenvalueFinder(Grid grid, BarePotential barePotential, double mass, Integrator integrator) {
		this.grid = grid;
		this.barePotential = barePotential;
		this.mass = mass;
		this.integrator = integrator;
	}

	public EnergyLevel findEigenvalueByBisection(int nOfDesiredNodes) {

		EnergyLevel level = this.computeApproximateEnergyLevel(nOfDesiredNodes);

		// now we can bisect the energy
		for (level.numberOfBisections = 1; level.numberOfBisections <= MAXIMUM_NUMBER_OF_BISECTIONS; level.numberOfBisections++) {
			level.energy = (level.upperBound + level.lowerBound) * 0.5d;
			propagatePsiLeftToRight(level);
			if ((level.upperBound - level.lowerBound) / Math.abs(level.energy) < TARGET_RELATIVE_ERROR) {
				break;
			}
			if (level.countNumberOfNodes() > nOfDesiredNodes) {
				level.upperBound = level.energy;
			} else {
				level.lowerBound = level.energy;
			}
		}

		DressedPotential dressedPotential = new DressedPotential(barePotential, mass, level.energy, grid);
		level.normalizePsi();
		level.perturbativeCorrectionToEnergy = integrator.computePerturbativeCorrection(level, dressedPotential);

		return level;
	}

	private void propagatePsiLeftToRight(EnergyLevel level) {
		// first (leftmost) point
		level.psi[0] = 0;

		// second point
		level.psi[1] = 0.0000001d; // arbitrary initial value

		DressedPotential dressedPotential = new DressedPotential(barePotential, mass, level.energy, grid);
		for (int n = 1; n < grid.getNumberOfPoints() - 1; n++) {
			level.psi[n + 1] = integrator.propagateForward(level.psi, n, dressedPotential);
		}
	}

	private EnergyLevel computeApproximateEnergyLevel(int nOfDesiredNodes) {
		EnergyLevel result = new EnergyLevel(grid);
		result.numberOfNodes = nOfDesiredNodes;
		result.upperBound = -Double.MAX_VALUE / 10.d;
		result.lowerBound = Double.MAX_VALUE / 10.d;
		double y = grid.getFirstYValue();

		while (true) {
			if (y >= grid.getLastYValue()) {
				break;
			}
			double r = grid.mappingFunctionRofY(y);
			double potentialValues = barePotential.value(r);
			if (potentialValues > result.upperBound) {
				result.upperBound = potentialValues;
			}
			if (potentialValues < result.lowerBound) {
				result.lowerBound = potentialValues;
			}
			y += grid.getStepSizeYCoordinate();
		}

		result.lowerBound = result.lowerBound - 0.05d * (result.upperBound - result.lowerBound);
		result.upperBound = result.upperBound + 0.05d * (result.upperBound - result.lowerBound);
		result.energy = (result.upperBound + result.lowerBound) / 2.d;

		return result;
	}

}
