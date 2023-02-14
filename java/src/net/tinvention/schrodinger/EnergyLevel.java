package net.tinvention.schrodinger;

import net.tinvention.schrodinger.grid.Grid;

public class EnergyLevel {
	public int numberOfNodes;
	public double energy;
	public double upperBound;
	public double lowerBound;
	public int numberOfBisections;
	public double[] psi;
	public double perturbativeCorrectionToEnergy;
	private Grid grid;
	private boolean isPsiNormalized = false;

	public EnergyLevel(Grid grid) {
		this.grid = grid;
		this.psi = new double[grid.getNumberOfPoints()];
	}

	// Note: because the wavefunctions are exponentially decreasing (or faster), the
	// trapezoidal rule is very quickly convergent (exponentially or so) until the
	// truncation error
	// at the borders takes over the global error.
	// The rectangle rule is practically the same.
	public double normalizePsi() {
		double sum = 0.0d;
		for (int i = 1; i < grid.getNumberOfPoints() - 1; i++) {
			sum += psiTimesGSquared(i);
		}
		sum += 0.5 * (psiTimesGSquared(0) + psiTimesGSquared(grid.getNumberOfPoints() - 1));

		sum = sum * grid.getStepSizeYCoordinate();

		double normalizationFactor = 1. / Math.sqrt(sum);
		for (int i = 0; i < grid.getNumberOfPoints(); i++) {
			psi[i] = psi[i] * normalizationFactor;
		}
		isPsiNormalized = true;
		return normalizationFactor;
	}

	private double psiTimesGSquared(int i) {
		return Math.pow(psi[i] * grid.mappingFunctionGofY(i), 2);
	}

	public boolean isPsiNormalized() {
		return isPsiNormalized;
	}

}
