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

	public EnergyLevel(Grid grid) {
		this.grid = grid;
		this.psi = new double[grid.getNumberOfPoints()];
	}

	public double normalizePsiTrapezoidalRule() {

		double sum = 0.5 * psiTimesG(0);
		for (int i = 1; i < grid.getNumberOfPoints() - 1; i++) {
			sum += psiTimesG(i);
		}
		sum += 0.5 * psiTimesG(grid.getNumberOfPoints() - 1);
		sum = sum * grid.getStepSizeYCoordinate();

		double normalizationFactor = 1. / Math.sqrt(sum);
		for (int i = 0; i < grid.getNumberOfPoints(); i++) {
			psi[i] = psi[i] * normalizationFactor;
		}
		return normalizationFactor;
	}

	public double normalizePsiSimpsonsOneThirdRule() {
		if (grid.getNumberOfPoints() % 2 == 0) {
			System.out.println("WARNING! At the moment simpson's rule works only for odd number of points");
		}

		double sum = psiTimesG(0);
		for (int i = 1; i < grid.getNumberOfPoints() - 1; i += 2) {
			sum += 4.0d * psiTimesG(i);
		}
		for (int i = 2; i < grid.getNumberOfPoints() - 1; i += 2) {
			sum += 2.0d * psiTimesG(i);
		}
		sum += psiTimesG(grid.getNumberOfPoints() - 1);
		sum = sum * grid.getStepSizeYCoordinate() / 3.0d;

		double normalizationFactor = 1. / Math.sqrt(sum);
		for (int i = 0; i < grid.getNumberOfPoints(); i++) {
			psi[i] = psi[i] * normalizationFactor;
		}
		return normalizationFactor;
	}

	private double psiTimesG(int i) {
		return (i < 0 || i > grid.getNumberOfPoints() - 1) ? 0.0d : psi[i] * grid.mappingFunctionGofY(i);
	}
}
