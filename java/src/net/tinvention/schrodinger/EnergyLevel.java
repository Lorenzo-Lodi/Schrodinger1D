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

	public void normalizePsi() {

		double sum = 0.5 * psi[0] * grid.mappingFunctionGofY(grid.getFirstYValue());
		for (int i = 1; i < grid.getNumberOfPoints() - 1; i++) {
			sum += psi[i] * grid.mappingFunctionGofY(grid.getYValue(i));
		}
		sum += 0.5 * psi[grid.getNumberOfPoints() - 1] * grid.mappingFunctionGofY(grid.getLastYValue());
		sum = sum * grid.getStepSizeYCoordinate();

		for (int i = 0; i < grid.getNumberOfPoints(); i++) {
			psi[i] = psi[i] / Math.sqrt(sum);
		}
	}
}
