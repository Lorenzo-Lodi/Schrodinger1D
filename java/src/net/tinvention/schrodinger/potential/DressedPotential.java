package net.tinvention.schrodinger.potential;

import net.tinvention.schrodinger.grid.Grid;

public class DressedPotential {
	private BarePotential barePotential;
	private double mass;
	private double energy;
	private Grid grid;

	public DressedPotential(BarePotential barePotential, double mass, double energy, Grid grid) {
		this.barePotential = barePotential;
		this.mass = mass;
		this.energy = energy;
		this.grid = grid;
	}

	public double value(double y) {
		double r = this.grid.mappingFunctionRofY(y);
		double potentialQ = 2.d * this.mass * (this.energy - this.barePotential.value(r));
		return potentialQ * Math.pow(this.grid.mappingFunctionGofY(y), 2) + this.grid.mappingFunctionFofY(y);
	}

	public BarePotential getPotential() {
		return barePotential;
	}

	public double getMass() {
		return mass;
	}

	public double getEnergy() {
		return energy;
	}

	public Grid getGrid() {
		return grid;
	}

}
