package net.tinvention.schrodinger.potential;

import net.tinvention.schrodinger.grid.Grid;

public class DressedPotential {
	private ClampedNucleiPotential potential;
	private double mass;
	private double energy;
	private Grid grid;

	public DressedPotential(ClampedNucleiPotential v, double mass, double energy, Grid grid) {
		this.potential = v;
		this.mass = mass;
		this.energy = energy;
		this.grid = grid;
	}

	public double value(double y) {
		double r = this.grid.mappingFunctionRofY(y);
		double potentialQ = -2.d * this.mass * (this.energy - this.potential.value(r));
		return potentialQ * Math.pow(this.grid.mappingFunctionGofY(y), 2) + this.grid.mappingFunctionFofY(y);
	}

}
