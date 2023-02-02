package net.tinvention.schrodinger.potential;

import net.tinvention.schrodinger.grid.Grid;

public class DressedPotential {
	private ClampedNucleiPotential v;
	private double mass;
	private double energy;
	private Grid grid;

	public DressedPotential(ClampedNucleiPotential v, double mass, double energy, Grid grid) {
		this.v = v;
		this.mass = mass;
		this.energy = energy;
		this.grid = grid;
	}

	public double value(double y) {
		double r = this.grid.mappingFunctionRofY(y);
		double potentialQ = -2.d * this.mass * (this.energy - this.v.value(r));
		return potentialQ * Math.pow(this.grid.mappingFunctionGofY(y), 2) + this.grid.mappingFunctionFofY(y);
	}

}
