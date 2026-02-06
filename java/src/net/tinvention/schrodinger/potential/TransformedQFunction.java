package net.tinvention.schrodinger.potential;

import net.tinvention.schrodinger.grid.Grid;

/**
 * Transformed Q-function Q̃(y) for the mapped equation: ϕ''(y) = -Q̃(y)ϕ(y)
 * MSL Eq. (8): Q̃(y) = g²(y)·Q(r(y)) + F(y)
 */
public class TransformedQFunction {
	private final PhysicalPotential physicalPotential;
	private final double mass;
	private final double energy;
	private final Grid grid;

	public TransformedQFunction(PhysicalPotential physicalPotential, double mass, double energy, Grid grid) {
		this.physicalPotential = physicalPotential;
		this.mass = mass;
		this.energy = energy;
		this.grid = grid;
	}

	public double value(double y) {
		double r = this.grid.r(y);
		double potentialQ = 2.d * this.mass * (this.energy - this.physicalPotential.value(r));
		return potentialQ * Math.pow(this.grid.g(y), 2) + this.grid.F(y);
	}

	public PhysicalPotential getPhysicalPotential() {
		return physicalPotential;
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
