package net.tinvention.schrodinger.potential;

public class DressedPotential {
	Potential v;
	double mass;
	double energy;

	public DressedPotential(Potential v, double mass, double energy) {
		this.v = v;
		this.mass = mass;
		this.energy = energy;
	}

	public double value(double x) {
		return -2.d * this.mass * (this.energy - this.v.value(x));
	}

}
