package net.tinvention.schrodinger.potential;

public class HarmonicPotential implements ClampedNucleiPotential {

	public double value(double x) {
		return x * x;
	}

}
