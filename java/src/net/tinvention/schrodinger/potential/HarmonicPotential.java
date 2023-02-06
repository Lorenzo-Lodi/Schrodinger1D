package net.tinvention.schrodinger.potential;

public class HarmonicPotential implements BarePotential {
	private double r0;
	private double alpha;

	public HarmonicPotential(double r0, double alpha) {
		this.r0 = r0;
		this.alpha = alpha;

	}

	public double value(double r) {
		return alpha * (r - r0) * (r - r0);
	}

}
