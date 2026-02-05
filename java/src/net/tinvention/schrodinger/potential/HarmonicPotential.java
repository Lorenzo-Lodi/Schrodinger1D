package net.tinvention.schrodinger.potential;

public class HarmonicPotential implements PhysicalPotential {
	private final double r0;
	private final double alpha;

	public HarmonicPotential(double r0, double alpha) {
		this.r0 = r0;
		this.alpha = alpha;

	}

	public double value(double r) {
		return alpha * (r - r0) * (r - r0);
	}

}
