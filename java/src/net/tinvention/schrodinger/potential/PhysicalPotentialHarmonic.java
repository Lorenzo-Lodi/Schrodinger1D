package net.tinvention.schrodinger.potential;

public class PhysicalPotentialHarmonic implements PhysicalPotential {
	private final double r0;
	private final double alpha;

	public PhysicalPotentialHarmonic(double r0, double alpha) {
		this.r0 = r0;
		this.alpha = alpha;

	}

	public double value(double r) {
		return alpha * (r - r0) * (r - r0);
	}

}
