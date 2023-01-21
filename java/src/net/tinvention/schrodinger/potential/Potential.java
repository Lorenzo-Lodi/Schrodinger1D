package net.tinvention.schrodinger.potential;

public interface Potential {
	/**
	 * 
	 * @param x position in bohrs
	 * @return The energy at point x, in hartrees
	 */
	public double value(double x);

	public default double dressedValue(double x, double mass, double energy) {
		return 2.0d * mass * (value(x) - energy);
	}

}
