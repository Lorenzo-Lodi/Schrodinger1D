package net.tinvention.schrodinger.potential;

public interface Potential {
/**
 * 
 * @param x position in bohrs 
 * @return The energy at point x, in hartrees 
 */
	public double value(double x);

}
