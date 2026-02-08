package schrodinger.potential;

public interface PhysicalPotential {
	/**
	 * 
	 * @param x position in bohrs
	 * @return The energy at point x, in hartrees
	 */
	double value(double x);

}
