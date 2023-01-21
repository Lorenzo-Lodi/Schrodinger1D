package net.tinvention.schrodinger.grid;

public class UniformGrid {
	public double xmin = -5.5d;
	public double xmax = -xmin;
	public double step = 0.1d;
	public int numberOfPoints;

	public UniformGrid(double xmin, double xmax, double h) {
		this.xmin = xmin;
		this.xmax = xmax;
		this.numberOfPoints = (int) ((xmax - xmin) / h);
		this.step = (xmax - xmin) / ((double) (numberOfPoints + 1));
	}

	/**
	 * 
	 * @param i Index going from 0 to numberOfPoints-1
	 * @return Value of the i-th grid point
	 */
	public double getGridValue(int i) {
		return xmin + ((double) i) * step;
	}
	
}
