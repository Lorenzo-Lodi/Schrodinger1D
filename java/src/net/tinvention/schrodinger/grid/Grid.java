package net.tinvention.schrodinger.grid;

public interface Grid {

	public double getFirstGridPointValue();

	public double getLastGridPointValue();

	public double getStepSize();

	public int getNumberOfPoints();

	/**
	 * 
	 * @param i Index going from 0 to numberOfPoints-1
	 * @return Value of the i-th grid point
	 */
	public default double getGridValue(int i) {
		return getFirstGridPointValue() + getStepSize() * i;
	}

	/**
	 * The mapping function y(r), mapping the r coordinate to the transformed y
	 * coordinate
	 * 
	 * @param r Value of the r coordinate
	 * @return Corresponding value y(r)
	 */
	public double mappingFunctionYofR(double r);

	/**
	 * The inverse mapping function r(y), mapping the transformed y coordinate to r
	 * 
	 * @param y Value of the mapped coordinate
	 * @return Corresponding value r(y)
	 */
	public double mappingFunctionRofY(double y);

	/**
	 * It's the derivative of the r(y) function, referred to as g(y) in Meshkov2008
	 * 
	 * @param y Value of the mapped coordinate
	 * @return Value of g(y)=r'(y)
	 */
	public double mappingFunctionGofY(double y);

	/**
	 * It's the function F(y) of eq. (9) in Meshkov2008, the additional term to the
	 * dressed potential \tilde{Q}
	 * 
	 * @param y
	 * @return
	 */
	public double mappingFunctionFofY(double y);

	public default String printGrid() {
		StringBuilder sb = new StringBuilder();
		sb.append("i           Y                R\n");
		for (int i = 0; i < getNumberOfPoints(); i++) {
			sb.append(String.format("%8d", i));
			sb.append("   ");
			double r = getGridValue(i);
			double y = mappingFunctionYofR(r);
			if (y >= 0.d) {
				sb.append(" ");
			}
			sb.append(String.format("%,.012f", y) + "  ");
			if (r >= 0.d) {
				sb.append(" ");
			}
			sb.append(String.format("%,.012f", r) + "\n");
		}
		return sb.toString();
	}

}
