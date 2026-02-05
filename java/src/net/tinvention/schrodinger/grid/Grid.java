package net.tinvention.schrodinger.grid;

public abstract class Grid {
	double rMin;
	double rMax;
	double yMin;
	double yMax;
	double stepSizeYCoordinate;
	int numberOfPoints;

	// Don't make this a constructor because we have to run it AFTER having done
	// initializations necessary
	// for functions mappingFunctionYofR and mappingFunctionYofR
	void initializeClassVariables(double rMin, double rMax, int numberOfPoints) {
		this.rMin = Math.min(rMin, rMax);
		this.rMax = Math.max(rMin, rMax);
		this.yMax = this.mappingFunctionYofR(this.rMax);
		this.yMin = this.mappingFunctionYofR(this.rMin);
		this.numberOfPoints = Math.max(numberOfPoints, 2);
		this.stepSizeYCoordinate = (yMax - yMin) / (numberOfPoints - 1);
	}

	public double getFirstYValue() {
		return getYValue(0);
	}

	public double getLastYValue() {
		return yMax;
	}

	public double getStepSizeYCoordinate() {
		return stepSizeYCoordinate;
	}

	public int getNumberOfPoints() {
		return numberOfPoints;
	}

	/**
	 * 
	 * @param i Index going from 0 to numberOfPoints-1
	 * @return Value of the i-th grid point
	 */
	public double getYValue(int i) {
		return yMin + getStepSizeYCoordinate() * i;
	}

	public double getRValue(int i) {
		return mappingFunctionRofY(getYValue(i));
	}

	/**
	 * The mapping function y(r), mapping the r coordinate to the transformed y
	 * coordinate
	 * 
	 * @param r Value of the r coordinate
	 * @return Corresponding value y(r)
	 */
	public abstract double mappingFunctionYofR(double r);

	/**
	 * The inverse mapping function r(y), mapping the transformed y coordinate to r
	 * 
	 * @param y Value of the mapped coordinate
	 * @return Corresponding value r(y)
	 */
	public abstract double mappingFunctionRofY(double y);

	/**
	 * It's the derivative of the r(y) function, referred to as g(y) in Meshkov2008
	 * 
	 * @param y Value of the mapped coordinate
	 * @return Value of g(y)=r'(y)
	 */
	public abstract double mappingFunctionGofY(double y);

	public double mappingFunctionGofY(int i) {
		return mappingFunctionGofY(getYValue(i));
	}

	/**
	 * It's the function F(y) of eq. (9) in Meshkov2008, the additional term to the
	 * transformed function \tilde{Q}
	 * 
	 * @param y
	 * @return
	 */
	public abstract double mappingFunctionFofY(double y);

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("i           R                Y                g(y)             F(y)\n");
		for (int i = 0; i < getNumberOfPoints(); i++) {
			sb.append(String.format("%8d", i));
			sb.append("   ");
			double y = getYValue(i);
			double r = mappingFunctionRofY(y);
			if (r >= 0.d) {
				sb.append(" ");
			}
			sb.append(String.format("%.012f", r) + "  ");
			if (y >= 0.d) {
				sb.append(" ");
			}
			sb.append(String.format("%.012f", y));
			sb.append("  ");
			sb.append(String.format("%.012f", mappingFunctionGofY(y)));
			sb.append("  ");
			sb.append(String.format("%.012f", mappingFunctionFofY(y)) + "\n");
		}
		return sb.toString();
	}

	public static Grid generateUniformGrid(double rMin, double rMax, int numberOfPoints) {
		return new UniformGrid(rMin, rMax, numberOfPoints);
	}

	public static Grid generateLogarithmicGrid(double rMin, double rMax, int numberOfPoints, double rRef) {
		return new LogarithmicGrid(rMin, rMax, numberOfPoints, rRef);
	}

	public static Grid generateSurkusGrid(double rMin, double rMax, int numberOfPoints, double rRef, double alpha) {
		return new SurkusGrid(rMin, rMax, numberOfPoints, rRef, alpha);
	}

	public static Grid generateSqrtGrid(double rMin, double rMax, int numberOfPoints, double rRef) {
		return new SqrtGrid(rMin, rMax, numberOfPoints, rRef);
	}

}
