package net.tinvention.schrodinger.grid;

public class SurkusGrid extends Grid {
	private double rRef;
	private double alpha;

	public SurkusGrid(double rMin, double rMax, int numberOfPoints, double rRef, double alpha) {
		this.rRef = rRef;
		this.alpha = alpha;
		initializeClassVariables(rMin, rMax, numberOfPoints);
	}

	@Override
	public double mappingFunctionYofR(double r) {
		return (Math.pow(r / rRef, alpha) - 1.0d) / (Math.pow(r / rRef, alpha) + 1.0d);
	}

	@Override
	public double mappingFunctionRofY(double y) {
		return rRef * Math.pow((1 + y) / (1 - y), 1.0d / alpha);
	}

	@Override
	public double mappingFunctionGofY(double y) {
		return (2.0d * rRef * Math.pow(1 + y, -1.0d + 1.0d / alpha))
				/ (alpha * Math.pow(1.0d - y, 1.0d + 1.0d / alpha));
	}

	@Override
	public double mappingFunctionFofY(double y) {
		return (1.0d - 1.0d / Math.pow(alpha, 2)) / Math.pow(1 - y * y, 2);
	}

}
