package net.tinvention.schrodinger.grid;

public class SqrtGrid extends Grid {
	double rRef;

	public SqrtGrid(double rMin, double rMax, int numberOfPoints, double rRef) {
		this.rRef = rRef;
		initializeClassVariables(rMin, rMax, numberOfPoints);
	}

	@Override
	public double mappingFunctionYofR(double r) {
		return 2.d * rRef * (Math.sqrt(1 + r / rRef) - 1.d);
	}

	@Override
	public double mappingFunctionRofY(double y) {
		return y * (1. + y / (4.0d * rRef));
	}

	@Override
	public double mappingFunctionGofY(double y) {
		return 1.0d + y / (2.0d * rRef);
	}

	@Override
	public double mappingFunctionFofY(double y) {
		return -3.0 / Math.pow(8.0d * rRef + 4.0d * y, 2);
	}

}
