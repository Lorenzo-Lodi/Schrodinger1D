package net.tinvention.schrodinger.grid;

public class QuadraticGrid extends Grid {
	double rRef;

	public QuadraticGrid(double rMin, double rMax, int numberOfPoints, double rRef) {
		this.rRef = rRef;
		initializeClassVariables(rMin, rMax, numberOfPoints);
	}

	@Override
	public double mappingFunctionYofR(double r) {
		return r + Math.pow(r / rRef, 2);
	}

	@Override
	public double mappingFunctionRofY(double y) {
		return 0.5d * rRef * (-rRef + Math.sqrt(rRef * rRef + 4.0d * y));
	}

	@Override
	public double mappingFunctionGofY(double y) {
		return rRef * Math.pow(rRef * rRef + 4.0d * y, -0.5d);
	}

	@Override
	public double mappingFunctionFofY(double y) {
		return 3.0d*Math.pow(rRef * rRef + 4.0d * y, -2.0d);
	}

}
