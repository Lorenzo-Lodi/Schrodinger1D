package net.tinvention.schrodinger.grid;

public class LogarithmicGrid extends Grid {
	double rRef;

	public LogarithmicGrid(double rMin, double rMax, int numberOfPoints, double rRef) {
		this.rRef = rRef;
		initializeClassVariables(rMin, rMax, numberOfPoints);
	}

	@Override
	public double mappingFunctionYofR(double r) {
		return rRef * Math.log(1 + r / rRef);
	}

	@Override
	public double mappingFunctionRofY(double y) {
		return rRef * (-1. + Math.exp(y / rRef));
	}

	@Override
	public double mappingFunctionGofY(double y) {
		return Math.exp(y / rRef);
	}

	@Override
	public double mappingFunctionFofY(double y) {
		return -1. / (4. * rRef * rRef);
	}

}
