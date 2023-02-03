package net.tinvention.schrodinger.grid;

public class UniformGrid extends Grid {
	
	public UniformGrid(double rMin, double rMax, int numberOfPoints) {
		initializeClassVariables(rMin, rMax, numberOfPoints);
	}

	public UniformGrid(double rMin, double rMax, double desiredStepSize) {
		initializeClassVariables(rMin, rMax, (int) (Math.abs(rMax - rMin) / desiredStepSize));
	}

	@Override
	public double mappingFunctionYofR(double r) {
		return r;
	}

	@Override
	public double mappingFunctionRofY(double y) {
		return y;
	}

	@Override
	public double mappingFunctionGofY(double y) {
		return 1.d;
	}

	@Override
	public double mappingFunctionFofY(double y) {
		return 0;
	}

}
