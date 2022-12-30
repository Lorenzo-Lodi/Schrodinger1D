package net.tinvention;

public class ThreePointProcedural {
	private static double xmin = -5.5d;
	private static double xmax = -xmin;
	private static double h = 0.1d;
	private static double mass = 2.0d;

	public static void main(String[] args) {

		double eHigh = -Double.MAX_VALUE / 1000.d;
		double eLow = Double.MAX_VALUE / 1000.d;
		double x = xmin;

		while (true) {
			if (x >= xmax) {
				break;
			}
			if (V(x) > eHigh) {
				eHigh = V(x);
			}
			if (V(x) < eLow) {
				eLow = V(x);
			}
			x += h;
		}

		eLow = eLow - 0.05d * (eHigh - eLow);
		eHigh = eHigh + 0.05d * (eHigh - eLow);
		double eTrial = (eHigh + eLow) / 2.d;

		for (int i = 0; i < 20; i++) {
			h = Math.pow(10.d, -1.d - ((double) i) * 0.25d);
			for (int nOfDesiredNodes = 0; nOfDesiredNodes < 1; nOfDesiredNodes++) {
				double eExact = 0.5d + nOfDesiredNodes;
				double ek = findEigenvalue(nOfDesiredNodes, eLow, eHigh, eTrial);
				System.out.println(i + " " + h + " " + ek);
			}
		}
	}

	private static double findEigenvalue(int nOfDesiredNodes, double eLowInput, double eHighInput, double eTrial) {
		double target_relative_error = Math.ulp(1.d); // change for single-precision float
		double e = eTrial;
		double eLow = eLowInput;
		double eHigh = eHighInput;

		// now we can bisect the energy
		int imax = 100; // maximum of 100 bisection, reduces by 2**100
		for (int i = 1; i <= imax; i++) {
			e = (eHigh + eLow) * 0.5d;
			if ((eHigh - eLow) / Math.abs(e) < target_relative_error) {
				break;
			}
			if (countNodes(xmin, xmax, h, e) > nOfDesiredNodes) {
				eHigh = e;
			} else {
				eLow = e;
			}
		}
		return e;
	}

	private static int countNodes(double xmin2, double xmax2, double h2, double e) {

		// first (leftmost) point
		double x = xmin;
		double f0 = 0.d;

		// second point
		x += h;
		double f1 = 0.0000001d; // arbitrary initial value
		int nOfNodes = 0;
		while (true) {
			x = x + h;
			double f2 = 2.0d * f1 * (1.0d - h * h * mass * (e - V(x))) - f0;
			if (f1 * f2 <= 0.d) {
				nOfNodes++;
			}
			f0 = f1;
			f1 = f2;
			if (x >= xmax) {
				break;
			}
		}

		return nOfNodes;
	}

	private static double V(double x) {
		return x * x;
	}

}
