package net.tinvention.schrodinger.integrator;

import java.math.BigDecimal;

import net.tinvention.schrodinger.potential.DressedPotential;

public class TaylorThreePoints implements Integrator {
	public double propagate(double yOfN, double psiOfNMinusOne, double psiOfN, double stepSize,
			DressedPotential qTilde) {

//		BigDecimal qTildeBigDecimal = BigDecimal.valueOf(qTilde.value(yOfN));
//		BigDecimal stepSizeSquaredBigDecimal = BigDecimal.valueOf(stepSize * stepSize);
//		BigDecimal result = stepSizeSquaredBigDecimal.multiply(qTildeBigDecimal);
//		result = result.multiply(BigDecimal.valueOf(-1));
//		result = result.add(BigDecimal.valueOf(2));
//		result = result.multiply(BigDecimal.valueOf(psiOfN));
//		result = result.subtract(BigDecimal.valueOf(psiOfNMinusOne));
//		return result.doubleValue();
		return psiOfN * (2.0d - stepSize * stepSize * qTilde.value(yOfN)) - psiOfNMinusOne;
	}

}
