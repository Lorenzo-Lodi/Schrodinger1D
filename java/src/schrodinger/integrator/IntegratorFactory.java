package schrodinger.integrator;

import schrodinger.pt_correction.NumerovPTCorrection;
import schrodinger.pt_correction.PTCorrection;
import schrodinger.pt_correction.TaylorThreePointsPTCorrection;

/**
 * Factory for creating integrators paired with their corresponding perturbative corrections.
 * This ensures that the correct correction method is always used with each integrator.
 */
public class IntegratorFactory {

    /**
     * Returns an integrator paired with its perturbative correction.
     * 
     * @param integratorType The type of integrator to create
     * @return A pair containing the integrator and its corresponding correction
     */
    public static IntegratorPair createIntegrator(String integratorType) {
        switch (integratorType) {
            case "Numerov":
                return new IntegratorPair(new Numerov(), new NumerovPTCorrection());
            case "TaylorThreePoints":
                return new IntegratorPair(new TaylorThreePoints(), new TaylorThreePointsPTCorrection());
            case "PredictorCorrector6":
                return new IntegratorPair(new PredictorCorrector6(), null);
            case "PredictorCorrector8":
                return new IntegratorPair(new PredictorCorrector8(), null);
            case "PredictorCorrector8Alt":
                return new IntegratorPair(new PredictorCorrector8Alt(), null);
            case "RaptisAllison":
                return new IntegratorPair(new RaptisAllison(), null);
            case "Stormer5":
                return new IntegratorPair(new Stormer5(), null);
            case "Stormer6":
                return new IntegratorPair(new Stormer6(), null);
            case "ExponentiallyFitted":
                return new IntegratorPair(new ExponentiallyFitted(), null);
            case "Vignoli":
                return new IntegratorPair(new Vignoli(), null);
            default:
                throw new IllegalArgumentException("Unknown integrator type: " + integratorType);
        }
    }

    // Backward-compatible factory methods (return just the integrator)
    public static Integrator getNumerov() {
        return new Numerov();
    }

    public static Integrator getTaylorThreePoints() {
        return new TaylorThreePoints();
    }

    public static Integrator getPredictorCorrector6() {
        return new PredictorCorrector6();
    }

    public static Integrator getPredictorCorrector8() {
        return new PredictorCorrector8();
    }

    public static Integrator getPredictorCorrector8Alt() {
        return new PredictorCorrector8Alt();
    }

    public static Integrator getRaptisAllison() {
        return new RaptisAllison();
    }

    public static Integrator getStormer7() {
        return new Stormer5();
    }

    public static Integrator getStormer8() {
        return new Stormer6();
    }

    public static Integrator getExponentiallyFitted() {
        return new ExponentiallyFitted();
    }

    public static Integrator getVignoli() {
        return new Vignoli();
    }

    /**
     * Pair containing an integrator and its perturbative correction.
     */
    public static class IntegratorPair {
        private final Integrator integrator;
        private final PTCorrection correction;

        public IntegratorPair(Integrator integrator, PTCorrection correction) {
            this.integrator = integrator;
            this.correction = correction;
        }

        public Integrator getIntegrator() {
            return integrator;
        }

        public PTCorrection getCorrection() {
            return correction;
        }
    }
}