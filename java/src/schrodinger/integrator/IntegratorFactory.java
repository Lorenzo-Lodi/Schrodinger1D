package schrodinger.integrator;

import schrodinger.integrator.expfitted.ExponentiallyFitted;
import schrodinger.integrator.expfitted.Numerov;
import schrodinger.integrator.expfitted.EfnFixedBeta;
import schrodinger.integrator.predcorr.PredictorCorrector6;
import schrodinger.integrator.predcorr.PredictorCorrector8NumerovIter1;
import schrodinger.integrator.predcorr.PredictorCorrector8NumerovIter2;
import schrodinger.integrator.stormer.Stormer5;
import schrodinger.integrator.stormer.Stormer6;
import schrodinger.integrator.stormer.Stormer8;
import schrodinger.pt_correction.PTCorrection;

import java.util.ArrayList;
import java.util.List;

public class IntegratorFactory {

    public static Integrator getTaylorThreePoints() {
        return new TaylorThreePoints();
    }

    public static Integrator getNumerov() {
        return new Numerov();
    }

    public static Integrator getEfnFixedBeta() {
        return new EfnFixedBeta();
    }

    public static Integrator getExponentiallyFitted() {
        return new ExponentiallyFitted();
    }

    public static Integrator getStormer5() {
        return new Stormer5();
    }

    public static Integrator getStormer6() {
        return new Stormer6();
    }

    public static Integrator getStormer8() {
        return new Stormer8();
    }

    public static Integrator getPredictorCorrector6() {
        return new PredictorCorrector6();
    }

    public static Integrator getPredictorCorrector8NumerovIter1() {
        return new PredictorCorrector8NumerovIter1();
    }

    public static Integrator getPredictorCorrector8NumerovIter2() {
        return new PredictorCorrector8NumerovIter2();
    }

    public static List<Integrator> getAll() {
        List<Integrator> list = new ArrayList<>();
        list.add(getTaylorThreePoints());
        list.add(getNumerov());
        list.add(getEfnFixedBeta());
        list.add(getExponentiallyFitted());
        list.add(getStormer5());
        list.add(getStormer6());
        list.add(getStormer8());
        list.add(getPredictorCorrector6());
        list.add(getPredictorCorrector8NumerovIter1());
        list.add(getPredictorCorrector8NumerovIter2());
        return list;
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