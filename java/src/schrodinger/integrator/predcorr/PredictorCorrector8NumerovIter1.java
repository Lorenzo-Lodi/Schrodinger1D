package schrodinger.integrator.predcorr;

import schrodinger.integrator.expfitted.Numerov;

public class PredictorCorrector8NumerovIter1 extends PredictorCorrector8Abstract {

    public PredictorCorrector8NumerovIter1() {
        super(new Numerov(), 1);
    }

}
