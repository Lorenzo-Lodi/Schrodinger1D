package schrodinger.integrator.predcorr;

import schrodinger.integrator.stormer.Stormer6;

public class PredictorCorrector8NumerovIter1 extends PredictorCorrector8Abstract {

    public PredictorCorrector8NumerovIter1() {
        super(new Stormer6(), 1);
    }

}
