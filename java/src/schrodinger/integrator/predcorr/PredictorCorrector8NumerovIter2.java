package schrodinger.integrator.predcorr;

import schrodinger.integrator.stormer.Stormer6;

public class PredictorCorrector8NumerovIter2 extends PredictorCorrector8Abstract {

    public PredictorCorrector8NumerovIter2() {
        super(new Stormer6(), 2);
    }

}
