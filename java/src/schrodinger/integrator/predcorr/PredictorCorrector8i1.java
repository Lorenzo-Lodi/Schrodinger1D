package schrodinger.integrator.predcorr;

import schrodinger.integrator.stormer.Stormer6;

public class PredictorCorrector8i1 extends PredictorCorrector8Abstract {

    public PredictorCorrector8i1() {
        super(new Stormer6(), 1);
    }

}
