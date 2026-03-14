package schrodinger.integrator.predcorr;

import schrodinger.integrator.stormer.Stormer6;

public class PredictorCorrector8i2 extends PredictorCorrector8Abstract {

    public PredictorCorrector8i2() {
        super(new Stormer6(), 2);
    }

}
