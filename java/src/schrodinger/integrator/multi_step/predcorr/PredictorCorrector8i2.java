package schrodinger.integrator.multi_step.predcorr;

import schrodinger.integrator.multi_step.Stormer6;

public class PredictorCorrector8i2 extends PredictorCorrector8Abstract {

    public PredictorCorrector8i2() {
        super(new Stormer6(), 2);
    }

}
