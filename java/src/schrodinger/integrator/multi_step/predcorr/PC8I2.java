package schrodinger.integrator.multi_step.predcorr;

import schrodinger.integrator.multi_step.Stormer6;

public class PC8I2 extends PC8Abstract {

    public PC8I2() {
        super(new Stormer6(), 2);
    }

}
