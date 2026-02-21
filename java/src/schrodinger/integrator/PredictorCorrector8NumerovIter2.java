package schrodinger.integrator;

public class PredictorCorrector8NumerovIter2 extends PredictorCorrector8Abstract {

    public PredictorCorrector8NumerovIter2() {
        super(new Numerov(), 2);
    }

}
