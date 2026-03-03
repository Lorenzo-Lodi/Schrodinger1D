package schrodinger.grid;

import org.junit.jupiter.api.Test;

import java.util.Arrays;

public class MappingStrategyTest {

    @Test
    public void checkGrid() {
        double rRef = 2.5;
        MappingStrategy s = new MappingSqrt(rRef);

        int rMaxFactor = 15;

        double rmin = 0.;
        double rmax = rMaxFactor * rRef;

        double ymin = s.y(rmin);
        double ymax = s.y(rmax);


        int np = 100;
        double h = (ymax - ymin) / (np - 1);

        double[] bucket = new double[rMaxFactor + 1];
        Arrays.fill(bucket, 0.0d);

        System.out.println("           y                    r");
        for (int i = 0; i < np; i++) {
            double y = ymin + h * i;
            int k = (int) (s.r(y) / rRef);
            bucket[k]++;
            String msg = String.format("%20.10f %20.10f  %d", y, s.r(y), k);
            System.out.println(msg);
        }

        for (int i = 0; i <= rMaxFactor; i++) {
            System.out.println(i + "  " + bucket[i]/np);
        }


    }

}
