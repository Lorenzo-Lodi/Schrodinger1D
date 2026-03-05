package schrodinger.grid;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class MappingAbstractTest {


    void mappingFunctionGofY_check_is_first_derivative_of_R_of_Y(MappingStrategy mapping) {
        double eps = 1e-6;
        int numberOfPoints = 23;
        double rmin = 0.;
        double rmax = 10.;
        double step = (rmax - rmin) / (numberOfPoints - 1);

        for (int i = 0; i < numberOfPoints; i++) {
            double y = mapping.y(0.1 + step * i);
            double der1 = (mapping.r(y + eps) - mapping.r(y - eps)) / (2.d * eps);
            assertEquals(der1, mapping.g(y), 1e-8 * der1);
        }
    }

}
