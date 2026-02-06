package net.tinvention.schrodinger.grid;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

class MappingSurkusTest {

    @Test
    void mappingFunctionYofR_check_is_zero_at_rref() {
        // Try at many values of rRef and alpha that gives zero at Rref
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                double rRef = 0.01 + i;
                double alpha = 0.01 + j;
                MappingSurkus mapping = new MappingSurkus(rRef, alpha);
                assertEquals(0, mapping.y(rRef));
            }
        }
    }

    @Test
    void mappingFunctionYofR_check_is_minus_one_at_zero() {
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                double rRef = 0.01 + i;
                double alpha = 0.01 + j;
                MappingSurkus mapping = new MappingSurkus(rRef, alpha);
                assertEquals(-1.d, mapping.y(0));
            }
        }
    }

    @Test
    void mappingFunctionYofR_check_selected_values() {
        double rRef = 1.5;
        double alpha = 2;
        MappingSurkus mapping = new MappingSurkus(rRef, alpha);
        assertEquals(0.6d, mapping.y(2 * rRef), 1e-15);
        assertEquals(-0.6d, mapping.y(0.5d * rRef), 1e-15);

        rRef = 0.7;
        alpha = 0.8;
        mapping = new MappingSurkus(rRef, alpha);
        assertEquals(0.2703662113749116d, mapping.y(2 * rRef), 1e-15);
        assertEquals(-0.2703662113749117d, mapping.y(0.5d * rRef), 1e-15);
    }

    @Test
    void mappingFunctionRofY_check_selected_values() {
        double rRef = 1.5;
        double alpha = 2;
        MappingSurkus mapping = new MappingSurkus(rRef, alpha);
        assertEquals(0.8660254037844386d, mapping.r(-0.5), 1e-15);
        assertEquals(2.598076211353316d, mapping.r(0.5), 1e-15);

        rRef = 0.7;
        alpha = 0.8;
        mapping = new MappingSurkus(rRef, alpha);
        assertEquals(0.1772949933187049d, mapping.r(-0.5), 1e-15);
        assertEquals(2.763755427200234d, mapping.r(0.5), 1e-15);
    }

    void mappingFunctionRofY_and_Y_of_R() {
        double rRef = 1.5;
        double alpha = 2;
        MappingSurkus mapping = new MappingSurkus(rRef, alpha);
        for (int i = 0; i < 10; i++) {
            double r = 0.01d + i / 0.5d;
            assertEquals(r, mapping.r(mapping.y(r)), 1e-15);
        }
    }

    @Test
    void mappingFunctionGofY_check_selected_values() {
        double rRef = 1.5;
        double alpha = 2;
        MappingSurkus mapping = new MappingSurkus(rRef, alpha);
        assertEquals(1.1547005383792515d, mapping.g(-0.5), 1e-15);
        assertEquals(3.4641016151377544d, mapping.g(0.5), 1e-15);
    }

    @Test
    void mappingFunctionGofY_check_is_first_derivative_of_R_of_Y() {
        double rRef = 1.5;
        double alpha = 2;
        MappingSurkus mapping = new MappingSurkus(rRef, alpha);

        double eps = 1e-6;
        int numberOfPoints = 20;
        double rmin = 0.1;
        double rmax = 10;
        double step = (rmax - rmin) / (numberOfPoints - 1);

        for (int i = 0; i < numberOfPoints; i++) {
            double y = 0.1 + step * i;
            double der1 = (mapping.r(y + eps) - mapping.r(y - eps)) / (2.d * eps);
            assertEquals(der1, mapping.g(y), 1e-8);
        }
    }

    @Test
    void mappingFunctionFofY_check_selected_values() {
        double rRef = 1.5;
        double alpha = 2;
        MappingSurkus mapping = new MappingSurkus(rRef, alpha);
        assertEquals(1.3333333333333333d, mapping.F(-0.5), 1e-15);
        assertEquals(1.3333333333333333d, mapping.F(0.5), 1e-15);
    }

}
