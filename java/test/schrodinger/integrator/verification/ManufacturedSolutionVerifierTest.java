package schrodinger.integrator.verification;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class ManufacturedSolutionVerifierTest {

    @Test
    void test_gridParamsToArray_basic_test01() {
        int[] points = ManufacturedSolutionVerifier.gridParamsToArray(100, 200, 100);
        assertEquals(2, points.length);
        assertEquals(100, points[0]);
        assertEquals(200, points[1]);
    }

    @Test
    void test_gridParamsToArray_basic_test02() {
        int[] points = ManufacturedSolutionVerifier.gridParamsToArray(100, 220, 100);
        assertEquals(2, points.length);
        assertEquals(100, points[0]);
        assertEquals(200, points[1]);
    }

    @Test
    void test_gridParamsToArray_basic_test03() {
        int[] points = ManufacturedSolutionVerifier.gridParamsToArray(100, 100, 100);
        assertEquals(1, points.length);
        assertEquals(100, points[0]);
    }

    @Test
    void test_gridParamsToArray_basic_test04() {
        int[] points = ManufacturedSolutionVerifier.gridParamsToArray(100, 100, 1);
        assertEquals(1, points.length);
        assertEquals(100, points[0]);
    }

    @Test
    void test_gridParamsToArray_basic_test05() {
        int[] points = ManufacturedSolutionVerifier.gridParamsToArray(100, 1234, 75);
        assertEquals(16, points.length);
        assertEquals(100, points[0]);
        assertEquals(175, points[1]);
        assertEquals(1225, points[15]);
    }


}