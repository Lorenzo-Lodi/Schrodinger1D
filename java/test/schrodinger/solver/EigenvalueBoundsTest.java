package schrodinger.solver;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;


public class EigenvalueBoundsTest {

    private static final double EPS = 1e-10;

    // Let us use harmonic-oscillator energies as example
// Energy nodes
// 0.25   0   // below ground state
// 0.5    0   GROUND STATE
// 1.5    1   FIRST EXCITED
// 2.5    2   SECOND EXCITED
// 3.5    3   THIRD EXCITED
// 4.5    4
// 5.5    5
// 6.5    6
    @Test
    public void test01() {
        EigenvalueBounds bounds = new EigenvalueBounds(5);
        double newValue = 0.25;
        bounds.updateBounds(newValue, 0);
        for (int v = 0; v < 5; v++) {
            assertEquals(newValue, bounds.getLowerBound(v), EPS);
        }
        newValue = 0.4;
        bounds.updateBounds(newValue, 0);
        for (int v = 0; v < 5; v++) {
            assertEquals(newValue, bounds.getLowerBound(v));
        }

        newValue = 0.51;
        bounds.updateBounds(newValue, 1);
        for (int v = 1; v < 5; v++) {
            assertEquals(newValue, bounds.getLowerBound(v));
        }
        assertEquals(newValue, bounds.getUpperBound(0));

        newValue = 7.5;
        bounds.updateBounds(newValue, 7);
        for (int v = 1; v < 5; v++) {
            assertEquals(newValue, bounds.getUpperBound(v));
        }
        assertEquals(0.51, bounds.getUpperBound(0));


    }

}
