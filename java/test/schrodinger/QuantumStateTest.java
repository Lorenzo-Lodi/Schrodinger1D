package schrodinger;

import static org.junit.jupiter.api.Assertions.assertEquals;

import schrodinger.grid.GridFactory;
import org.junit.jupiter.api.Test;

import schrodinger.grid.Grid;
import schrodinger.potential.SchrodingerSystem;

class QuantumStateTest {
    private static final double SMALL_EPS = 1e-15;
    private static final double LARGE_EPS = 8e-4;
    private static final SchrodingerSystem system = new SchrodingerSystem(null, 1, null);

    @Test
    void normalizePsiTrapezoidalRuleConstantIntegrand() {
        Grid grid = GridFactory.generateUniformGrid(0, 1, 50);
        QuantumState level = new QuantumState(system, grid);
        for (int i = 0; i < grid.getNumberOfPoints(); i++) {
            level.psi[i] = 1;
        }
        double normalizationFactor = level.normalizePsi();
        assertEquals(1.0d, normalizationFactor, SMALL_EPS);
    }

    @Test
    void normalizePsiTrapezoidalRuleLinearIntegrand() {
        Grid grid = GridFactory.generateUniformGrid(0, 1, 50);
        QuantumState level = new QuantumState(system, grid);
        for (int i = 0; i < grid.getNumberOfPoints(); i++) {
            level.psi[i] = grid.getRValue(i);
        }
        double normalizationFactor = level.normalizePsi();
        assertEquals(1.73187048891682771453, normalizationFactor, SMALL_EPS);
    }

    @Test
    void normalizePsiTrapezoidalRuleQuadraticIntegrand() {
        Grid grid = GridFactory.generateUniformGrid(0, 1, 50);
        QuantumState level = new QuantumState(system, grid);
        for (int i = 0; i < grid.getNumberOfPoints(); i++) {
            level.psi[i] = Math.pow(grid.getRValue(i), 2);
        }
        double normalizationFactor = level.normalizePsi();
        assertEquals(2.2352923244791710742, normalizationFactor, LARGE_EPS);
    }

    @Test
    void normalizePsiTrapezoidalRuleCubicIntegrand() {
        Grid grid = GridFactory.generateUniformGrid(0, 1, 50);
        QuantumState level = new QuantumState(system, grid);
        for (int i = 0; i < grid.getNumberOfPoints(); i++) {
            level.psi[i] = Math.pow(grid.getRValue(i), 3);
        }
        double normalizationFactor = level.normalizePsi();
        assertEquals(2.6438252937568609246, normalizationFactor, SMALL_EPS);
    }

    @Test
    void normalizePsiTrapezoidalRuleQuarticIntegrand() {
        Grid grid = GridFactory.generateUniformGrid(0, 1, 50);
        QuantumState level = new QuantumState(system, grid);
        for (int i = 0; i < grid.getNumberOfPoints(); i++) {
            level.psi[i] = Math.pow(grid.getRValue(i), 4);
        }
        double normalizationFactor = level.normalizePsi();
        assertEquals(2.9962596611853944811, normalizationFactor, LARGE_EPS);
    }

}
