package schrodinger;

import static org.junit.jupiter.api.Assertions.assertEquals;

import schrodinger.grid.GridFactory;
import org.junit.jupiter.api.Test;

import schrodinger.grid.Grid;
import schrodinger.potential.PhysicalPotential;

class QuantumLevelTest {
    private static final double SMALL_EPS = 1e-15;
    private static final double LARGE_EPS = 8e-4;

    private final Grid grid;
    private final SchrodingerSystem system;
    private final QuantumLevel level;

    public QuantumLevelTest() {
        grid = GridFactory.generateUniformGrid(0, 1, 50);
        system = new SchrodingerSystem(x -> 0, 1, grid);
        level = new QuantumLevel(system);
    }

    @Test
    void normalizePsiTrapezoidalRuleConstantIntegrand() {
        for (int i = 0; i < grid.getNumberOfPoints(); i++) {
            level.psi[i] = 1;
        }
        double normalizationFactor = level.normalizePsi();
        assertEquals(1.0d, normalizationFactor, SMALL_EPS);
    }

    @Test
    void normalizePsiTrapezoidalRuleLinearIntegrand() {
        for (int i = 0; i < grid.getNumberOfPoints(); i++) {
            level.psi[i] = grid.rAtGridPoint(i);
        }
        double normalizationFactor = level.normalizePsi();
        assertEquals(1.73187048891682771453, normalizationFactor, SMALL_EPS);
    }

    @Test
    void normalizePsiTrapezoidalRuleQuadraticIntegrand() {
        for (int i = 0; i < grid.getNumberOfPoints(); i++) {
            level.psi[i] = Math.pow(grid.rAtGridPoint(i), 2);
        }
        double normalizationFactor = level.normalizePsi();
        assertEquals(2.2352923244791710742, normalizationFactor, LARGE_EPS);
    }

    @Test
    void normalizePsiTrapezoidalRuleCubicIntegrand() {
        for (int i = 0; i < grid.getNumberOfPoints(); i++) {
            level.psi[i] = Math.pow(grid.rAtGridPoint(i), 3);
        }
        double normalizationFactor = level.normalizePsi();
        assertEquals(2.6438252937568609246, normalizationFactor, SMALL_EPS);
    }

    @Test
    void normalizePsiTrapezoidalRuleQuarticIntegrand() {
        for (int i = 0; i < grid.getNumberOfPoints(); i++) {
            level.psi[i] = Math.pow(grid.rAtGridPoint(i), 4);
        }
        double normalizationFactor = level.normalizePsi();
        assertEquals(2.9962596611853944811, normalizationFactor, LARGE_EPS);
    }

}
