package net.tinvention.schrodinger;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import net.tinvention.schrodinger.grid.UniformGrid;
import net.tinvention.schrodinger.grid.Grid;

class EnergyLevelTest {
	private static final double SMALL_EPS = 5e-16;
	private static final double LARGE_EPS = 8e-4;

	@Test
	void normalizePsiTrapezoidalRuleConstantIntegrand() {
		Grid grid = new UniformGrid(0, 1, 50);
		EnergyLevel level = new EnergyLevel(grid);
		for (int i = 0; i < grid.getNumberOfPoints(); i++) {
			level.psi[i] = 1;
		}
		double normalizationFactor = level.normalizePsi();
		assertEquals(1.0d, normalizationFactor, SMALL_EPS);
	}

	@Test
	void normalizePsiTrapezoidalRuleLinearIntegrand() {
		Grid grid = new UniformGrid(0, 1, 50);
		EnergyLevel level = new EnergyLevel(grid);
		for (int i = 0; i < grid.getNumberOfPoints(); i++) {
			level.psi[i] = grid.getRValue(i);
		}
		double normalizationFactor = level.normalizePsi();
		assertEquals(1.4142135623730950488, normalizationFactor, SMALL_EPS);
	}

	@Test
	void normalizePsiTrapezoidalRuleQuadraticIntegrand() {
		Grid grid = new UniformGrid(0, 1, 50);
		EnergyLevel level = new EnergyLevel(grid);
		for (int i = 0; i < grid.getNumberOfPoints(); i++) {
			level.psi[i] = Math.pow(grid.getRValue(i), 2);
		}
		double normalizationFactor = level.normalizePsi();
		assertEquals(1.7320508075688772935, normalizationFactor, LARGE_EPS);
	}

	@Test
	void normalizePsiTrapezoidalRuleCubicIntegrand() {
		Grid grid = new UniformGrid(0, 1, 50);
		EnergyLevel level = new EnergyLevel(grid);
		for (int i = 0; i < grid.getNumberOfPoints(); i++) {
			level.psi[i] = Math.pow(grid.getRValue(i), 3);
		}
		double normalizationFactor = level.normalizePsi();
		assertEquals(2.0d, normalizationFactor, LARGE_EPS);
	}

	@Test
	void normalizePsiTrapezoidalRuleQuarticIntegrand() {
		Grid grid = new UniformGrid(0, 1, 50);
		EnergyLevel level = new EnergyLevel(grid);
		for (int i = 0; i < grid.getNumberOfPoints(); i++) {
			level.psi[i] = Math.pow(grid.getRValue(i), 4);
		}
		double normalizationFactor = level.normalizePsi();
		assertEquals(2.2360679774997896964, normalizationFactor, LARGE_EPS);
	}

}
