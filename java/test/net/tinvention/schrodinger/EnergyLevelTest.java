package net.tinvention.schrodinger;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import net.tinvention.schrodinger.grid.UniformGrid;
import net.tinvention.schrodinger.grid.Grid;

class EnergyLevelTest {

	@Test
	void normalizePsiTrapezoidalRule() {
		Grid grid = new UniformGrid(0, 2, 50);
		EnergyLevel level = new EnergyLevel(grid);
		for (int i = 0; i < grid.getNumberOfPoints(); i++) {
			level.psi[i] = 1;
		}
		double normalizationFactor = level.normalizePsiTrapezoidalRule();
		assertEquals(0.70710678118654752440, normalizationFactor, 1e-16);
		
		for (int i = 0; i < grid.getNumberOfPoints(); i++) {
			level.psi[i] = i;
		}
		normalizationFactor = level.normalizePsiTrapezoidalRule();
		assertEquals(0.70710678118654752440, normalizationFactor, 1e-16);


	}

}
