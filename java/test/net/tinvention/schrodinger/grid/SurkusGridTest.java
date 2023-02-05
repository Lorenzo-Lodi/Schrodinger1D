package net.tinvention.schrodinger.grid;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class SurkusGridTest {

	@Test
	void mappingFunctionYofR_check_is_zero_at_rref() {
		// Try at many values of rRef and alpha that gives zero at Rref
		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 10; j++) {
				double rRef = 0.01 + i;
				double alpha = 0.01 + j;
				SurkusGrid grid = new SurkusGrid(0.3, 5.5, 10, rRef, alpha);
				assertEquals(0, grid.mappingFunctionYofR(rRef));
			}
		}
	}

	@Test
	void mappingFunctionYofR_check_is_minus_one_at_zero() {
		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 10; j++) {
				double rRef = 0.01 + i;
				double alpha = 0.01 + j;
				SurkusGrid grid = new SurkusGrid(0.3, 5.5, 10, rRef, alpha);
				assertEquals(-1.d, grid.mappingFunctionYofR(0));
			}
		}
	}

	@Test
	void mappingFunctionYofR_check_selected_values() {
		double rRef = 1.5;
		double alpha = 2;
		SurkusGrid grid = new SurkusGrid(0.3, 5.5, 10, rRef, alpha);
		assertEquals(0.6d, grid.mappingFunctionYofR(2 * rRef), 1e-15);
		assertEquals(-0.6d, grid.mappingFunctionYofR(0.5d * rRef), 1e-15);

		rRef = 0.7;
		alpha = 0.8;
		grid = new SurkusGrid(0.3, 5.5, 10, rRef, alpha);
		assertEquals(0.2703662113749116d, grid.mappingFunctionYofR(2 * rRef), 1e-15);
		assertEquals(-0.2703662113749117d, grid.mappingFunctionYofR(0.5d * rRef), 1e-15);
	}

	@Test
	void mappingFunctionRofY_check_selected_values() {
		double rRef = 1.5;
		double alpha = 2;
		SurkusGrid grid = new SurkusGrid(0.3, 5.5, 10, rRef, alpha);
		assertEquals(0.8660254037844386d, grid.mappingFunctionRofY(-0.5), 1e-15);
		assertEquals(2.598076211353316d, grid.mappingFunctionRofY(0.5), 1e-15);

		rRef = 0.7;
		alpha = 0.8;
		grid = new SurkusGrid(0.3, 5.5, 10, rRef, alpha);
		assertEquals(0.1772949933187049d, grid.mappingFunctionRofY(-0.5), 1e-15);
		assertEquals(2.763755427200234d, grid.mappingFunctionRofY(0.5), 1e-15);
	}

	void mappingFunctionRofY_and_Y_of_R() {
		double rRef = 1.5;
		double alpha = 2;
		SurkusGrid grid = new SurkusGrid(0.3, 5.5, 10, rRef, alpha);
		for (int i = 0; i < 10; i++) {
			double r = 0.01d + i / 0.5d;
			assertEquals(r, grid.mappingFunctionRofY(grid.mappingFunctionYofR(r)), 1e-15);
		}
	}

	@Test
	void mappingFunctionGofY_check_selected_values() {
		double rRef = 1.5;
		double alpha = 2;
		SurkusGrid grid = new SurkusGrid(0.3, 5.5, 10, rRef, alpha);
		assertEquals(1.1547005383792515d, grid.mappingFunctionGofY(-0.5), 1e-15);
		assertEquals(3.4641016151377544d, grid.mappingFunctionGofY(0.5), 1e-15);
	}

	@Test
	void mappingFunctionFofY_check_selected_values() {
		double rRef = 1.5;
		double alpha = 2;
		SurkusGrid grid = new SurkusGrid(0.3, 5.5, 10, rRef, alpha);
		assertEquals(1.3333333333333333d, grid.mappingFunctionFofY(-0.5), 1e-15);
		assertEquals(1.3333333333333333d, grid.mappingFunctionFofY(0.5), 1e-15);
	}

	
}
