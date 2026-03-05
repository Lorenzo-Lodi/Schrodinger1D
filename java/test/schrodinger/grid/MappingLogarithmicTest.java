package schrodinger.grid;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class MappingLogarithmicTest extends MappingAbstractTest {

    @Test
    public void mappingFunctionGofY_check_is_first_derivative_of_R_of_Y() {
        double rRef = 1.5;
        MappingStrategy mapping = new MappingLogarithmic(rRef);
        super.mappingFunctionGofY_check_is_first_derivative_of_R_of_Y(mapping);
    }

}
