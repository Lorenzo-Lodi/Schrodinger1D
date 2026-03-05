package schrodinger.grid;

import org.junit.jupiter.api.Test;

class MappingSqrtTest extends MappingAbstractTest {

    @Test
    public void mappingFunctionGofY_check_is_first_derivative_of_R_of_Y() {
        double rRef = 1.5;
        MappingStrategy mapping = new MappingSqrt(rRef);
        super.mappingFunctionGofY_check_is_first_derivative_of_R_of_Y(mapping);
    }

}
