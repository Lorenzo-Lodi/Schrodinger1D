package schrodinger.grid;

public interface MappingStrategy {

    default double y(double r) {
        return r;
    }

    default double r(double y) {
        return y;
    }

    default double g(double y) {
        return 1.d;
    }

    default double F(double y) {
        return 0;
    }
}
