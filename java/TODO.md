1. Add more logs (writes) thoughout, especially to write the level-specific values
2. Investigate while in DissociationTest sometimes we have NaN using secant.
3. Investigate while in DissociationTest sometimes we have convergence to wrong state using secant.
4. Investigate while in DissociationTest regula falsi often fails
5. Improve the code to evaluate the energy scale; potentially move to SchrodingerSystem
6. Write boostrapping for methods that need more than 2 points. Richardson-extrapolated numerov (or EFN) with finer, standard, coarse grids should be enough. Order, in theory becomes 8.
7. Write more integration tests, for a combination of states, grids, integrators.
8. Make integrator pass in: 1. An half-integer index (allow for that)
9. Make integrator pass in: 2. Derivative(s) of QTilde(y)
10. Make integrator pass in: 3. psiPrime[i]
11. Implement some Runge-Kutta-Nystrom methods (just for fun)
12. Implement numerov using derivatives
13. Implement general polynomial potential
14. Make the integrator pre-compute the points for potential/ grid transformation
15. Implement new grid mappings linear-switchover-linear, using arctan and using rational function
