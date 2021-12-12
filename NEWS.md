# cvalr 0.5.2

- Realign with `rmo@v0.8` because of renamed algorithms

# cvalr 0.5.1

- Realign with `rmo@v0.7.1` because of renamed algorithms

# cvalr 0.5

- Realign with `rmo@v0.7.0` because of breaking changes in the `PoissonBernsteinFunction`-class, see hsloot/rmo#81.

# cvalr 0.4

- Rename `"CuadrasAuge*"` classes to `"Armageddon*"` because of `rmo@v0.6`.
- Add additional tests and refactor tests for simulating calibration parameters.

# cvalr 0.3.3

- Refactor backend using different, linear representations.

# cvalr 0.3.2

- Change logic of `attrs` argument for `expected_value`.

# cvalr 0.3.1

- Fix H2-ext. CA-sampling and problems with Miwa (`copula` and `mvtnorm` algorithm)

# cvalr 0.3

- Internal refactoring of `expected_value` methods and user interface.
- Improved sampling of *average default counting process* for ext. MO families
- Implement missing `probability_distribution` method for `H2ExMarkovParam`-class
- Implement missing `expected_pcds_loss` methods
- Include continuous benchmark integration on Github
- Fix validity methods

# cvalr 0.2

- Code refactoring and tests added

# cvalr 0.1

- Calibration parameter classes (extendible and H2-extendible)
- Auxiliary methods for portfolio CDS and CDO pricing
