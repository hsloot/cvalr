#+ r load-package
library(cvalr)
options("cvalr.enable_messages" = FALSE, "cvalr.enable_warnings" = FALSE)

#+ r config-and-setup
times <- seq(25e-2, 5, by = 25e-2)
discount_factors <- rep(1, length(times))
recovery_rate <- 4e-1
lower <- c(0, 0.1, 0.2, 0.35)
upper <- c(0.1, 0.2, 0.35, 1)
coupon <- c(rep(5e-2, 3), 0)
upfront <- c(8e-1, 5e-1, 1e-1, 0)

composition <- c(28, 24, 16, 4, 3)
dim <- sum(composition)
lambda <- 8e-2
tau1 <- 4e-1
tau2 <- c(3e-2, 6e-2)

parm_caextmo2f <- ArmageddonExtMO2FParam(dim = dim, lambda = lambda, tau = tau1)
parm_asextmo2f <- AlphaStableExtMO2FParam(dim = dim, lambda = lambda, tau = tau1)
parm_poextmo2f <- PoissonExtMO2FParam(dim = dim, lambda = lambda, tau = tau1)
parm_exextmo2f <- ExponentialExtMO2FParam(dim = dim, lambda = lambda, tau = tau1)

parm_extga2f <- ExtGaussian2FParam(dim = dim, lambda = lambda, tau = tau1)

parm_clextar2f <- ClaytonExtArch2FParam(dim = dim, lambda = lambda, tau = tau1)
parm_frextar2f <- FrankExtArch2FParam(dim = dim, lambda = lambda, tau = tau1)
parm_guextar2f <- GumbelExtArch2FParam(dim = dim, lambda = lambda, tau = tau1)
parm_joextar2f <- JoeExtArch2FParam(dim = dim, lambda = lambda, tau = tau1)


parm_cah2extmo3f <- ArmageddonH2ExtMO3FParam(
  composition = composition, lambda = lambda, tau = tau2)
parm_ash2extmo3f <- AlphaStableH2ExtMO3FParam(
  composition = composition, lambda = lambda, tau = tau2)
parm_poh2extmo3f <- PoissonH2ExtMO3FParam(
  composition = composition, lambda = lambda, tau = tau2)
parm_exh2extmo3f <- ExponentialH2ExtMO3FParam(
  composition = composition, lambda = lambda, tau = tau2)

parm_h2extga3f <- H2ExtGaussian3FParam(composition = composition, lambda = lambda, tau = tau2)

parm_clh2extar3f <- ClaytonH2ExtArch3FParam(composition = composition, lambda = lambda, tau = tau2)
parm_frh2extar3f <- FrankH2ExtArch3FParam(composition = composition, lambda = lambda, tau = tau2)
parm_guh2extar3f <- GumbelH2ExtArch3FParam(composition = composition, lambda = lambda, tau = tau2)
parm_joh2extar3f <- JoeH2ExtArch3FParam(composition = composition, lambda = lambda, tau = tau2)

bench::mark(
  ArmageddonExtMO2FParam = expected_cdo_equation(
    parm_caextmo2f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  AlphaStableExtMO2FParam = expected_cdo_equation(
    parm_asextmo2f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  PoissonExtMO2FParam = expected_cdo_equation(
    parm_poextmo2f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  ExponentialExtMO2FParam = expected_cdo_equation(
    parm_exextmo2f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  ExtGaussian2FParam = expected_cdo_equation(
    parm_extga2f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  ClaytonExtArch2FParam = expected_cdo_equation(
    parm_clextar2f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  FrankExtArch2FParam = expected_cdo_equation(
    parm_frextar2f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  GumbelExtArch2FParam = expected_cdo_equation(
    parm_guextar2f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  JoeExtArch2FParam = expected_cdo_equation(
    parm_joextar2f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  ArmageddonH2ExtMO3FParam = expected_cdo_equation(
    parm_cah2extmo3f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  AlphaStableH2ExtMO3FParam = expected_cdo_equation(
    parm_ash2extmo3f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  PoissonH2ExtMO3FParam = expected_cdo_equation(
    parm_poh2extmo3f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  ExponentialH2ExtMO3FParam = expected_cdo_equation(
    parm_exh2extmo3f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  H2ExtGaussian3FParam = expected_cdo_equation(
    parm_h2extga3f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  ClaytonH2ExtArch3FParam = expected_cdo_equation(
    parm_clh2extar3f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  FrankH2ExtArch3FParam = expected_cdo_equation(
    parm_frh2extar3f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  GumbelH2ExtArch3FParam = expected_cdo_equation(
    parm_guh2extar3f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  JoeH2ExtArch3FParam = expected_cdo_equation(
    parm_joh2extar3f, times, discount_factors, recovery_rate, lower, upper, coupon, upfront),
  min_time = 1, min_iterations = 1L, check = FALSE
)
