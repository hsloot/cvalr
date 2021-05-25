#+ r load-package
library(cvalr)

#+ r config-and-setup
pd_args <- list(method = "CalibrationParam", sim_args = list(n_sim = 1e4))
times <- seq(0.25, 5, by = 0.25)
discount_factors <- rep(1, length(times))
recovery_rate <- 4e-1
spread <- 6e-2

composition <- c(28, 24, 16, 4, 3)
dim <- sum(composition)
lambda <- 8e-2
tau1 <- 4e-1
tau2 <- c(3e-2, 6e-2)

parm_caextmo2f <- CuadrasAugeExtMO2FParam(dim = dim, lambda = lambda, tau = tau1)
parm_asextmo2f <- AlphaStableExtMO2FParam(dim = dim, lambda = lambda, tau = tau1)
parm_poextmo2f <- PoissonExtMO2FParam(dim = dim, lambda = lambda, tau = tau1)
parm_exextmo2f <- ExponentialExtMO2FParam(dim = dim, lambda = lambda, tau = tau1)

parm_extga2f <- ExtGaussian2FParam(dim = dim, lambda = lambda, tau = tau1)

parm_clextar2f <- ClaytonExtArch2FParam(dim = dim, lambda = lambda, tau = tau1)
parm_frextar2f <- FrankExtArch2FParam(dim = dim, lambda = lambda, tau = tau1)
parm_guextar2f <- GumbelExtArch2FParam(dim = dim, lambda = lambda, tau = tau1)
parm_joextar2f <- JoeExtArch2FParam(dim = dim, lambda = lambda, tau = tau1)


parm_cah2extmo3f <- CuadrasAugeH2ExtMO3FParam(
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
  CuadrasAugeExtMO2FParam = expected_pcds_equation(
    parm_caextmo2f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  AlphaStableExtMO2FParam = expected_pcds_equation(
    parm_asextmo2f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  PoissonExtMO2FParam = expected_pcds_equation(
    parm_poextmo2f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  ExponentialExtMO2FParam = expected_pcds_equation(
    parm_exextmo2f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  ExtGaussianMO2FParam = expected_pcds_equation(
    parm_extga2f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  ClaytonExtArch2FParam = expected_pcds_equation(
    parm_clextar2f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  FrankExtArch2FParam = expected_pcds_equation(
    parm_frextar2f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  GumbelExtArch2FParam = expected_pcds_equation(
    parm_guextar2f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  JoeExtArch2FParam = expected_pcds_equation(
    parm_joextar2f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  CuadrasAugeH2ExtMO3FParam = expected_pcds_equation(
    parm_cah2extmo3f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  AlphaStableH2ExtMO3FParam = expected_pcds_equation(
    parm_ash2extmo3f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  PoissonH2ExtMO3FParam = expected_pcds_equation(
    parm_poh2extmo3f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  ExponentialH2ExtMO3FParam = expected_pcds_equation(
    parm_exh2extmo3f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  H2ExtGaussianMO3FParam = expected_pcds_equation(
    parm_h2extga3f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  ClaytonH2ExtArch3FParam = expected_pcds_equation(
    parm_clh2extar3f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  FrankH2ExtArch3FParam = expected_pcds_equation(
    parm_frh2extar3f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  GumbelH2ExtArch3FParam = expected_pcds_equation(
    parm_guh2extar3f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  JoeH2ExtArch3FParam = expected_pcds_equation(
    parm_joh2extar3f, times, discount_factors, recovery_rate, spread, 0,
    method = "CalibrationParam", pd_args = pd_args),
  min_time = 1, min_iterations = 1L, check = FALSE
)
