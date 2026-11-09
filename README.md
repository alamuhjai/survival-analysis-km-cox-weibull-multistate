# Survival Analysis: Kaplan–Meier, Cox, Weibull and Multi-State Models in R

This repository holds three R demos of time-to-event methods. The first two are single-event survival analyses, and the third is a multi-state model.

- **Kaplan-Meier, RMST and stratified Cox regression (VA lung cancer trial):** non-parametric survival curves, 1-year restricted mean survival time, log-rank and Peto-Peto tests, and a stratified Cox model with a treatment-by-prior-therapy interaction (`R/nonparam.R`, `R/coxDiag.R`).
- **Weibull PH vs. AFT and Poisson regression (NCCTG lung cancer):** one Weibull model read in both proportional-hazards and accelerated-failure-time parameterisations, compared against a constant-hazard Poisson GLM with a log person-time offset (`R/weibullAft.R`).
- **Reversible illness-death model (synthetic person-year cohort):** a three-state continuous-time model with annual transition rates, a likelihood-ratio test for time-homogeneity, and 3-year sojourn and occupation times (`R/multistate.R`).

---

## Kaplan-Meier, RMST and stratified Cox regression (VA lung cancer trial)

The `survival::veteran` dataset comes from Kalbfleisch and Prentice (1980) and follows n = 137 patients with advanced, inoperable lung cancer. During follow-up, 128 patients died and 9 were right-censored. Patients were randomised to standard therapy (`trt = 1`, n = 69, 64 deaths) or test chemotherapy (`trt = 2`, n = 68, 64 deaths). Baseline covariates are histological cell type (`squamous`, `smallcell`, `adeno`, `large`), Karnofsky performance score (`karno`, split at `karno <= 50` vs. `karno > 50` in `karno1`), and prior therapy (`prior`, coded `0` for no and `10` for yes).

### Non-parametric survival and 1-year RMST

I fit Kaplan-Meier survival curves $\hat{S}(t)$ and Nelson-Aalen cumulative hazard curves $\hat{\Lambda}(t)$ with `survfit` using `ctype = 1` and complementary log-log 95% confidence intervals (`conf.type = "log-log"`). Because survival curves cross across cell types and violate proportional hazards, I report the 1-year restricted mean survival time up to $\tau = 365$ days next to the median:

$$
\mathrm{RMST}(\tau) = \int_0^{\tau} S(t)\mathrm{d}t.
$$

I test group differences with `survdiff` using both the unweighted log-rank test (`rho = 0`) and the Peto-Peto test (`rho = 1`), which weights early event times by $\hat{S}(t^-)$.

| Group | Stratum | n | Events | Median (log-log 95% CI, days) | 1-yr RMST (SE, days) | RMST Wald 95% CI | Log-rank p (`rho = 0`) | Peto-Peto p (`rho = 1`) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `celltype` | `squamous` | 35 | 31 | 118.0 (44, 242) | 170.64 (24.23) | (123.14, 218.14) | 1.27e-05 | 1.95e-04 |
| `celltype` | `smallcell` | 48 | 45 | 51.0 (24, 61) | 77.64 (13.97) | (50.25, 105.02) | 1.27e-05 | 1.95e-04 |
| `celltype` | `adeno` | 27 | 26 | 51.0 (24, 90) | 65.56 (9.93) | (46.09, 85.02) | 1.27e-05 | 1.95e-04 |
| `celltype` | `large` | 27 | 26 | 156.0 (100, 216) | 162.23 (20.46) | (122.13, 202.34) | 1.27e-05 | 1.95e-04 |
| `karno1` | `karno <= 50` | 52 | 50 | 20.5 (16, 33) | 61.86 (12.63) | (37.11, 86.60) | 1.08e-07 | 2.98e-11 |
| `karno1` | `karno > 50` | 85 | 78 | 112.0 (95, 139) | 147.97 (12.62) | (123.24, 172.70) | 1.08e-07 | 2.98e-11 |
| `trt` | `Standard` | 69 | 64 | 103.0 (54, 126) | 118.97 (13.02) | (93.45, 144.49) | 0.928 | 0.351 |
| `trt` | `Test` | 68 | 64 | 52.5 (43, 90) | 112.40 (14.87) | (83.25, 141.56) | 0.928 | 0.351 |

Survival differs strongly by cell type (χ² = 25.40 on 3 df, p = 1.27e-05) and baseline Karnofsky status (log-rank χ² = 28.22 on 1 df, p = 1.08e-07; Peto-Peto χ² = 44.19 on 1 df, p = 2.98e-11). Over the first year, patients with `karno > 50` live 147.97 days on average compared with 61.86 days for `karno <= 50`, a difference of 86.11 days. The unadjusted treatment comparison shows no overall difference (log-rank χ² = 0.01 on 1 df, p = 0.928).

![Kaplan-Meier Survival and Nelson-Aalen Cumulative Hazard by Cell Type](figures/kmCurve.png)

### Stratified Cox regression and residual diagnostics

In an unstratified Cox model (`coxUnstrat`), both the complementary log-log curves and the Schoenfeld residual test (`cox.zph`) show non-proportional hazards for `celltype` (χ² = 12.591 on 3 df, p = 0.0056) and `karno1` (χ² = 15.788 on 1 df, p = 7.1e-05). The global Schoenfeld test is χ² = 27.974 on 6 df (p = 9.5e-05).

I therefore stratify on the 8 combinations of `celltype` and `karno1`, indexed by $s = 1, \dots, 8$ with stratum-specific baseline hazard $\lambda _ {0s}(t)$ (per day), and test whether prior therapy modifies the treatment effect (`fitCox`):

$$
\lambda_s(t \mid x) = \lambda_{0s}(t)\exp\left(\beta_1 x_{\mathrm{trt}} + \beta_2 x_{\mathrm{prior}} + \beta_3 x_{\mathrm{trt}} x_{\mathrm{prior}}\right).
$$

Here $x _ {\mathrm{trt}}$ (`trt`, coded `1` for standard and `2` for test) increases by 1 unit between arms, while $x _ {\mathrm{prior}}$ (`prior`) is coded `0` vs. `10`. Because $\beta_3$ is per unit of `prior`, the log hazard ratio for test vs. standard therapy is $\beta_1$ when `prior = 0` and $\beta_1 + 10\beta_3$ when `prior = 10`. I compute the variance of the pre-treated contrast from `vcov(fitCox)`:

$$
\mathrm{Var}\left(\hat{\beta}_1 + 10\hat{\beta}_3\right) = \mathrm{Var}(\hat{\beta}_1) + 100\mathrm{Var}(\hat{\beta}_3) + 20\mathrm{Cov}(\hat{\beta}_1, \hat{\beta}_3).
$$

| Model | Contrast | $\log \mathrm{HR}$ | SE | $\mathrm{HR}$ | Wald 95% CI | Wald p | Schoenfeld p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Stratified, main effects | Test vs. Standard (`trt`) | 0.1189 | 0.2037 | 1.1262 | (0.7555, 1.6789) | 0.5595 | 0.5257 |
| Stratified, main effects | Prior yes vs. no (`prior`, 10 vs. 0) | 0.1024 | 0.2156 | 1.1078 | (0.7259, 1.6905) | 0.6350 | 0.2258 |
| Stratified, interaction | Test vs. Standard, no prior therapy (`prior = 0`) | 0.3036 | 0.2396 | 1.3548 | (0.8470, 2.1670) | 0.2051 | - |
| Stratified, interaction | Test vs. Standard, prior therapy (`prior = 10`) | -0.3153 | 0.3628 | 0.7296 | (0.3583, 1.4854) | 0.3848 | - |
| Stratified, interaction | Interaction: ratio of HRs (`prior = 10` vs. `0`) | -0.6189 | 0.4266 | 0.5385 | (0.2334, 1.2426) | 0.1468 | - |

After stratifying on `celltype` and `karno1`, the global Schoenfeld test shows no remaining violation (χ² = 1.846 on 2 df, p = 0.40). Among patients with no prior therapy, test chemotherapy has HR = 1.3548 (Wald 95% CI: 0.8470, 2.1670; p = 0.2051). Among previously treated patients, test chemotherapy has HR = 0.7296 (Wald 95% CI: 0.3583, 1.4854; p = 0.3848). The interaction ratio is $\exp(10\hat{\beta} _ 3)$ = 0.5385 (Wald 95% CI: 0.2334, 1.2426; p = 0.1468).

![Complementary Log-Log Curves and Scaled Schoenfeld Residuals](figures/phDiag.png)

I plot null-model Martingale residuals against continuous `karno` to check the split at `karno = 50`, and `dfbeta` residuals to check case influence on the treatment slope. First-order `dfbeta` residuals shift $\hat{\beta} _ {\mathrm{trt}}$ by at most 0.046, less than a quarter of its standard error (0.2037). An exact leave-one-out refit gives a maximum shift of 0.068, about a third of the standard error.

![Martingale Functional Form and dfbeta Influence Diagnostics](figures/coxResid.png)

---

## Weibull PH vs. AFT and Poisson regression (NCCTG lung cancer)

The `survival::lung` dataset from the North Central Cancer Treatment Group (Loprinzi et al., 1994) has 228 patients in the raw table. Dropping 1 patient with a missing `ph.ecog` leaves n = 227 complete cases with 164 deaths and 63 right-censored observations over 69,522 person-days. I model survival time $T$ (in days) using `age` (in years), `sex` (`1` = male, `2` = female), and ECOG performance score (`ph.ecog`).

### Connecting Weibull PH, Weibull AFT and Poisson regression

The Weibull distribution is both a proportional-hazards (PH) model and an accelerated-failure-time (AFT) model. With shape $\gamma > 0$, baseline scale $\lambda_0 > 0$, and covariate vector $x_i$,

$$
\begin{aligned}
\text{PH:}\quad \lambda(t \mid x_i) &= \gamma\lambda_0 t^{\gamma - 1}\exp\left(x_i^{\top}\beta_{\mathrm{PH}}\right), \\
\text{AFT:}\quad \log T_i &= \mu + x_i^{\top}\beta_{\mathrm{AFT}} + \sigma W_i, \qquad \sigma = \frac{1}{\gamma},
\end{aligned}
$$

where $\mu = -\log(\lambda_0)/\gamma$ and $W_i$ follows the standard minimum extreme-value distribution. I fit the AFT form with `survreg` (`weib2`) and convert to the PH parameterisation (`weib1`) using

$$
\beta_{\mathrm{AFT}} = -\frac{\beta_{\mathrm{PH}}}{\gamma}.
$$

A negative $\beta _ {\mathrm{PH}}$ (lower death hazard, $\mathrm{HR} = \exp(\beta _ {\mathrm{PH}}) < 1$) corresponds to a positive $\beta _ {\mathrm{AFT}}$ (longer survival time, $\exp(\beta _ {\mathrm{AFT}}) > 1$).

When $\gamma = 1$, the hazard is constant. With event indicator $d_i \in \{0, 1\}$, follow-up time $t_i$ (in days), and rate $\lambda_i = \exp(\alpha + x_i^{\top}\beta _ {\mathrm{PH}})$, the exponential survival log-likelihood is

$$
\ell(\alpha, \beta_{\mathrm{PH}}) = \sum_{i=1}^{n}\left(d_i\log\lambda_i - \lambda_i t_i\right).
$$

This has the same kernel as a Poisson GLM for $d_i \sim \mathrm{Poisson}(\mu_i)$ with mean $\mu_i = t_i\lambda_i$ and offset $\log t_i$, which I fit with `glm(event ~ age + sex + ph.ecog, offset = log(time), family = poisson)` (`model2`). Because the Poisson model fixes $\gamma = 1$, its hazard ratios differ slightly from the Weibull PH hazard ratios where $\hat{\gamma}$ = 1.3678.

For the homogeneous constant hazard $\hat{\lambda} = D / \sum_i t_i$ (164 deaths over 69,522 person-days, or 2.3590 deaths per 1,000 person-days), I also compare the log-scale delta-method 95% Wald interval $\hat{\lambda}\exp(\pm 1.96/\sqrt{D})$ = (2.0242, 2.7491) against the linear Wald interval (1.9979, 2.7200) per 1,000 person-days.

| Covariate | $\hat{\beta} _ {\mathrm{PH}}$ (SE) | PH $\widehat{\mathrm{HR}}$ (Wald 95% CI) | PH Wald p | $\hat{\beta} _ {\mathrm{AFT}}$ (SE) | AFT factor $\exp(\hat{\beta} _ {\mathrm{AFT}})$ (Wald 95% CI) | $-\hat{\beta} _ {\mathrm{PH}}/\hat{\gamma}$ | Poisson GLM $\hat{\beta}$ ($\widehat{\mathrm{HR}}$) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `age` (per year) | 0.0102 (0.0092) | 1.0103 (0.9922, 1.0287) | 0.268 | -0.0075 (0.0068) | 0.9926 (0.9795, 1.0058) | -0.0075 | 0.0102 (1.0103) |
| `sex` (female vs. male) | -0.5486 (0.1673) | 0.5778 (0.4162, 0.8020) | 0.00104 | 0.4011 (0.1237) | 1.4935 (1.1718, 1.9033) | 0.4011 | -0.5091 (0.6011) |
| `ph.ecog` (per 1 unit) | 0.4646 (0.1137) | 1.5913 (1.2735, 1.9884) | 4.38e-05 | -0.3396 (0.0835) | 0.7120 (0.6046, 0.8386) | -0.3396 | 0.4050 (1.4993) |

With estimated Weibull shape $\hat{\gamma}$ = 1.3678 ($\hat{\sigma}$ = 0.7311, log-likelihood = -1132.4), female patients have a 42.2% lower death hazard (HR = 0.5778, p = 0.00104). In the AFT parameterisation, this corresponds to a 49.4% longer survival time ($\exp(\hat{\beta} _ {\mathrm{AFT}})$ = 1.4935). Each 1-point increase in `ph.ecog` raises the death hazard by 59.1% (HR = 1.5913, p = 4.38e-05) and shortens survival time by 28.8% ($\exp(\hat{\beta} _ {\mathrm{AFT}})$ = 0.7120).

![Empirical KM vs Fitted Weibull Survival and PH-AFT Duality](figures/weibullFit.png)

---

## Reversible illness-death model (10,080 person-years)

This demo (`R/multistate.R`) models a continuous-time process $Z(t) \in \{1, 2, 3\}$ with state 1 (healthy), state 2 (ill), and absorbing state 3 (dead) over 3 years. Both setups in the script come from textbook exercises. First, I plot closed-form direct mortality ($1\to 3$, with $\lambda _ {12}$ = 0.2 per year and $\lambda _ {13}$ = 1.0 per year) against illness-mediated mortality ($1\to 2\to 3$, with $\lambda _ {23} = \alpha \in \{1, 2, 5\}$ per year). Second, I analyse a synthetic 3-year cohort with 2,480 person-years at risk in state 1 (980, 800, and 700 across years 1 to 3) and 7,600 person-years in state 2 (2,500, 2,600, and 2,500 across years 1 to 3).

### Testing time-homogeneity and estimating 3-year sojourn times

Let $N _ {kl,r}$ be the transition count from state $k$ to state $l$ and $R _ {k,r}$ the person-years at risk in state $k$ during year $r \in \{1, 2, 3\}$. The year-specific and pooled maximum-likelihood transition rates (per year) are $\hat{\lambda} _ {kl,r} = N _ {kl,r}/R _ {k,r}$ and $\hat{\lambda} _ {kl} = N _ {kl,\cdot}/R _ {k,\cdot}$. Within year $r$, the process has generator matrix

$$
Q_r = \begin{pmatrix}
-(\lambda_{12,r} + \lambda_{13,r}) & \lambda_{12,r} & \lambda_{13,r} \\
\lambda_{21,r} & -(\lambda_{21,r} + \lambda_{23,r}) & \lambda_{23,r} \\
0 & 0 & 0
\end{pmatrix}.
$$

I test whether the transition rates are constant across the three years (8 df) with the Poisson likelihood-ratio statistic

$$
W = 2\sum_{(k,l)}\sum_{r=1}^{3} N_{kl,r}\log\left(\frac{\hat{\lambda}_{kl,r}}{\hat{\lambda}_{kl}}\right) \sim \chi^2_{8}.
$$

Because time-homogeneity is rejected (W = 460.503 on 8 df, p = 2.08e-94), the constant-rate infinite-horizon sojourn formula $1/\hat{\lambda} _ {k\cdot}$ (2.3006 years in state 1 and 4.0729 years in state 2) does not fit the data. Letting $\lambda _ {k\cdot,r} = \sum _ {l\ne k}\lambda _ {kl,r}$ denote the total exit rate from state $k$ in year $r$, I integrate the piecewise-exponential survival curve $S_k(t)$ of the holding time $T_k$ over $[0, 3]$ to get the 3-year restricted mean sojourn time (in years):

$$
\mu_k(3) = \mathbb{E}\left[\min(T_k, 3)\right] = \int_0^{3} S_k(t)\mathrm{d}t
= \sum_{r=1}^{3}\exp\left(-\sum_{j < r}\lambda_{k\cdot,j}\right)\frac{1 - \exp(-\lambda_{k\cdot,r})}{\lambda_{k\cdot,r}}.
$$

I compute its variance by the delta method across the three independent annual exit rates, using $\mathrm{Var}(\hat{\lambda} _ {k\cdot,r}) = N _ {k\cdot,r}/R _ {k,r}^{2}$:

$$
\mathrm{Var}\left(\hat{\mu}_k(3)\right) = \sum_{r=1}^{3}\left(\frac{\partial\mu_k(3)}{\partial\lambda_{k\cdot,r}}\right)^{2}\frac{N_{k\cdot,r}}{R_{k,r}^{2}}.
$$

To allow re-entry after recovery ($2\to 1$) or relapse ($1\to 2$), I also compute the total 3-year occupation time $\Psi _ {kk}(3)$ (in years) from the piecewise matrix exponential $P(r - 1, r - 1 + u) = \exp(u Q_r)$:

$$
\Psi_{kk}(3) = \int_0^{3} P_{kk}(0, t)\mathrm{d}t.
$$

| Transition / quantity | Year 1 | Year 2 | Year 3 | Pooled / 3-year estimate | SE | 95% CI | Interval method | LRT $W$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| $1\to 2$ (illness onset, per year) | 0.4031 | 0.5000 | 0.4000 | 0.4335 (pooled) | 0.0132 | (0.4079, 0.4602) | Exact Poisson (Garwood) | 11.779 |
| $1\to 3$ (death from healthy, per year) | 0.0010 | 0.0013 | 0.0014 | 0.0012 (pooled) | 0.0007 | (0.0002, 0.0035) | Exact Poisson (Garwood) | 0.058 |
| $2\to 1$ (recovery, per year) | 0.4160 | 0.1923 | 0.1280 | 0.2447 (pooled) | 0.0057 | (0.2337, 0.2561) | Exact Poisson (Garwood) | 447.539 |
| $2\to 3$ (death from ill, per year) | 0.0008 | 0.0004 | 0.0012 | 0.0008 (pooled) | 0.0003 | (0.0003, 0.0017) | Exact Poisson (Garwood) | 1.127 |
| **Sojourn $\hat{\mu} _ 1(3)$, state 1 (years)** | - | - | - | **1.6808** | **0.0293** | **(1.6234, 1.7381)** | Delta-method Wald | - |
| **Sojourn $\hat{\mu} _ 2(3)$, state 2 (years)** | - | - | - | **1.9273** | **0.0205** | **(1.8871, 1.9675)** | Delta-method Wald | - |
| Occupation $\hat{\Psi} _ {11}(3)$, with re-entry (years) | - | - | - | 1.8853 | - | - | - | - |
| Occupation $\hat{\Psi} _ {22}(3)$, with re-entry (years) | - | - | - | 2.2818 | - | - | - | - |

Because the two mortality rates are based on only 3 deaths ($1\to 3$) and 6 deaths ($2\to 3$), Wald intervals are unreliable and give a negative lower bound (-0.0002 per year) for $1\to 3$. I therefore report exact Poisson (Garwood) 95% intervals for the pooled transition rates, and delta-method Wald 95% intervals for the sojourn times.

Almost all of the non-homogeneity (W = 447.539 out of 460.503, or 97.2%) comes from the recovery rate $2\to 1$, which drops from 0.4160 per year in year 1 to 0.1923 in year 2 and 0.1280 in year 3. Over the 3-year window, the restricted single-spell sojourn time is 1.6808 years in state 1 and 1.9273 years in state 2. Including re-entries, total 3-year occupation time is 1.8853 years in state 1 and 2.2818 years in state 2.

![Illness-Death Mortality Decomposition and Piecewise-Constant Transition Rates](figures/stateOcc.png)

---

## Files and running the code

* `run.R`: runs all five scripts in `R/` and saves the CSV tables and PNG plots.
* `R/prep.R`: loads `survival::veteran` and `survival::lung` and sets up the factor variables (`karno1`, `priorFac`, `trtFac`, `sexFac`).
* `R/nonparam.R`: Kaplan-Meier, Nelson-Aalen, 1-year RMST, log-rank, and Peto-Peto tests (`figures/kmCurve.png`, `results/rmstSum.csv`).
* `R/coxDiag.R`: stratified Cox models, interaction contrasts, and Schoenfeld, Martingale, and `dfbeta` diagnostics (`figures/phDiag.png`, `figures/coxResid.png`, `results/coxCont.csv`).
* `R/weibullAft.R`: Weibull PH vs. AFT models and Poisson regression with a log-time offset (`figures/weibullFit.png`, `results/weibullDual.csv`).
* `R/multistate.R`: illness-death mortality trajectories, homogeneity likelihood-ratio test, exact Poisson (Garwood) intervals, and 3-year sojourn and occupation times (`figures/stateOcc.png`, `results/stateSoj.csv`).

To run everything from the terminal (requires R and the `survival` package):

```bash
Rscript run.R
```
