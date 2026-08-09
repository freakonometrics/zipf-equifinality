# Simulation protocol — v0.6

## Fixed primary design

- Corpus size: `n = 100000` tokens.
- Common finite support for i.i.d. Zipf, Markov and SSR: `V = 5000`.
- Primary rank window: ranks 10--500, with empirical count at least 5.
- Held-out paper replications: 100, beginning at seed 2001.
- Within-sequence shuffle surrogates for excess NMI: 20 per replication.
- Primary collapsed alphabet: top 50 frequent types + `OTHER`.
- NMI sensitivity: top 25, 50 and 100 types + `OTHER`.
- Bootstrap rank-effect replicates: 2000.
- Markov persistence in the primary benchmark: `rho = 0.35`.
- Common i.i.d./Markov block length: 250.

## Exact-marginal primary controls

1. **Finite i.i.d. Zipf:** `p_j = 1/(j H_V)`.
2. **Sticky Markov Zipf:** `P_ij = rho 1{i=j} + (1-rho) p_j`; the same `p` is
   exactly stationary and the chain is reversible.
3. **Canonical SSR:** from `i > 1`, move uniformly to `{1,...,i-1}`; after state
   1 restart uniformly on `{1,...,V}`. Its exact stationary distribution is
   `pi_j = 1/(j H_V)`. The observed chain is initialized in stationarity.
4. **Latent mixture:** conditionally non-Zipfian sequence distributions pooled
   over a calibrated distribution of latent scales.

SSR is not calibrated to the fitted exponent in v0.6.

## Calibration / evaluation separation

Random segmentation and Simon reinforcement use calibration seeds beginning at
1001. The latent mixture uses a 50-seed final calibration panel, 1001:1050 by
default, after a smaller prespecified screening panel.

No held-out seed beginning at 2001 is used to select calibration parameters.
The exact i.i.d., Markov and SSR controls require no exponent calibration.

## Primary diagnostics

- OLS Zipf exponent;
- restricted-domain MLE on the same empirical rank window;
- final realized vocabulary size;
- global and early/late Heaps slopes;
- excess lag-1 normalized mutual information;
- transition direction `P(next rank is more frequent | rank changes)`.

## Review-driven robustness checks

- `r_max` sensitivity: 300, 500, 1000, 2000, 3000;
- `r_min` sensitivity: 1, 5, 10, 20, 50, 100;
- Markov persistence: 0, .10, .20, .35, .50, .65, .80;
- effective-information comparison to i.i.d. samples of size `N_eff(rho)`;
- exact expected realized vocabulary under the implemented Markov block design;
- SSR sequence diagnostics under natural cycles and common fixed blocks;
- NMI alphabet collapse at K = 25, 50, 100;
- latent-scale conditioning;
- random segmentation and Simon reinforcement under changing observation windows.

## Interpretation

The benchmark is a model-discrimination exercise under specified observation
designs. `I_excess` detects order dependence relative to the shuffle null; it is
not SSR-specific. High within-cycle `D` is a prespecified positive control for
canonical SSR contraction, not a unique mechanistic signature in arbitrary
data.
