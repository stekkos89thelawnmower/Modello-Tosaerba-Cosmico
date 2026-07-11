# Assembly bias, splashback radius, and AGN feedback: a test with real data

*Author: Stefano Boi*

*Informal project name: "Cosmic Lawnmower"*

This repository documents a two-phase journey:

1. **`halo_simulation.py`** — a toy model that formalizes a hypothesis:
   two-mode AGN feedback (radiative/kinetic) leaves an independent
   imprint on the splashback radius of dark matter halos, beyond what
   is already explained by the known correlation between mass
   accretion rate (Gamma) and splashback radius.
2. **`tng_real_data_test.py`** — the test of that hypothesis against
   real data from the IllustrisTNG simulation (TNG100-1), including the
   debugging process that required several fixes before reaching a
   reliable result.

## Final result

**The hypothesis is not testable with the method used in this
repository — not just "unsupported", but in principle incapable of
being confirmed or refuted with these tools.**

The path here was longer than expected, and it's worth reporting in
full because the final conclusion only emerged after several rounds of
testing (see "Debugging journey" below). In summary:

1. A first test (central halos, kinetic/radiative feedback vs
   splashback residual) initially appeared to show a signal (p<0.001),
   which disappeared after controlling for M_BH (p=0.139) — a mass
   confound.
2. A second test (satellites, baryon stripping vs the central's
   feedback mode) showed a signal (p=0.037) that disappeared once
   correctly aggregated by host instead of by satellite
   (pseudoreplication: 600 satellites are not 600 independent
   observations, only ~30 are).
3. A third test (cumulative black hole activity, quiescent vs grown,
   controlling for Gamma and halo mass) showed a very strong signal
   (p<0.0001) that turned out to be **structurally impossible to
   interpret as a feedback effect**: `R_sp/R200m` in this entire
   pipeline was never actually measured from the simulation — it was
   always computed analytically by Colossus as a deterministic function
   of only two quantities, Gamma and halo mass (`nu200m`). No channel
   exists, anywhere in this calculation chain, through which AGN
   feedback (or any other black hole physics) could influence the R_sp
   value used here. The signal found was therefore necessarily a
   residual of a non-linear correlation between M_BH and (Gamma, mass)
   not fully removed by the linear control — a statistical-method
   artifact, not feedback physics.

**Honest conclusion**: none of the tests in this repository can, by
construction, confirm or refute the original AGN feedback/splashback
hypothesis, because the tool used to compute R_sp (Colossus, from Gamma
and mass) does not contain that information. A genuine test requires a
DIRECT measurement of R_sp from the simulation (see "Why we didn't
simply measure R_sp directly" below).

## Why we didn't simply measure R_sp directly

Directly measuring the splashback radius requires tracking individual
particle orbits around each halo (the SPARTA algorithm, Diemer 2017)
or a full density profile — not just catalog-level properties (mass,
virial radius) obtainable from the search API used in this repository.
A direct precedent exists — one study applied SPARTA to IllustrisTNG to
measure splashback for 812 cluster-mass halos (~1e14 Msun), studying
the relationship between splashback and the accretion shock — but it
is not a public catalog downloadable via API, it is limited to cluster
scale (not the 1e11-1e12 Msun band used here), and Diemer's public
catalogs with measured splashback (benediktdiemer.com/data) cover the
Erebos N-body simulations, not TNG.

Doing this properly would require downloading full snapshots with
particle positions over time (tens/hundreds of GB, not the few hundred
lightweight HTTP requests used so far) and running SPARTA independently
— a jump in scale that makes the Colab/mobile workflow used in this
repository no longer adequate.

## Debugging journey: an honest account

The code didn't work on the first try. I list the problems encountered,
in chronological order, because each one is instructive:

1. **Wrong API endpoint.** The assumed "bulk field download" endpoint
   (`.../halos/?fields=...`) does not exist in the current TNG API; it
   returned the homepage with a 200 status, masking the error. Fixed by
   switching to a search+detail pattern verified empirically against
   the live API (`/subhalos/?primary_flag=1&mass__gt=...` followed by
   per-subhalo detail and then the parent group's `info.json`).

2. **Fragile authentication on Colab.** The API key was lost between
   session restarts or remained an unreplaced placeholder, producing
   misleading 403 errors ("Failed session auth").

3. **Missing parameter in Colossus.** `splashback.splashbackModel()`
   requires an explicit `rspdef` parameter (e.g. `sp-apr-mn`), not
   documented as optional in the first call.

4. **A ~1e15 factor error in the sigma_T constant.** The conversion of
   the Thomson cross-section from cm^2 to kpc^2 contained a typo in the
   exponent, making the lambda_Edd calculation systematically wrong for
   every halo. Fixed by deriving the Eddington constant directly in CGS
   and verifying it numerically against the known Salpeter time (~45
   Myr for eta=0.1).

5. **Wrong feedback transition threshold.** A fixed threshold
   (lambda_Edd >= 0.01) was replaced with the one actually implemented
   in TNG (Weinberger et al. 2017): a threshold that depends on black
   hole mass, chi(M_BH) = min(0.002*(M_BH/1e8 Msun)^2, 0.1).

6. **Sampling bug.** The `mass__gt=threshold` filter with `limit=N`
   always returned the most massive halos above the threshold (because
   TNG IDs are ordered by decreasing mass across the whole box), not a
   representative sample. Fixed by adding an upper mass bound
   (`mass__lt`) as well, to force the inclusion of smaller halos closer
   to the true radiative/kinetic transition (M_BH ~ 10^8.2 Msun).

7. **The first confound check.** The residual test on the central halo
   showed a signal (p<0.001) that a check controlling for M_BH
   dismantled (p=0.139) — a mass confound.

8. **Pseudoreplication in the satellite test.** An extension of the
   test to satellites (baryon stripping as a function of distance from
   the center) showed a signal (p=0.037) across 600 satellites, which
   disappeared once correctly aggregated by host (effective N=30, not
   600).

9. **The final structural discovery.** A test on cumulative black hole
   activity (quiescent vs. grown) showed a very strong signal
   (p<0.0001) that turned out to be physically uninterpretable:
   R_sp/R200m in this pipeline was always computed analytically by
   Colossus as a function of only Gamma and mass, so it cannot, by
   construction, contain any information about AGN feedback. The signal
   was a residual of a non-linear correlation between M_BH and (Gamma,
   mass), not feedback physics. See "Why we didn't simply measure R_sp
   directly" for details.

## Stated methodological limitations

- **Approximate definition of Gamma**: computed between z=0 and z~0.5
  from the main progenitor branch of the SubLink merger tree; not
  identical to the definitions used in all reference papers.
- **R_sp not measured directly — this is the central limitation, not a
  minor one** (see "Final result" and "Why we didn't simply measure
  R_sp directly" above): it is derived from Gamma and mass via
  Colossus, so it cannot contain information about AGN feedback by
  construction.
- **Only z=0, only TNG100-1**: no robustness check against other
  redshifts or other simulations (EAGLE, SIMBA) with different AGN
  feedback implementations — moot anyway until the limitation above is
  resolved.

## Next steps, if pursuing this further

The only next step that could genuinely answer the original question is
replacing the computed R_sp (Colossus, from Gamma and mass) with R_sp
directly measured from the simulation (SPARTA, from particle orbits).
Every other extension discussed during development (different
redshift, larger sample, different simulation) would not fix the
structural problem: they would keep measuring the residual of a
deterministic function of Gamma and mass, not a physical feedback
effect.

If this is pursued in the future, the computational cost is
substantially higher (downloading full snapshots with particle data,
not just group catalogs) and likely requires an environment other than
Colab/mobile — a compute cluster or at least a machine with adequate
storage and RAM.

### Why even a "lightweight" attempt at SPARTA via the API doesn't work

After checking the official API documentation, the problem isn't just
one of scale — it's structural. TNG's "cutout" mechanism (which we
would use to fetch a halo's particles) only returns **particles
currently bound to a subhalo/group at that specific snapshot**. But the
particles relevant to splashback are exactly those that, after being
inside the halo, have since bounced back out — so in more recent
snapshots they are no longer "members" of any group, and a
membership-based cutout excludes them right at the moment they become
interesting. There is no documented endpoint for requesting "all
particles within a sphere around a point," independent of group
membership — only membership-based cutouts, or downloading the entire
snapshot (with offset tables to organize hundreds of GB-TB of data per
box).

The only realistic alternative found is the **JupyterLab environment
hosted directly by the TNG project** (tng-project.org/data/lab/), where
the raw files are already available locally on their cluster, with no
need to download them. This would be the concrete path for a serious
SPARTA attempt, but it's a different environment, likely requires a
separate access request, and isn't something that can be orchestrated
from within this working session.

### Who can realistically access this final step

Access to TNG's JupyterLab is not automatic with basic registration
(the kind used for this repository's API key). According to the
official documentation, the service is offered **on an experimental
(beta) basis, and only to active academic users**, following a
separate, specific access request. In practice, this means
researchers, PhD students, students supervised by a research group, or
participants in schools/events officially recognized by the project (a
concrete example found: a winter school where participants had to
register with their institutional email and specify the event name in
the request).

This makes the final step (directly measuring R_sp via SPARTA) a
limitation not just of computation but also of access: it requires an
institutional affiliation, not just time and compute resources. A user
with a plain public API account (like this repository) does not
currently have a direct path to complete it.

## Files in this repository

- `halo_simulation.py` — original synthetic toy model
- `tng_real_data_test.py` — first real-data test on TNG100-1 (central halo, feedback vs splashback)
- `tng_unified_test.py` — extension to satellites (baryon stripping)
- `tng_bh_activity_test.py` — test on cumulative black hole activity (the one that revealed the structural limitation)
- `LICENSE` — MIT

## Acknowledgments

This project was developed with the assistance of Claude (Anthropic)
for code writing, debugging, literature checks, and drafting this
documentation. All key conceptual decisions — insisting on every
robustness check, requesting a unified comparison, recognizing when a
technical limitation was also a scientific one — were driven by the
Stefano Boi, the repository's author.

## License

This project is distributed under the MIT license — see `LICENSE` for
the full text. In short: anyone can use, modify, and redistribute this
code, including for commercial purposes, provided the copyright notice
is retained. No warranty is provided regarding the scientific
correctness of the results (see "Stated methodological limitations"
above).
