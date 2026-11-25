# resolution-gap-every-state

**Every state keeps its arrow: an unconditional lower bound on the resolution gap of the arrow of time**

Code, stored results, interval-arithmetic checks and a verification script for the paper by Taishi Namba (2026).

**Paper:** T. Namba, *Every state keeps its arrow: an unconditional lower bound on the resolution gap of the arrow of time*, preprint, Zenodo (2026), doi: to be added.

[![DOI](https://img.shields.io/badge/DOI-to%20be%20added-lightgrey.svg)](#citation)
[![verify](https://github.com/neginyan/resolution-gap-every-state/actions/workflows/verify.yml/badge.svg)](https://github.com/neginyan/resolution-gap-every-state/actions/workflows/verify.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## The result in brief

An observer sees a state of the pendulum near its upper equilibrium through Gaussian noise of width σ; the resolution gap g measures how much of the true expansion rate the blurred description loses. [Paper E](https://doi.org/10.5281/zenodo.23617343) conjectured that g never falls below c′/√(ℓ + C′), ℓ = ln(1/σ), whatever the state. Papers F, H, I and J proved the lower side for states built from a class of windows. **This paper proves a lower bound for every state, unconditionally: no condition on the shape of the state.**

- **Theorem 1 (every state keeps a gap).** The only assumptions are a finite first moment, E|X| < ∞, and the confinement of the problem itself, Var(X_z) ≤ π² (the observer looks at the neighbourhood of the equilibrium). For every such state and every σ ≤ 1.5×10⁻¹⁴,
  g ≥ ¾ α/(ασ + 1268 k̄),  k̄² = 289 L̄^{5/2}/(ln 2L̄)^{3/2} + 4L̄ + 3.57,  α = 1/(2(π² + σ²)),  L̄ = 2ℓ + ln(1/α);
  for σ ≤ 10⁻¹⁵ this gives g ≥ 1.7×10⁻⁶ (ln 2L̄)^{3/4}/L̄^{5/4}. For σ ≤ 10⁻⁴, every state outside a small open set O also has the larger bound of part (a).
- **The present power, range and constants.** The power ℓ^{−5/4}(ln ℓ)^{3/4} is the present one; Paper E's conjecture has ℓ^{−1/2}. The range σ ≤ 1.5×10⁻¹⁴ is narrow and the constants are small. The paper states all three plainly, as what the next steps must improve; within its range the theorem holds for every state.
- **Three new steps.** (1) A three-circle argument in complex time bounds the effective wavenumber of every state: k_eff² ≤ 289 L^{5/2}/(ln 2L)^{3/2} + 4L + 3.57, L = ln(1/A) (Theorem 8). (2) Once the information exceeds 112σ², the remainder of Paper E's gap formula is at most half of its main terms (Lemma 11), so every remainder estimate is needed only on an open set O where the information is of order σ². (3) On O the gap is expanded around the observed point instead of the posterior mean; the new remainder carries √(A + B), and the one second-order quantity left is controlled by the same circles in the isotropic flow (Lemma 16, Theorem 19).
- **Every constant by interval arithmetic.** No step of the proof is computer-assisted, but all 101 constants were recomputed with outward-rounded 60-digit decimal interval arithmetic, and the 34 joints of the chain (the conclusion of one lemma inserted into the hypothesis of the next) were evaluated on 6062 recorded states.
- **The printed numbers are checked.** Every number printed in the text that comes from a computation (208 entries, with an inventory of the 196 numbers of the running text) is recomputed from the stored records or the interval arithmetic and checked for the direction of its rounding; ten deliberately wrong versions are rejected.

Every statement in the paper carries its status: proof, computer-assisted proof, derivation, numerical, or conjecture.

## Verify

```bash
pip install -r requirements.txt
python verification/quick.py                     # under a minute; no other repository needed
python verification/run_all.py --e-repo ../resolution-gap-log-law --f-repo ../resolution-gap-two-sided   # about half an hour
```

`quick.py` reruns every constant of the chain by interval arithmetic (App. C), the 34 joints from the stored records (App. D), the range of σ as a function of the constants (Table 2), the independent checks on three point masses and two Gaussian blobs (Lemmas 4, 6, 16, 18), our recomputation of the three-point check made in the review of the manuscript, the generated tables of the paper, and the check of the printed numbers. It is the check of the GitHub Actions workflow.

`run_all.py` runs `quick.py` and then the three circles step by step on the recorded states (two-dimensional, isotropic, one-variable), the real-time conclusions on the 2449 two-dimensional states, engine states of Paper F's code re-evaluated afresh without cache, and checks that the formulas quoted from Paper E are in its text. The engine is imported from the published repository of Paper F (v1.5.0), which is not copied here: pass its location with `--f-repo` (or `F_REPO`); Paper E's repository with `--e-repo` (or `E_REPO`). With `--full` it also re-evaluates all 692 engine states (about an hour).

Every rerun is made in a temporary copy and compared with the stored record entry by entry with the comparison of the published `common.py` of Papers I and J (flags, counts, strings and interval bounds exactly; double-precision floats up to a relative 10⁻⁹); the stored records are never overwritten.

## What is where

| result | status | script → results |
|---|---|---|
| Theorem 1 (assembly; parts (a), (b), (c)) | proof | `paper/sections/07_assembly.tex`; constants: `src/ei_constants.py` → `results_EI/ei_constants.json` |
| Theorem 8 (k_eff² ≤ ψ for every state; three circles) | proof; every step numerical on recorded states | `src/eg_circ.py 2d`, `1d`, `1ds` → `results_EG/eg_circ_*.json`; `src/ef_q2.py 2d` → `results_EF/ef_q2_2d.json` |
| Lemmas 10, 11, 14, Proposition 15 (outside O) | proof | joints: `src/eh_joints.py` → `results_EH/eh_joints.json` |
| Lemma 16, Lemmas 17, 18, Theorem 19 (inside O) | proof | `src/eg_circ.py iso`, `src/eg_point.py`, `src/ei_point.py`, `src/eg_direct.py` → `results_EG/`, `results_EI/ei_point.json` |
| every constant of the chain (101 claims) | interval arithmetic | `src/eh_iv.py`, `src/eh_constants.py`, `src/ei_constants.py` → `results_EH/eh_constants.json`, `results_EI/ei_constants.json` |
| the 34 joints on 6062 recorded states | numerical | `src/eh_joints.py` → `results_EH/eh_joints.json` |
| the 692 engine states (Paper F's code) | numerical | `src/eh_states.py`, `src/eh_fresh.py`, `src/eg_recheck.py` (with `src/ed_engine.py`, `src/ee_engine.py`, `src/ef_wide.py`) → `results_EH/`, `results_EG/` |
| the range of σ (Table 2) | derivation | `src/eh_range.py` → `results_EH/eh_range.json` |
| the tables of App. C and D, the chain figure | generated from the records | `src/ei_tables.py`, `src/ei_figure.py` |
| the printed numbers of the text | checked against the records | `verification/text_numbers.py` → `results_EI/ei_text_numbers.json` |

| folder | content |
|---|---|
| `src/` | the code (Python 3, NumPy, SciPy, Matplotlib) of the lab stages ED–EI behind this paper |
| `results_ED/` … `results_EG/` | the stored records of the recorded states and of the step-by-step checks (stages ED–EG) |
| `results_EH/`, `results_EI/` | the constants, the joints, the engine states, the range, the point checks, the printed numbers |
| `verification/` | `quick.py`, `run_all.py`, `text_numbers.py`, `common.py` |
| `docs/` | `chain_of_theorem_EG.md`: the whole proof as one chain, item by item, with the joints and the corrections of the external review |
| `paper/` | LaTeX source of the paper (REVTeX 4.2), its figure, and a check PDF built with a stand-in class (`paper/paperK_check_standin.pdf`; the final PDF is built with REVTeX) |

The working record of the proofs (the stage reports of the lab sessions) is kept by the author.

## Known issues

The two notes of 2026-11-23 and 2026-11-24 in the README of Paper F's repository concern its Conjecture 18 and the step from it to Paper E's conjecture in case 3. This paper uses neither: it replaces them by Theorem 8 (k_eff² ≤ ψ for every state), the window inequality of Paper F's Theorem AT2 with its factor σ², and the Cramér–Rao floor, so the notes do not affect Theorem 1.

## Related work by the author

- T. Namba, *Is the resolution gap of the arrow of time uniform?*, Zenodo (2026), [doi:10.5281/zenodo.23617343](https://doi.org/10.5281/zenodo.23617343) (paper E; code [resolution-gap-log-law](https://github.com/neginyan/resolution-gap-log-law)).
- T. Namba, *How fine can a hidden ripple be? A two-sided logarithmic law for the resolution gap of the arrow of time*, Zenodo (2026), [doi:10.5281/zenodo.23746971](https://doi.org/10.5281/zenodo.23746971) (paper F, all versions; code [resolution-gap-two-sided](https://github.com/neginyan/resolution-gap-two-sided), [doi:10.5281/zenodo.23746555](https://doi.org/10.5281/zenodo.23746555)).
- T. Namba, *Beyond few harmonics: the lower side of the logarithmic law of the resolution gap for windows with a weak ripple*, v1.3, Zenodo (2026), [doi:10.5281/zenodo.24369012](https://doi.org/10.5281/zenodo.24369012) (paper H; code [resolution-gap-weak-ripple](https://github.com/neginyan/resolution-gap-weak-ripple)).
- T. Namba, *Wherever the observer focuses: the lower side of the logarithmic law of the resolution gap at every focus and every resolution*, Zenodo (2026), [doi:10.5281/zenodo.24371173](https://doi.org/10.5281/zenodo.24371173) (paper I; code [resolution-gap-every-focus](https://github.com/neginyan/resolution-gap-every-focus), [doi:10.5281/zenodo.24370921](https://doi.org/10.5281/zenodo.24370921)).
- T. Namba, *Wherever the observer focuses, whatever the ripple: the lower side of the logarithmic law of the resolution gap for every window of the class*, Zenodo (2026), [doi:10.5281/zenodo.24418628](https://doi.org/10.5281/zenodo.24418628) (paper J, all versions; code [resolution-gap-every-ripple](https://github.com/neginyan/resolution-gap-every-ripple), [doi:10.5281/zenodo.24418549](https://doi.org/10.5281/zenodo.24418549)).

Papers E, F, H, I and J and their repositories are not changed by this work.

## Use of AI assistants

The study was carried out with AI assistants (Anthropic Claude, Google Gemini and OpenAI ChatGPT). Their roles are described in the acknowledgments of the paper. The author directed the research and takes full responsibility for it.

## Citation

Please cite the paper (doi: to be added) and this repository; see `CITATION.cff`. The code will be archived at Zenodo (doi: to be added).

## License

MIT © 2026 Taishi Namba
