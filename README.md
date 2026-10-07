# Meson Transition Atlas

**Site:** https://mmikhasenko.github.io/MesonTransitionAtlas/

Level diagrams for ten meson flavor sectors of the Godfrey–Isgur relativized
quark model, with an arrow for every decay that
[GIModel.jl](https://github.com/mmikhasenko/GIModel.jl) can compute, labelled
by its partial width in MeV.

- **Model:** S. Godfrey and N. Isgur, *Mesons in a relativized quark model with
  chromodynamics*, [Phys. Rev. D 32, 189 (1985)](https://doi.org/10.1103/PhysRevD.32.189).
- **Reproduction of the model:** M. Mikhasenko, *A Full Reproduction of the
  Godfrey–Isgur Relativized Quark Model*,
  [arXiv:2609.37716](https://arxiv.org/abs/2609.37716)
  ([INSPIRE](https://inspirehep.net/literature/3209266)).
- **Computation:** [GIModel.jl](https://github.com/mmikhasenko/GIModel.jl)
  v0.4.2 at commit
  [`cc7bf22`](https://github.com/mmikhasenko/GIModel.jl/commit/cc7bf22efc7611bb1fddb9fefa759184b1e67d0b),
  pinned in [`Project.toml`](Project.toml) and [`Manifest.toml`](Manifest.toml).
  [Documentation](https://mmikhasenko.github.io/GIModel.jl/dev/).

## Sectors

| tab | quarks | tab | quarks |
|---|---|---|---|
| Light isovector | u d̄ | Charmonium | c c̄ |
| Light isoscalar | n n̄ + s s̄ | Bottom | q b̄ |
| Strange | u s̄ | Bottom-strange | s b̄ |
| Charm | c q̄ | Bottom-charm | c b̄ |
| Charm-strange | c s̄ | Bottomonium | b b̄ |

Levels: S waves with n ≤ 3, P and D waves with n ≤ 2, `OscillatorSolver()`,
the paper's parameters. Isoscalars use the annihilation mixing of the paper
(`PaperP2Annihilation` for pseudoscalars, Table III amplitudes for vectors and
tensors).

## What is computed

| class | operator | where |
|---|---|---|
| strong: π, K, η, η′ emission | `PseudoscalarEmission` | all sectors with a light quark |
| radiative E1, M1, M2 | `PhotonEmission` | every sector, per charge state |
| gg, ggg | `GluonicAnnihilation` | c c̄, b b̄ |
| γγ, e⁺e⁻ | `TwoPhotonAnnihilation`, `LeptonicCurrent` | neutral hidden flavor |
| ℓν | `LeptonicCurrent` | charged ground-state pseudoscalars |

Every pair of states is tried with every operator. A refused call is kept,
grouped by message, and listed on the site under "Not computed in this sector".
The site explains the conventions: strong couplings, charge-channel sums, the
S+P fallback for D-wave admixtures, and exact zeros.

## Repository layout

| path | what it is |
|---|---|
| [`compute.jl`](compute.jl) | the precompute: solves all sectors, evaluates all operators, writes the data |
| [`run.sh`](run.sh) | full run: spectra, parallel strong pass, assembly |
| [`data/transitions.json`](data/transitions.json) | the result: levels, widths, refused calls, metadata |
| [`site/index.html`](site/index.html) | the renderer, a single HTML file without dependencies |
| [`site/data.js`](site/data.js) | the same data as a script, so the page also opens from disk |
| [`.github/workflows/pages.yml`](.github/workflows/pages.yml) | deploys `site/` to GitHub Pages |

## Reproduce

Requires Julia 1.11. The environment is fully pinned, including GIModel.jl.

```sh
git clone https://github.com/mmikhasenko/MesonTransitionAtlas.git
cd MesonTransitionAtlas
sh run.sh 3
```

`run.sh` instantiates the environment, solves the spectra, runs the strong
pass in 3 worker processes, and writes `data/transitions.json` and
`site/data.js`. It takes about 45 minutes on a laptop. Each Julia process needs
1–2 GB of memory; on an 8 GB machine do not use more than 3 workers.

Intermediate results are cached in `cache/` (not tracked):

| cache | holds | delete after |
|---|---|---|
| `cache/states.jls` | all solved states | a change to spectra or parameters |
| `cache/strong/` | one file per decaying parent | a change to strong decays |

With both caches present, `julia --project=. compute.jl` redoes only the
radiative and annihilation widths and the assembly, in about 4 minutes. The
caches do not record the package version, so delete them yourself after
updating GIModel.jl.

To view the result, open `site/index.html` in a browser. No server is needed.

## Numerical notes

- The strong pass samples the oscillator waves onto a 0.005 GeV⁻¹ mesh and
  drops mixed-state components with |c| < 10⁻⁴. Widths agree with the direct
  oscillator overlaps to 10⁻⁶ relative and are 28–67 times faster. Without
  this, the strong pass takes hours
  ([GIModel.jl#21](https://github.com/mmikhasenko/GIModel.jl/issues/21),
  [#20](https://github.com/mmikhasenko/GIModel.jl/issues/20)).
- Radiative and gluonic widths of states with a D-wave admixture are computed
  on the S+P states and marked †
  ([#22](https://github.com/mmikhasenko/GIModel.jl/issues/22)).
- Channels forbidden by parity, G parity or isospin are counted, not drawn
  ([#23](https://github.com/mmikhasenko/GIModel.jl/issues/23)).
- Kinematics use model masses. Widths are not separately convergence-checked;
  the oscillator solver certifies energies only.

## Data format

`data/transitions.json` has `meta` (versions, solver, couplings) and
`sectors`. Each sector has `levels` (key such as `"1^3P_2"`, name, mass in
GeV, composition), `transitions` (class, `from`, `to` or `toSector`/`toKey`,
emitted particle, `width` in MeV, charge realization, basis `full` or `SP`,
charge channels), `gaps` (refused calls by message) and `strong_zeros`.

## Citing

If you use the atlas, please cite the reproduction paper and the original model:

- M. Mikhasenko, *A Full Reproduction of the Godfrey–Isgur Relativized Quark
  Model*, [arXiv:2609.37716](https://arxiv.org/abs/2609.37716)
  ([INSPIRE](https://inspirehep.net/literature/3209266)).
- S. Godfrey and N. Isgur, Phys. Rev. D 32, 189 (1985),
  [doi:10.1103/PhysRevD.32.189](https://doi.org/10.1103/PhysRevD.32.189).
