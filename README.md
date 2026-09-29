# Ultrasonics by Magnetostriction

Python analysis of continuous-wave (CW) resonances and pulse propagation in a magnetostrictive rod. This code-only repository demonstrates resonance extraction, linear regression, sound-speed estimation and pulse-decay analysis.

**Measurements and results are private.** No original measurement arrays, rod dimensions, experimental uncertainty values, result tables, images, reports or teaching PDFs are distributed. The CSV templates have headers only, and configuration values must be supplied locally.

## Run locally

Requires Python 3.10 or later.

```bash
python -m pip install -r requirements.txt
python src/analysis.py
```

Without data, the script prints setup instructions. To calculate results, create `private_data/` and add completed copies of `settings.json`, `cw_scans.csv` and `pulses.csv` from `templates/`. The optional observed-mode CSV should only be copied if populated. See [input format](docs/input_format.md).

```bash
jupyter notebook notebooks/ultrasonics_analysis.ipynb
```

Results are saved locally in `private_outputs/`, which is ignored by Git. The notebook contains no outputs and does not automatically display measurements or results.

## Structure

| File or directory | Purpose |
|---|---|
| `src/models.py` | Deduplicated original numerical helper functions |
| `src/analysis.py` | Input validation, configurable workflow and local output |
| `notebooks/ultrasonics_analysis.ipynb` | Output-free entry point |
| `templates/` | Empty data and configuration templates |
| `docs/input_format.md` | Columns, units and settings |
| `docs/methods.md` | Equations, preserved choices and limitations |

## Methods retained

Two CW methods are available: original three-point parabolic peak refinement with interpolated half-power crossings, and original Lorentzian-like amplitude fitting. Neither is silently substituted for the other; choose a method in local settings. Pulse arrival times use weighted linear regression. Pulse amplitudes use the original offset-grid search with logarithmic regression.

**CW attenuation needs a width-convention review.** The source code multiplies the full half-power bandwidth by 4 pi / c. That expression appears to treat a full width as a half width. The calculation is preserved under explicitly named `alpha_cw_legacy` outputs, not represented as a validated attenuation result. See [methods](docs/methods.md).

## About this package

Based on the author's analysis code, with ChatGPT-assisted deduplication, private-input separation and documentation. Helper functions are preserved; orchestration is refactored. Mode exclusions, pulse fit cutoff, rod lengths and uncertainty assumptions are supplied locally rather than disclosed or selected automatically. This is not a byte-for-byte copy of the source notebook.

The supplied laboratory teaching material was consulted but is not redistributed. No numerical experimental findings are claimed in this public package.
