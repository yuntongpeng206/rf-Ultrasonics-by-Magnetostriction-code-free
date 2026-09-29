# Private input specification

## Required files

Place these in `private_data/` in the project root:

| File | Columns |
|---|---|
| cw_scans.csv | mode, frequency_kHz, amplitude_mV |
| pulses.csv | pulse_number, delay_us, amplitude_mV |

Each CW mode contains one selected scan, at least five distinct frequencies, and a bracketed peak. Crossing analysis must cover both half-power crossings. Provide enough modes that at least three distinct modes remain after exclusion. Do not pool different runs into one mode.

Pulses must have distinct pulse numbers, increasing delays and positive amplitudes. At least three are needed. Use one consistent amplitude convention per dataset.

Optionally supply `observed_modes_optional.csv` with columns `mode,frequency_kHz`. This replaces fitted centers for the CW speed regression; frequency uncertainties are interpolated from detailed linewidths using the original half-width rule. Outside the measured mode range, NumPy interpolation holds endpoint widths constant. Do not submit an empty optional CSV.

## Settings

The template deliberately supplies no measured values and no automatic scientific choices. Required null entries must be filled.

| Setting | Meaning |
|---|---|
| cw_method | `lorentzian` or `crossings` |
| L_cw_m, L_pulse_m | Positive rod lengths in metres for the two measurements |
| u_L_cw_m, u_L_pulse_m | Nonnegative length uncertainties in metres |
| u_delay_us | Positive common pulse arrival-time uncertainty in microseconds |
| fit_upto | Number of initial pulses used for decay fitting; null uses all |
| excluded_modes | Mode numbers excluded from CW speed and selected attenuation summary; [] excludes none |
| frequency_width_fraction | Positive fraction of full linewidth used for CW center uncertainty when no separate observed-mode list is supplied |
| frequency_error_floor_kHz | Positive lower bound for that assigned center uncertainty |
| linewidth_fractional_error | Nonnegative fractional linewidth-error assumption |
| pulse_decay_fractional_error | Nonnegative assumed fractional uncertainty in fitted pulse decay rate |

The source notebook used hand-chosen uncertainty floors and excluded modes in some versions. The template makes these explicit local choices rather than interpreting them as instrument specifications. `fit_upto` selects by sorted pulse order, not an amplitude threshold.

## Local results

`private_outputs/` contains a CW table, a JSON summary and a figure. These reveal measurements/results and must remain local. Do not upload an executed notebook with outputs. Git ignore rules do not remove previously committed files and do not prevent manual website uploads.
