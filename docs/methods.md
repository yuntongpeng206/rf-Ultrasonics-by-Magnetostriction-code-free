# Methods and source-code review

## Deduplication and choices

The original notebook contains an interpolation workflow, a pulse-only workflow and a later Lorentzian workflow. Shared weighted regression, ordinary regression and offset-exponential functions are stored once. Unique peak interpolation and Lorentzian functions are both retained. All helper function source bodies are extracted from the original code; the later definition is used where a name repeats. Empty cells, measured inputs, plots and saved outputs are not copied.

The driver is rewritten to read private files, validate inputs and export local results. It retains separate CW and pulse rod lengths, the original uncertainty propagation, optional observed-mode interpolation and user-selected exclusions. It reports selected-mode attenuation using the configurable exclusion list instead of hard-coding one specific mode. Original plot layouts and automatic verbal conclusions are not reproduced. Coil-position metadata and an unused amplitude-reading uncertainty are not needed by the calculations and are not distributed.

## CW speed and linewidth

For standing-wave mode n in a rod of length L:

$$f_n=\frac{nc}{2L}.$$

Weighted regression includes a free intercept. Frequency is fitted in kHz; c = 2 L times slope times 1000 converts to m/s. Assigned frequency errors are either a configurable linewidth fraction with a configurable floor, or the original interpolated half-width values when a separate observed-mode list is supplied. These are practical assumptions, not covariance-derived instrument errors.

The fit model is

$$V(f)=B+\frac{A}{\sqrt{1+((f-f_0)/g)^2}}.$$

Here g is a frequency half-width. The full background-corrected half-power bandwidth is W = 2g and Q = f0/W. The alternative crossing method estimates background from scan edges, refines the peak with three points and interpolates crossings around the observed peak. It does not assume a full nonlinear line-shape fit. Edge background and a three-point parabola can be biased by asymmetry and noise.

## CW attenuation convention requiring review

The source calculates `alpha_cw = 4*pi*W_Hz/c`, where W is already the full half-power width. This is preserved and explicitly labeled `alpha_cw_legacy`.

The supplied teaching material uses energy attenuation alpha, amplitude decay exp(-alpha*c*t/2), and an angular-frequency **half width** approximately alpha*c/2. Consequently, a full bandwidth W in Hz implies alpha = 2*pi*W/c under that approximation. Applying 4*pi to the full width gives a factor-of-two discrepancy. The public code does not silently correct the source calculation or claim agreement between the attenuation methods. This issue should be resolved before interpreting legacy attenuation values physically.

## Pulse speed and decay

Fit delay against pulse number with a free intercept. Successive echoes travel an additional 2L, so c = 2L/T with T the fitted time slope in seconds.

The original decay routine fits A(t)=C+A0 exp(-Gamma*t). It searches 500 nonnegative offset candidates below the smallest selected amplitude, performs log-linear fits and selects by amplitude-space squared residuals. It is not the same as jointly optimizing all parameters using nonlinear least squares. The selected initial-pulse count affects the estimate, and late noisy pulses can distort it. The reported R squared is evaluated across all supplied pulses, including any outside the fit subset.

For time in microseconds, the original energy attenuation relation is alpha_pulse = 2*Gamma*10^6/c. Its rate uncertainty is an assumed configurable fraction, not a fitted confidence interval.

## Interpretation limits

CW linewidth errors combine a fit-derived contribution (where available) and an assigned fractional floor. Shared fitted speed makes attenuation estimates correlated. Across-mode averages and residual-based dispersion indicators are descriptive; they are not proof of absence of dispersion. The velocity difference normalized by quadrature errors assumes independent errors and ignores any shared length calibration. Reflection, coupling and support losses may contribute to observed attenuation.

Method reference: University of Manchester second-year laboratory, *Ultrasonics by Magnetostriction*, updated 2023. The PDF, diagrams, teaching instructions and experimental results are not included.
