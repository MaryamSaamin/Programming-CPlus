**Black Scholes Model**

**Please open each step by clicking on**
![](media/image1.png){width="0.12501093613298336in"
height="0.19168307086614172in"} **and see the explanations.**

**First calculations**

I computed how an option's price changes when the stock price changes,
using three different methods, and compared their accuracy.

option's price: a financial contract that gives us the right to buy or
sell a stock later.

The Black--Scholes model gives the fair price of an option. It uses
inputs like stock price (S), strike (K), volatility (σ), interest rate
(r), and time (T).

From this price, I could calculate **Greeks**, mainly **Delta (Δ)** and
**Gamma (Γ)**.

**Delta (Δ):** How much the option price changes if the stock moves by
\$1

**Gamma (Γ):** How much Delta changes if the stock moves by \$1

**1. Analytic (Exact Formula)**

This is the math formula directly from Black--Scholes.

Example:\
If the stock = 100, strike = 100, σ = 0.2, T = 1\
Δ = 0.58685, Γ = 0.01895

Exact and used as the "truth" to compare other methods.

**2. Finite Difference**

I bumped the stock price slightly and saw how the option price changes.

Example:\
If h = 0.01,\
Δ_fd = 0.58686\
Γ_fd = 0.01894

It works well for medium step sizes(h), but unstable if h is too small
(rounding errors).

**3. Complex-Step Differentiation**

This is a smarter formula, instead of a small real step, I took a tiny
**imaginary** step (i·h) in the input. This avoids rounding errors
completely.

Example:\
Δ_cs = 0.58685,\
Γ_45° = 0.01895,\
Γ_cs_real about 0

It is very accurate, even for extremely small h.

**4. Testing in different step sizes**

I tested many step sizes from 1e−16 to 1e−4 and saved results. This
shows how each method's error changes when h gets smaller.

Result pattern:

-   Finite difference → beat at medium h

-   Complex-step → stable for all h

-   45° complex-step → best for Gamma

  **Method**         **Δ**      **Γ**       **Error (vs exact)**
  ------------------ ---------- ----------- ----------------------
  Analytic           0.586851   0.0189506   ---
  FD                 0.586861   0.0189463   small
  Complex-Step       0.586851   0.000000    Δ exact, Γ unstable
  Complex-Step 45°   ---        0.018951    almost exact

It is learned that the **complex-step method** is the most stable and
accurate way to calculate derivatives numerically. Finite differences
are easier but can lose accuracy for tiny step sizes due to rounding
errors.

**Difference of two scenarios**

  **Scenario**                              **Meaning**                                                                    **Parameters**
  ----------------------------------------- ------------------------------------------------------------------------------ --------------------------------------------------------
  **Scenario 1 --- ATM reference**          A **normal**, well-behaved option: 1 year to expiry and moderate volatility.   S = 100, K = 100, σ = 0.20, T = 1 year
  **Scenario 2 --- Near-expiry, low-vol**   An **extreme** case: almost expired and barely moves.                          S = 100, K = 100, σ = 0.01, T about 1 day (1/365 year)

What This Means Physically

  **Feature**               **Scenario 1 (normal)**                     **Scenario 2 (near-expiry + low vol)**
  ------------------------- ------------------------------------------- ----------------------------------------------
  **Option price**          \~ \$7--8 (has time and volatility value)   Very small (\~ a few cents)
  **Price sensitivity**     Smooth and continuous                       Very sharp near strike
  **Gamma**                 Moderate (about 0.019)                      Extremely large or unstable numerically
  **Delta**                 Smooth between 0--1                         Almost jumps from 0 to 1 around S = K
  **Numerical behaviour**   Stable, easy to differentiate               Sensitive, tiny rounding errors matter a lot

**Scenario 1 -- Normal (ATM)**

-   The option price changes **smoothly** with S.

-   Finite-difference and complex-step methods all behave well.

-   Errors are small, predictable, and form the typical "U-shape" when
    plotted vs h.

All methods look stable; complex-step just slightly better.

**Scenario 2 -- Near-expiry + low vol**

-   The option's price surface is **almost a step function** around S =
    K.\
    (If the stock is even slightly above K, the option is worth S − K;
    otherwise almost 0.)

-   That sharp transition makes numerical derivatives explode or
    fluctuate.

-   The price itself is **very small**, so dividing by tiny differences
    amplifies rounding error.

**Result:**

-   **Finite-difference** Δ and Γ swing wildly --- huge errors.

-   **Complex-step Δ** stays accurate because it doesn't subtract nearly
    equal numbers.

-   **Complex-step 45° Γ** remains the most stable, but still shows
    noise for extreme h.

-   **Real-part CS Γ** fails completely --- returns large constant bias.

  **Scenario**     **Shape of the road**          **Why differentiation is hard**
  ---------------- ------------------------------ ---------------------------------
  **Scenario 1**   easy to measure slope          Small step gives good estimate
  **Scenario 2**   tiny change gives huge slope   Small step ⇒ big numerical jump

The two scenarios differ in volatility and time to expiry.\
In the normal one (1 year, 20% vol), the option price is smooth, so all
numerical methods work well.\
In the near-expiry case (1 day, 1% vol), the option price becomes
extremely sharp and tiny, so finite differences blow up while
complex-step methods remain stable.

**Scenario 1 (ATM Reference)**

**Parameters:** S = K = 100 r = q = 0 σ = 0.20 T = 1 year

  **h_rel**   **Δ_fd**   **Δ_cs**   **err Δ_fd**   **err Δ_cs**   **Γ_fd**   **Γ_cs real**   **Γ_45°**   **err Γ_fd**   **err Γ_cs real**   **err Γ_45°**
  ----------- ---------- ---------- -------------- -------------- ---------- --------------- ----------- -------------- ------------------- ---------------
  1e-2        0.5497     0.5398     9.9 × 10⁻³     0              0.01984    0.03969         0.01984     3.5 × 10⁻⁶     0.01985             3.5 × 10⁻⁶
  1e-4        0.53993    0.53983    1.0 × 10⁻⁴     0              0.01985    0.03970         0.01985     0              0.01985             0
  1e-6        0.53983    0.53983    1.0 × 10⁻⁶     0              0.01985    0.03970         0.01985     1.5 × 10⁻⁶     0.01985             0
  1e-8        0.53983    0.53983    0              0              0          0.04263         0.01985     0.01985        0.01985             0

(Analytic Δ = 0.53982784, Γ = 0.01984763)

**Observations**

**Δ (Delta)**

> Finite difference error decreases steadily as h → smaller.
>
> **Complex-step Δ_cs** = exact (0 error at all h).

**Γ (Gamma)**

> Finite difference Γ \_fd accurate for mid-range h (\~1e-4 to 1e-6).
>
> **Real-part CS Γ** → constant bias (\~0.02 error).
>
> **45° CS Γ** = analytic-level match (error about 0).

**Interpretation & Reasoning**

For moderate volatility (20%) and 1 year expiry, the option price is a
**smooth function** of S.

Finite differences give the typical "U-shape" error curve, truncation
dominates for large h, round-off for small h.

Complex-step avoids subtraction, so round-off → 0 ⇒ flat error line at
machine precision.

**Recommendation**

  **Greek**   **Preferred Method**          **Step Size (h_rel)**   **Reason**
  ----------- ----------------------------- ----------------------- -----------------------------
  Δ           Complex-Step (Δ_cs)           1 × 10⁻⁸                Exact, no sensitivity
  Γ           Complex-Step 45° (Γ_cs_45°)   1 × 10⁻⁶                Highest stability, accuracy
  Backup      Finite Diff (Δ_fd / Γ_fd)     1 × 10⁻⁴ -- 1 × 10⁻⁶    For quick checks only

In the normal ATM case, both Delta and Gamma from complex-step methods
match the analytic Greeks to machine precision; finite differences are
accurate for moderate h, while the real-part complex-step Gamma shows a
consistent bias.

<https://colab.research.google.com/drive/1JFA7ZONjqPX90OI7HIvkFF1vT_3jPNYY?usp=sharing>

![](media/image2.png){width="6.268055555555556in"
height="4.241666666666666in"}

**Observations:**

-   **Delta FD:**\
    Error is smallest near (h\_{rel} about 10\^{-6}).\
    For larger or smaller steps, the error grows, typical "U-shape"
    pattern caused by truncation vs rounding tradeoff.

-   **Delta CS:**\
    Flat line at the bottom, almost zero error across all step sizes.
    Most accurate and stable.

-   **Gamma FD:**\
    Accurate for mid-range h (\~1e−6) but quickly loses accuracy for too
    large or too small h. Sensitive to step size.

-   **Gamma CS real:**\
    Constant large error (\~1e−2). This variant is unstable for small
    option values.

-   **Gamma CS 45°:**\
    Stays near the bottom, best Gamma accuracy. Matches analytic Gamma
    almost exactly.

Delta CS and Gamma CS 45° give machine-level precision. Finite
difference methods are okay but depend heavily on h.

**Scenario 2 (Near-expiry + Low Volatility)**

**Parameters:** S = K = 100 r = q = 0 σ = 0.01 T = 1 / 365 (year)

  **h_rel**   **Δ_fd**   **Δ_cs**   **err Δ_fd**   **err Δ_cs**   **Γ_fd**   **Γ_cs real**   **Γ_45°**   **err Γ_fd**   **err Γ_cs real**   **err Γ_45°**
  ----------- ---------- ---------- -------------- -------------- ---------- --------------- ----------- -------------- ------------------- ---------------
  1e-2        0.9791     0.5002     0.4790         0.00012        0.9582     15.173          0.707       6.66           7.55                6.91
  1e-4        0.5381     0.5001     0.0380         0              7.599      15.244          7.599       0.023          7.622               0.023
  1e-6        0.5005     0.5001     0.00038        0              7.6218     15.244          7.6218      1.8 × 10⁻⁶     7.622               2.3 × 10⁻⁶
  1e-8        0.5001     0.5001     3.8 × 10⁻⁶     0              7.617      15.248          7.6218      0.0048         7.626               1.7 × 10⁻⁷

(Analytic Δ = 0.500104, Γ = 7.621781)

**Observations**

> **Δ (Delta)**
>
> Finite-difference Δ_fd becomes unstable when h is large → massive
> error (about 0.48).
>
> For smaller h, Δ_fd improves, but still noisy.
>
> **Complex-step Δ_cs** remains exact (0 error) for all h.
>
> **Γ (Gamma)**
>
> **Finite-difference Γ_fd** has error 6 -- 7 for large h, drops to
> about 0 for medium h, then grows again when h is too small.
>
> **Real-part complex Γ_cs_real** is consistently wrong (about 7 error
> offset).
>
> **45° complex-step Γ_cs_45°** tracks analytic Γ closely (error about 0
> for 1e-6 -- 1e-8).

**Interpretation & Scenario Effect**

Near expiry + low volatility → option price surface becomes a **sharp
kink** around S = K.\
That makes Δ and Γ change very steeply:

  **Aspect**          **Behaviour**
  ------------------- ---------------------------------------------------------------------
  Finite-difference   Highly sensitive to tiny perturbations in S → large numerical noise
  Complex-step        Bypasses subtractive round-off → stays stable even for very small T
  Real CS Γ           Fails due to incorrect formula scaling

**Reasoning (Truncation vs Round-off)**

**Finite-difference**:\
As h → small, subtraction of nearly equal prices creates round-off
noise.\
As h → large, the difference is too truncation error.\
→ Produces a **U-shaped error curve**.

**Complex-step**:\
Adds an imaginary perturbation → no subtraction of nearly equal
numbers.\
Round-off vanishes → flat, constant tiny error.

**Recommendation**

  **Greek**   **Preferred Method**          **Step Size (h_rel)**   **Reason**
  ----------- ----------------------------- ----------------------- ----------------------------------------------
  Δ           Complex-Step (Δ_cs)           1 × 10⁻⁸                Accurate even for sharp near-expiry profiles
  Γ           Complex-Step 45° (Γ_cs_45°)   1 × 10⁻⁶                Most stable and robust across all h
  Backup      Finite-Diff (FD)              1 × 10⁻⁴ -- 1 × 10⁻⁶    Usable only for moderate h values

In the near-expiry stress case, finite-difference methods break down
because the option price changes too abruptly.\
The complex-step Δ and especially the 45° Γ remain accurate and stable
for all step sizes, making them the best practical choice.

<https://colab.research.google.com/drive/1JFA7ZONjqPX90OI7HIvkFF1vT_3jPNYY?usp=sharing>

![](media/image3.png){width="6.268055555555556in"
height="4.241666666666666in"}

**Observations:**

-   **Delta CS:** still perfectly accurate, flat line.

-   **Delta FD:** again, shows U-shape pattern, error grows when h is
    too small.

-   **Gamma FD:** errors vary more strongly here because Gamma itself
    becomes small and numerically sensitive near expiry.

-   **Gamma CS real:** still large constant error (unstable).

-   **Gamma CS 45°:** very close to analytic across all step sizes.\
    The **most stable** even under extreme conditions (low volatility,
    short maturity).

Even when the option price is tiny, **complex-step methods** remain
accurate. Finite differences fluctuate, and the **real-part complex
Gamma** fails.

**Summery DF**

  **Scenario**      **Method**       **Max Error**   **Median Error**   **p99 Error**
  ----------------- ---------------- --------------- ------------------ ---------------
  **ATM**           err_D\_fd        0.009872        0.000050           0.009579
                    err_D\_cs        0.000000        0.000000           0.000000
                    err_G\_fd        0.019848        0.000002           0.019252
                    err_G\_cs_real   0.022785        0.019848           0.022697
                    err_G\_cs_45     0.000004        0.000000           0.000003
  **Near-expiry**   err_D\_fd        0.479014        0.019186           0.465783
                    err_D\_cs        0.000116        0.000000           0.000113
                    err_G\_fd        6.663544        0.013931           6.464331
                    err_G\_cs_real   7.626466        7.621782           7.626325
                    err_G\_cs_45     6.914675        0.011561           6.707928

**ATM Summary:**\
All methods behave well in this smooth scenario; complex-step Δ and 45°
Γ give exact analytic results.\
Real-part CS Γ shows constant bias. Finite-difference performs decently
for mid-range h.

**Near-expiry Summary:**\
Finite differences explode due to steep payoff curvature and tiny option
prices.\
Complex-step Δ remains perfect; 45° CS Γ still best but shows slight
high-h noise; real-part CS Γ fails completely.

I have one hands-on project in my GitHub (Time Series) with stock market
data, you can see in
[GitHub](https://github.com/MaryamSaamin/Data-Analysis.git).
