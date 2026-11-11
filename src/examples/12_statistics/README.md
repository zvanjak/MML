# Tourism statistics

This example loads [StatistikaTurizam.csv](../../../data/StatistikaTurizam.csv)
with MML's [DataLoader](../../../mml/tools/DataLoader.h) and computes descriptive
statistics using [Statistics](../../../mml/algorithms/Statistics.h).

The CSV is UTF-8 with a BOM and uses semicolons. The example excludes the `Ukupno`
annual totals row, analyzes the **21 complete years from 2005 through 2025**, and
compares available 2026 months with exactly the same months in 2025.

The example is split into [main.cpp](main.cpp) for calculations and reporting,
[TourismData.h](TourismData.h) for the year data model, and
[TourismData.cpp](TourismData.cpp) for count parsing and yearly data validation.
Report tables use MML's [ConsolePrinter](../../../mml/tools/ConsolePrinter.h)
`TablePrinter` with ASCII borders, automatic widths, and explicit per-column
precision. Count columns have zero decimal places, percentages and stay lengths
have two, and concentration has four. Descriptive-statistics values retain two
decimal places. Headings and interpretation notes remain text.

## Analyses

- The original 2025 monthly table and descriptive statistics: totals, means,
  medians, minima, maxima, sample standard deviations, and Pearson correlation.
- Annual arrivals and nights, year-over-year growth, total growth, and CAGR
  (20 annual intervals between 2005 and 2025).
- Domestic and foreign growth, visitor shares and average nights per arrival.
  Stay length is the ratio of totals, not an unweighted mean of monthly ratios.
- An exact symmetric decomposition of the change in overall stay length into
  within-segment stay changes and domestic/foreign arrival-mix changes.
- Seasonality: July–August share, June–September share, and shoulder-season share
  (April, May, September, October), all based on overnight stays.
- Monthly concentration `H = sum(p_m^2)`, with `p_m` the month's share of annual
  nights. Effective months are `1/H`; these are not literal season duration.
- Linear trends in peak share and effective months, including slope and R².
- The 2020 shock relative to 2019, first subsequent annual recovery to at least
  the 2019 level, and annual night indices relative to 2019.
- Annual growth volatility by visitor segment and season. Sample standard
  deviations use `n - 1`; growth SD is in percentage points. Segment growth CV
  is `SD / abs(mean)`, omitted for means smaller than 0.01 percentage points.
  A second SD excludes transitions ending in 2020–2022 to show how strongly
  the disruption/rebound affects the results.
- Same-month night growth: 2005–2025 CAGR and 2025 versus 2024.
- Annual linear fits using MML's curve-fitting API, their largest positive and
  negative residuals, and correlation of annual arrival/night growth rates.
- Matched-month 2026 comparisons for all six count series, without extrapolating
  missing months or comparing a partial year with a complete year.

## Reading the results

In the supplied data, July–August's share of nights falls from **63.00% to 54.50%**,
shoulder-season share increases from **21.80% to 25.92%**, and effective months
increase from **4.35 to 5.15**. This supports a broader distribution of nights
between the endpoints, but the linear trend fits are weak and include COVID.
It is not proof of a steadily lengthening season or a statistical significance test.

Average nights per arrival fall from **5.54 to 4.58**. Arrivals grow at **4.23% CAGR**
and nights at **3.24% CAGR**. Total nights first regain the 2019 level in **2023**;
arrivals do so in **2024**.

Correlations are descriptive, not evidence of causation. Annual fits include COVID
and are not forecasts. Growth CV is unstable near zero mean growth and does not
by itself establish resilience. Arrivals are visits, not unique people, and nights
are not revenue. Long-term comparisons assume consistent source coverage and
definitions; the example cannot establish that from the CSV alone.

## Data validation

DataLoader loads counts as text with type inference disabled, avoiding accidental
coercion. The source's `....` marker is handled explicitly as unavailable data.
Historical years must have all twelve chronological monthly rows. Each available
month must have all six positive finite counts, and domestic plus foreign counts
must equal total counts. Partial 2026 observations are selected by availability,
not by an assumed number of months.

Month labels are printed as supplied in the CSV, including any source typos.
Use a UTF-8-capable terminal to display Croatian month names.

From the repository root, using a double-precision build:

```powershell
cmake -B build -DMML_REAL_TYPE=double
cmake --build build --config Release --target Example12_Statistics --parallel
& .\build\src\examples\Release\Example12_Statistics.exe
```

An optional first argument supplies a different CSV path with the same schema:

```powershell
& .\build\src\examples\Release\Example12_Statistics.exe 'D:\Projects\MinimalMathLibrary\data\StatistikaTurizam.csv'
```

Expected 2025 totals: **20,698,963 arrivals** and **94,820,989 overnight stays**.
The default path is `D:\Projects\MinimalMathLibrary\data\StatistikaTurizam.csv`;
pass a path explicitly when using a different checkout.
Missing files, schema mismatches, reordered months, incomplete historical years,
inconsistent segment counts, and invalid counts produce an error message and a
nonzero exit code.
