#include "TourismData.h"

#include <mml/tools/DataLoader.h>
#include <mml/tools/ConsolePrinter.h>
#include <mml/algorithms/Statistics.h>
#include <mml/algorithms/CurveFitting.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

using namespace MML;
using namespace Tourism;

namespace
{
    ColumnFormat NumericColumn(const std::string& name, int precision = 2)
    {
        return ColumnFormat(name).autoWidth().format(FormatType::Fixed).precision(precision).align(Alignment::Right);
    }

    TablePrinter<std::string, Real> ReportTable(const std::string& tag, std::vector<ColumnFormat> columns)
    {
        TablePrinter<std::string, Real> table(ColumnFormat(tag).autoWidth().align(Alignment::Left), std::move(columns));
        table.style(TableStyle().border(BorderStyle::Simple));
        return table;
    }

    Real Change(Real current, Real previous) { return 100 * (current / previous - 1); }
    Real CAGR(Real first, Real last, int years) { return 100 * (std::pow(last / first, Real{1} / years) - 1); }

    Real MonthShare(const TourismYear& year, std::initializer_list<int> months)
    {
        Real sum = 0;
        for (int month : months)
            sum += year.counts[Nights][month - 1];
        return 100 * sum / year.Total(Nights);
    }

    Real Concentration(const TourismYear& year)
    {
        Real sum = 0;
        for (Real nights : year.counts[Nights]) {
            const Real share = nights / year.Total(Nights);
            sum += share * share;
        }
        return sum;
    }

    void PrintAnnualAnalysis(const std::vector<TourismYear>& years)
    {
        std::cout << std::setprecision(2) << "\n=== Annual growth, visitor mix and stay length ===\n";
        TablePrinter<std::string, std::string> annual(ColumnFormat("Year").autoWidth().align(Alignment::Left),
            {NumericColumn("Arrivals", 0), NumericColumn("Nights", 0), NumericColumn("Arr.YoY%"),
             NumericColumn("Night.YoY%"), NumericColumn("Stay"), NumericColumn("Foreign nights%")});
        annual.style(TableStyle().border(BorderStyle::Simple));
        const TourismYear* previous = nullptr;
        for (const auto& year : years) {
            if (!year.Complete()) continue;
            annual.addRow(std::to_string(year.year), {
                NumericColumn("", 0).formatValue(year.Total(Arrivals)),
                NumericColumn("", 0).formatValue(year.Total(Nights)),
                previous ? NumericColumn("").formatValue(Change(year.Total(Arrivals), previous->Total(Arrivals))) : "n/a",
                previous ? NumericColumn("").formatValue(Change(year.Total(Nights), previous->Total(Nights))) : "n/a",
                NumericColumn("").formatValue(year.Total(Nights) / year.Total(Arrivals)),
                NumericColumn("").formatValue(100 * year.Total(ForeignNights) / year.Total(Nights))
            });
            previous = &year;
        }
        annual.print();
        const auto& first = years.front();
        const auto& last = years[20];
        std::cout << "\n2005-2025 growth and CAGR (20 intervals):\n";
        auto growthTable = ReportTable("Metric", {NumericColumn("Growth%"), NumericColumn("CAGR%")});
        for (auto metric : {Arrivals, Nights, DomesticArrivals, DomesticNights, ForeignArrivals, ForeignNights})
            growthTable.addRow(columnSuffixes[metric], {Change(last.Total(metric), first.Total(metric)),
                CAGR(first.Total(metric), last.Total(metric), 20)});
        growthTable.print();
        for (const auto* year : {&first, &last})
            std::cout << year->year << " nights/arrival: total " << year->Total(Nights) / year->Total(Arrivals)
                      << ", domestic " << year->Total(DomesticNights) / year->Total(DomesticArrivals)
                      << ", foreign " << year->Total(ForeignNights) / year->Total(ForeignArrivals)
                      << "; foreign arrivals share " << 100 * year->Total(ForeignArrivals) / year->Total(Arrivals) << "%\n";
        Real withinSegment = 0, visitorMix = 0;
        // Symmetric decomposition of S=sum(arrival share * segment stay length).
        for (auto metric : {DomesticArrivals, ForeignArrivals}) {
            const Real firstWeight = first.Total(metric) / first.Total(Arrivals);
            const Real lastWeight = last.Total(metric) / last.Total(Arrivals);
            const auto nightMetric = metric == DomesticArrivals ? DomesticNights : ForeignNights;
            const Real firstStay = first.Total(nightMetric) / first.Total(metric);
            const Real lastStay = last.Total(nightMetric) / last.Total(metric);
            withinSegment += (firstWeight + lastWeight) / 2 * (lastStay - firstStay);
            visitorMix += (firstStay + lastStay) / 2 * (lastWeight - firstWeight);
        }
        std::cout << std::setprecision(4) << "Stay-length change decomposition (nights/arrival): within-segment "
                  << withinSegment << ", visitor-mix " << visitorMix
                  << ", total " << withinSegment + visitorMix << "\n" << std::setprecision(2);
        std::cout << "Stay length uses ratios of totals, not unweighted monthly averages.\n";
    }

    void PrintSeasonality(const std::vector<TourismYear>& years)
    {
        std::cout << "\n=== Is the season getting longer? (overnight stays) ===\n";
        auto table = ReportTable("Year", {NumericColumn("Jul-Aug%"), NumericColumn("Jun-Sep%"),
            NumericColumn("Apr-May+Sep-Oct%"), NumericColumn("H", 4), NumericColumn("Effective months")});
        Vector<Real> time, peak, effective;
        for (const auto& year : years) {
            if (!year.Complete()) continue;
            const Real h = Concentration(year);
            table.addRow(std::to_string(year.year), {MonthShare(year, {7, 8}), MonthShare(year, {6, 7, 8, 9}),
                MonthShare(year, {4, 5, 9, 10}), h, 1 / h});
            time.push_back(year.year - 2005);
            peak.push_back(MonthShare(year, {7, 8}));
            effective.push_back(1 / h);
        }
        table.print();
        const auto peakFit = LinearLeastSquaresDetailed(time, peak);
        const auto effectiveFit = LinearLeastSquaresDetailed(time, effective);
        const auto& first = years.front();
        const auto& last = years[20];
        std::cout << "2005-2025 peak-share change: " << peak[20] - peak[0] << " percentage points\n"
                  << "Shoulder-share change: " << MonthShare(last, {4, 5, 9, 10}) - MonthShare(first, {4, 5, 9, 10})
                  << " percentage points\nEffective-month change: " << effective[20] - effective[0] << "\n"
                  << std::setprecision(4) << "OLS peak-share slope: " << peakFit.a << " percentage points/year"
                  << "; R^2=" << peakFit.r_squared << "\n"
                  << "OLS effective-month slope: " << effectiveFit.a << " months/year"
                  << "; R^2=" << effectiveFit.r_squared << "\n"
                  << "H=sum(monthly share^2); effective months=1/H, not literal season duration.\n"
                  << "These are descriptive trends including the COVID disruption, not significance tests.\n";
    }

    void PrintRecoveryAndVolatility(const std::vector<TourismYear>& years)
    {
        std::cout << "\n=== COVID shock, recovery and growth stability ===\n";
        const auto& baseline = years[14];
        const auto& covid = years[15];
        for (auto metric : {Arrivals, Nights, DomesticArrivals, DomesticNights, ForeignArrivals, ForeignNights}) {
            Vector<Real> growth, ordinaryGrowth;
            for (std::size_t i = 1; i <= 20; ++i) {
                growth.push_back(Change(years[i].Total(metric), years[i - 1].Total(metric)));
                if (years[i].year < 2020 || years[i].year > 2022)
                    ordinaryGrowth.push_back(growth[growth.size() - 1]);
            }
            int recovery = 0;
            for (std::size_t i = 16; i <= 20; ++i)
                if (!recovery && years[i].Total(metric) >= baseline.Total(metric))
                    recovery = years[i].year;
            const Real mean = Statistics::Avg(growth);
            const Real sd = Statistics::SampleStdDev(growth);
            std::cout << std::setprecision(2) << columnSuffixes[metric]
                      << ": 2020 vs 2019 " << Change(covid.Total(metric), baseline.Total(metric)) << "%"
                      << "; first >=2019: " << (recovery ? std::to_string(recovery) : "not yet")
                      << "; growth mean " << mean << "%, sample SD " << sd << " pp, CV ";
            if (std::abs(mean) < Real{0.01})
                std::cout << "n/a (near-zero mean)";
            else
                std::cout << sd / std::abs(mean);
            std::cout << "; SD excluding transitions ending 2020-2022: "
                      << Statistics::SampleStdDev(ordinaryGrowth) << " pp\n";
        }
        std::cout << "Growth CV=SD/abs(mean): sensitive to near-zero means and the rebound; not a risk forecast.\n";
        std::cout << "Seasonal night-growth stability (all transitions, then excluding 2020-2022):\n";
        auto stability = ReportTable("Season", {NumericColumn("Growth mean%"), NumericColumn("Sample SD (pp)"),
            NumericColumn("SD excl. 2020-22 (pp)")});
        const std::array<std::vector<int>, 3> seasons{{{7, 8}, {4, 5, 9, 10}, {1, 2, 3, 11, 12}}};
        const std::array<std::string, 3> names{"Peak Jul-Aug", "Shoulder Apr-May/Sep-Oct", "Off-season Jan-Mar/Nov-Dec"};
        for (std::size_t season = 0; season < seasons.size(); ++season) {
            Vector<Real> growth, ordinaryGrowth;
            auto total = [&](const TourismYear& year) {
                Real sum = 0;
                for (int month : seasons[season]) sum += year.counts[Nights][month - 1];
                return sum;
            };
            for (std::size_t i = 1; i <= 20; ++i) {
                growth.push_back(Change(total(years[i]), total(years[i - 1])));
                if (years[i].year < 2020 || years[i].year > 2022)
                    ordinaryGrowth.push_back(growth[growth.size() - 1]);
            }
            stability.addRow(names[season], {Statistics::Avg(growth), Statistics::SampleStdDev(growth),
                Statistics::SampleStdDev(ordinaryGrowth)});
        }
        stability.print();
        std::cout << "Annual nights relative to 2019 (2019=100):\n";
        auto indices = ReportTable("Year", {NumericColumn("Total"), NumericColumn("Domestic"), NumericColumn("Foreign")});
        for (std::size_t i = 15; i <= 20; ++i)
            indices.addRow(std::to_string(years[i].year), {100 * years[i].Total(Nights) / baseline.Total(Nights),
                100 * years[i].Total(DomesticNights) / baseline.Total(DomesticNights),
                100 * years[i].Total(ForeignNights) / baseline.Total(ForeignNights)});
        indices.print();
    }

    void PrintMonthlyGrowth(const std::vector<TourismYear>& years, const Data::Dataset& months)
    {
        std::cout << "\n=== Same-month growth and seasonality (nights) ===\n";
        auto table = ReportTable("Month", {NumericColumn("2005 share%"), NumericColumn("2025 share%"),
            NumericColumn("2005-25 CAGR%"), NumericColumn("2025 vs 2024%")});
        for (int month = 0; month < 12; ++month)
            table.addRow(months["Mjesec"].GetAsString(month), {
                100 * years[0].counts[Nights][month] / years[0].Total(Nights),
                100 * years[20].counts[Nights][month] / years[20].Total(Nights),
                CAGR(years[0].counts[Nights][month], years[20].counts[Nights][month], 20),
                Change(years[20].counts[Nights][month], years[19].counts[Nights][month])});
        table.print();
    }

    void PrintTrends(const std::vector<TourismYear>& years)
    {
        std::cout << "\n=== Annual linear trends and unusual years ===\n";
        Vector<Real> time, arrivalGrowth, nightGrowth;
        for (int i = 0; i <= 20; ++i) {
            time.push_back(i);
            if (i > 0) {
                arrivalGrowth.push_back(Change(years[i].Total(Arrivals), years[i - 1].Total(Arrivals)));
                nightGrowth.push_back(Change(years[i].Total(Nights), years[i - 1].Total(Nights)));
            }
        }
        for (auto metric : {Arrivals, Nights}) {
            Vector<Real> totals;
            for (int i = 0; i <= 20; ++i) totals.push_back(years[i].Total(metric));
            const auto fit = LinearLeastSquaresDetailed(time, totals);
            int highest = 0, lowest = 0;
            Vector<Real> residuals;
            for (int i = 0; i <= 20; ++i) {
                residuals.push_back(totals[i] - (fit.a * time[i] + fit.b));
                if (residuals[i] > residuals[highest]) highest = i;
                if (residuals[i] < residuals[lowest]) lowest = i;
            }
            std::cout << columnSuffixes[metric] << ": slope " << std::setprecision(2) << fit.a
                      << " counts/year; R^2=" << std::setprecision(4) << fit.r_squared << "\n"
                      << std::setprecision(0) << "  Largest positive residual: " << years[highest].year << " (" << residuals[highest]
                      << "); largest negative: " << years[lowest].year << " (" << residuals[lowest] << ")\n";
        }
        std::cout << std::setprecision(4) << "Correlation of annual arrival and night growth rates: "
                  << Statistics::PearsonCorrelation(arrivalGrowth, nightGrowth)
                  << "\nLinear fits include COVID, are descriptive, and are not forecasts.\n";
    }

    void PrintPartialYear(const std::vector<TourismYear>& years, const Data::Dataset& months)
    {
        const auto& current = years.back();
        const auto& previous = years[20];
        std::cout << "\n=== 2026 vs 2025: matched available months only ===\n";
        auto monthly = ReportTable("Month", {NumericColumn("Arrival change%"), NumericColumn("Night change%")});
        std::array<Real, 6> currentTotals{}, previousTotals{};
        int count = 0;
        for (int month = 0; month < 12; ++month) {
            if (!current.available[month]) continue;
            ++count;
            for (std::size_t metric = 0; metric < columnSuffixes.size(); ++metric) {
                currentTotals[metric] += current.counts[metric][month];
                previousTotals[metric] += previous.counts[metric][month];
            }
            monthly.addRow(months["Mjesec"].GetAsString(month), {
                Change(current.counts[Arrivals][month], previous.counts[Arrivals][month]),
                Change(current.counts[Nights][month], previous.counts[Nights][month])});
        }
        if (!count) {
            std::cout << "No available 2026 observations; no comparison calculated.\n";
            return;
        }
        monthly.print();
        std::cout << std::setprecision(2);
        std::cout << "Matched months: " << count << " (not a full-year comparison)\n";
        auto totals = ReportTable("Metric", {NumericColumn("2026", 0), NumericColumn("Matched 2025", 0),
            NumericColumn("Change%")});
        for (std::size_t metric = 0; metric < columnSuffixes.size(); ++metric)
            totals.addRow(columnSuffixes[metric], {currentTotals[metric], previousTotals[metric],
                Change(currentTotals[metric], previousTotals[metric])});
        totals.print();
        std::cout << "Matched-period nights/arrival: " << currentTotals[Nights] / currentTotals[Arrivals]
                  << " vs " << previousTotals[Nights] / previousTotals[Arrivals] << "\n";
    }
}

void PrintStatistics(const std::string& label, const Vector<Real>& values)
{
    Real mean, sampleStdDev;
    Statistics::AvgStdDev(values, mean, sampleStdDev);

    std::cout << "\n" << label << "\n";
    auto table = ReportTable("Statistic", {NumericColumn("Value")});
    table.addRow("Total", {std::accumulate(values.begin(), values.end(), Real{0})});
    table.addRow("Monthly mean", {mean});
    table.addRow("Monthly median", {Statistics::Median(values)});
    table.addRow("Minimum", {*std::min_element(values.begin(), values.end())});
    table.addRow("Maximum", {*std::max_element(values.begin(), values.end())});
    table.addRow("Sample std. dev.", {sampleStdDev});
    table.print();
}

int main(int argc, char* argv[])
{
    if (argc > 2) {
        std::cerr << "Usage: Example12_Statistics [CSV path]\n";
        return 1;
    }

    const std::string filename = argc == 2 ? argv[1] : "D:\\Projects\\MinimalMathLibrary\\data\\StatistikaTurizam.csv";

    try {
        Data::CSVLoadOptions options;
        options.delimiter = ';';
        options.strictSchema = true;
        // Keep source counts as text so the dataset's "...." marker is handled explicitly.
        options.inferTypes = false;
        const Data::Dataset data = Data::LoadCSV(filename, options);

        // The first row is an annual aggregate, not another monthly observation.
        const Data::Dataset months = data.FilterRows("Mjesec", [](const std::string& month) {
            return month != "Ukupno";
        });
        if (months.NumRows() != 12)
            throw DataError("Expected 12 monthly observations after excluding Ukupno");

        const auto years = ReadYears(months);
        const auto& arrivals = years[20].counts[Arrivals];
        const auto& nights = years[20].counts[Nights];

        std::cout << "Tourism statistics - 2025\n"
                  << "Source: " << filename << "\n"
                  << "Loaded " << data.NumRows() << " rows and " << data.NumColumns() << " columns.\n"
                  << "Analyzing 12 months; excluding the Ukupno totals row.\n\n";

        auto monthly = ReportTable("Month", {NumericColumn("Arrivals", 0), NumericColumn("Overnight stays", 0)});
        for (std::size_t row = 0; row < months.NumRows(); ++row)
            monthly.addRow(months["Mjesec"].GetAsString(row), {arrivals[row], nights[row]});
        monthly.print();

        std::cout << std::fixed << std::setprecision(2);
        PrintStatistics("Arrivals", arrivals);
        PrintStatistics("Overnight stays", nights);

        const Real totalArrivals = std::accumulate(arrivals.begin(), arrivals.end(), Real{0});
        const Real totalNights = std::accumulate(nights.begin(), nights.end(), Real{0});
        std::cout << "\nAverage nights per arrival: " << totalNights / totalArrivals
                  << "\nPearson correlation (monthly arrivals vs. nights): "
                  << std::setprecision(4) << Statistics::PearsonCorrelation(arrivals, nights)
                  << "\nCorrelation describes shared seasonality, not causation.\n";

        PrintAnnualAnalysis(years);
        PrintSeasonality(years);
        PrintRecoveryAndVolatility(years);
        std::cout << std::setprecision(2);
        PrintMonthlyGrowth(years, months);
        PrintTrends(years);
        PrintPartialYear(years, months);
        std::cout << "\nArrivals are visits, not unique people; nights are not revenue.\n"
                  << "Long-term comparisons assume consistent coverage and definitions in the source data.\n";
    }
    catch (const std::exception& error) {
        std::cerr << "Tourism statistics failed: " << error.what() << "\n";
        return 1;
    }

    return 0;
}
