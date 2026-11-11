#include "TourismData.h"

#include <cmath>
#include <initializer_list>
#include <limits>
#include <utility>

using namespace MML;

namespace Tourism
{
    const std::array<std::string, 6> columnSuffixes{
        "ukupno dolasci", "ukupno no\u0107enja", "doma\u0107i dolasci",
        "doma\u0107i no\u0107enja", "strani dolasci", "strani no\u0107enja"
    };

    namespace
    {
        Real ReadCount(const Data::DataColumn& column, std::size_t row)
        {
            const std::string text = column.GetAsString(row);
            if (column.IsMissing(row) || text == "....")
                return std::numeric_limits<Real>::quiet_NaN();
            if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
                throw DataError("Invalid count in " + column.name + " at monthly row " + std::to_string(row + 1));
            const Real count = std::stod(text);
            if (!std::isfinite(count) || count <= 0)
                throw DataError("Expected a positive finite count in " + column.name);
            return count;
        }
    }

    std::vector<TourismYear> ReadYears(const Data::Dataset& months)
    {
        const std::array<std::string, 12> labels{
            "Sije\u010danj", "Velja\u010da", "O\u017eujak", "Travanj", "Svibanj", "Lipanj",
            "Srpanj", "Kolovoz", "Rujan", "Listopad", "Studeni", "Prosinac"
        };
        for (std::size_t month = 0; month < labels.size(); ++month) {
            const std::string label = months["Mjesec"].GetAsString(month);
            if (label != labels[month] && !(month == 2 && label == "O\u00fbujak"))
                throw DataError("Expected chronological month rows; unexpected label " + label);
        }

        std::vector<TourismYear> years;
        for (int year = 2005; year <= 2026; ++year) {
            TourismYear values{year, {}, {}};
            for (std::size_t metric = 0; metric < columnSuffixes.size(); ++metric) {
                const auto& column = months[std::to_string(year) + " " + columnSuffixes[metric]];
                values.counts[metric] = Vector<Real>(12);
                for (std::size_t month = 0; month < 12; ++month)
                    values.counts[metric][month] = ReadCount(column, month);
            }
            for (std::size_t month = 0; month < 12; ++month) {
                values.available[month] = std::isfinite(values.counts[Arrivals][month]);
                for (const auto& metric : values.counts)
                    if (std::isfinite(metric[month]) != values.available[month])
                        throw DataError("Inconsistent segment availability in " + std::to_string(year));
                if (values.available[month]) {
                    for (auto metric : {Arrivals, Nights}) {
                        const Real segments = values.counts[metric + 2][month] + values.counts[metric + 4][month];
                        if (segments != values.counts[metric][month])
                            throw DataError("Domestic plus foreign counts do not match total in " + std::to_string(year));
                    }
                }
            }
            if (year < 2026 && !values.Complete())
                throw DataError("Incomplete historical year: " + std::to_string(year));
            years.push_back(std::move(values));
        }
        return years;
    }
}
