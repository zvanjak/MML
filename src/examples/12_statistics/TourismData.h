#ifndef MML_EXAMPLE_TOURISM_DATA_H
#define MML_EXAMPLE_TOURISM_DATA_H

#include <mml/tools/data_loader/DataLoaderTypes.h>

#include <algorithm>
#include <array>
#include <numeric>
#include <string>
#include <vector>

namespace Tourism
{
    enum Metric { Arrivals, Nights, DomesticArrivals, DomesticNights, ForeignArrivals, ForeignNights };
    extern const std::array<std::string, 6> columnSuffixes;

    struct TourismYear {
        int year;
        std::array<MML::Vector<MML::Real>, 6> counts;
        std::array<bool, 12> available{};

        bool Complete() const {
            return std::all_of(available.begin(), available.end(), [](bool value) { return value; });
        }

        MML::Real Total(Metric metric) const {
            return std::accumulate(counts[metric].begin(), counts[metric].end(), MML::Real{0});
        }
    };

    std::vector<TourismYear> ReadYears(const MML::Data::Dataset& months);
}

#endif
