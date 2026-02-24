#pragma once

namespace sasktran2::twostream {
    enum class SourceType { ONLY_SOLAR, ONLY_THERMAL, SOLAR_THERMAL };

    // Meta programming
    template <SourceType source> constexpr int num_azimuth() {
        if constexpr (source == SourceType::ONLY_THERMAL) {
            return 1;
        } else {
            return 2;
        }
    }

    template <SourceType source> constexpr bool has_solar() {
        return source == SourceType::ONLY_SOLAR ||
               source == SourceType::SOLAR_THERMAL;
    }

    template <SourceType source> constexpr bool has_thermal() {
        return source == SourceType::ONLY_THERMAL ||
               source == SourceType::SOLAR_THERMAL;
    }
}; // namespace sasktran2::twostream
