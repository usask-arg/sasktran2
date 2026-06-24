#pragma once

#include <cstddef>

namespace sasktran2::rust::raytracer {
    struct RustTraceSummary;
    struct RustTraceLayer;

    class CppTraceResult;

    void prepare_trace_result(CppTraceResult& result,
                              const RustTraceSummary& summary);
    void set_trace_layer(CppTraceResult& result, std::size_t index,
                         const RustTraceLayer& layer);
} // namespace sasktran2::rust::raytracer
