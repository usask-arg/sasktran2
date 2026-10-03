# Frame-aligned successive-orders interpolation Implementation Plan

> **Note:** The implementation deviated from this plan, and the design spec
> (`docs/specs/2026-10-02-so-frame-aligned-interpolation-design.md`) is
> authoritative. For example, aligned grids apply the pole-avoiding
> pre-rotation to every aligned Lebedev rule rather than only to the
> reduced-horizon outgoing rule; aligned ground grids carry an extra 0.08 rad
> tilt so that no node lies on the horizon, while every ground hemisphere
> follows main's horizon rule (#306: nodes within `1e-12 * |location|` are
> excluded but keep half their weight) instead of the zero/`1e-12` rejection
> described in Task 3; and cubic LOS interpolation falls back to linear
> weights unless every stencil column is sunlit and the absolute stencil
> weights sum to at most 2.

> **For agentic workers:** REQUIRED SUB-SKILL: Use usask-foundation:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reach today's 11-column Geometry2D successive-orders accuracy with about 7 columns by aligning each column's angular grid to its local solar frame and interpolating the observer-LOS source cubically in horizontal angle.

**Architecture:** `SourceGeometry1D` builds per-column angular grids, each equal to today's grid rotated by the column's local frame `[solar horizontal, up x solar horizontal, up]`. It also builds a second, cubic-mode 2D location interpolator, used only for observer-LOS compilation. The scattering assemblers share scalar bases per rotation class (altitude). They give the vector operator explicit per-column synthesis groups instead of assuming one global synthesis. A `Config` switch restores today's behaviour.

**Tech Stack:** C++17 (Eigen, Catch2, OpenMP), Rust FFI (`sasktran2-sys`, `sasktran2-rs`, PyO3 `sasktran2-py-ext`), Python (`sasktran2`, pytest), pixi.

Spec: `docs/specs/2026-10-02-so-frame-aligned-interpolation-design.md`.

## Global Constraints

- Config switch name: `successive_orders_legacy_interpolation`, type bool, default `False`.
- Default (non-legacy) behaviour applies only to `geometrytype::spherical`; plane-parallel and pseudo-spherical keep today's grids.
- Aligned sphere = `R(p) * legacy_rotation` applied to the legacy sphere, where `R(p) = [x, up.cross(x), up]`, `up = p.normalized()`, `x = solar_horizontal_reference(up, geometry)`, and `legacy_rotation` is `pole_avoiding_rotation()` for reduced-horizon interior outgoing spheres and identity otherwise.
- Cubic LOS interpolation: 4-point Lagrange in horizontal angle, stencil start `clamp(lower - 1, 0, n - 4)`, used only strictly inside `(grid[0], grid[n-1])` and only when `n >= 4`; otherwise today's linear/extend behaviour.
- Diffuse incoming rays and ground forcing stay bilinear.
- Legacy mode (`successive_orders_legacy_interpolation = True`) must reproduce today's numerics bitwise.
- Rebuild commands: Python extension `pixi run build`; C++ tests `cd cpp && pixi run cmake --build build --target sasktran2_tests -j 8`; run C++ SO tests `cd cpp/build/lib/tests && ./sasktran2_tests "[successive_orders]"`.
- Python tests: `pixi run python -m pytest <path> -q -p no:cacheprovider`.
- Before every commit: `pixi run pre-commit` must pass. Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## File Map

| File | Responsibility |
|---|---|
| `cpp/include/sasktran2/config.h`, `cpp/include/c_api/config.h`, `cpp/c_api/config.cpp` | C++ and C API switch |
| `rust/sasktran2-sys/src/bindings.rs`, `rust/sasktran2-rs/src/bindings/config.rs`, `rust/sasktran2-py-ext/src/config.rs` | Rust FFI plumbing |
| `src/sasktran2/config.py`, `src/sasktran2/orbital.py` | Python property and orbital structural signature |
| `cpp/lib/successive_orders/geometry.h/.cpp` | Settings field, `SourcePoint::angular_class`, aligned spheres, LOS interpolator |
| `cpp/lib/successive_orders/horizontal_interpolation.h` (new) | Cubic Lagrange weight function |
| `cpp/lib/successive_orders/source.cpp` | Config to settings mapping |
| `cpp/lib/successive_orders/scattering.h/.cpp` | Vector operator synthesis groups |
| `cpp/lib/successive_orders/scattering_assembler.h/.cpp` | Class-shared scalar bases; per-column vector bases and groups |
| `cpp/lib/tests/successive_orders/test_*.cpp` | C++ tests |
| `tests/config/test_config_basic.py`, `tests/engine/test_successive_orders_2d.py` | Python tests |
| `docs/sphinx/source/components/source_terms/successive_orders.md` | User docs |
| `tools/benchmarks/so_interpolation_benchmark.py` (new), `tools/benchmarks/README.md` | Benchmarks |

---

### Task 1: `successive_orders_legacy_interpolation` switch

**Files:**
- Modify: `cpp/include/sasktran2/config.h` (after `set_successive_orders_reduced_horizon_quadrature`; member after `m_successive_orders_reduced_horizon_quadrature`)
- Modify: `cpp/include/c_api/config.h:87-90`, `cpp/c_api/config.cpp:408-425`
- Modify: `rust/sasktran2-sys/src/bindings.rs:879-890`
- Modify: `rust/sasktran2-rs/src/bindings/config.rs:809-845`
- Modify: `rust/sasktran2-py-ext/src/config.rs:548-562`
- Modify: `src/sasktran2/config.py:488-501`, `src/sasktran2/orbital.py:143`
- Modify: `cpp/lib/successive_orders/geometry.h` (`SourceGeometrySettings`), `cpp/lib/successive_orders/source.cpp` (`initialize_config`)
- Test: `tests/config/test_config_basic.py::test_cpp_successive_orders_controls_round_trip_and_validation`

**Interfaces:**
- Produces: C++ `Config::successive_orders_legacy_interpolation() const -> bool`, `Config::set_successive_orders_legacy_interpolation(bool)`; C `sk_config_get/set_successive_orders_legacy_interpolation`; Rust `Config::successive_orders_legacy_interpolation() -> Result<bool>`, `with_successive_orders_legacy_interpolation(bool)`; Python `sk.Config.successive_orders_legacy_interpolation`; `SourceGeometrySettings::legacy_interpolation` (bool, default `false`).

- [ ] **Step 1: Failing test.** In `tests/config/test_config_basic.py::test_cpp_successive_orders_controls_round_trip_and_validation` add `assert config.successive_orders_legacy_interpolation is False` beside the reduced-horizon default assertion, `config.successive_orders_legacy_interpolation = True` beside the reduced-horizon assignment, and `assert config.successive_orders_legacy_interpolation is True` beside its post-assignment assertion.
- [ ] **Step 2:** Run `pixi run python -m pytest tests/config/test_config_basic.py -q -p no:cacheprovider -k successive_orders_controls`; expect `AttributeError`.
- [ ] **Step 3: C++ Config.** Add to `config.h` after the reduced-horizon setter:

```cpp
        /** Restore the legacy successive-orders source interpolation: one
         * globally oriented angular grid shared by every source point and
         * bilinear observer-LOS interpolation. By default spherical geometries
         * align each source column's angular grid with its local solar frame
         * and Geometry2D interpolates the observer-LOS source cubically in
         * horizontal angle. */
        bool successive_orders_legacy_interpolation() const {
            return m_successive_orders_legacy_interpolation;
        }
        void set_successive_orders_legacy_interpolation(bool enabled) {
            m_successive_orders_legacy_interpolation = enabled;
        }
```

and the member `bool m_successive_orders_legacy_interpolation = false;` after `m_successive_orders_reduced_horizon_quadrature`.

- [ ] **Step 4: C API.** In `cpp/include/c_api/config.h` after the reduced-horizon declarations:

```c
int sk_config_get_successive_orders_legacy_interpolation(Config* config,
                                                         int* enabled);
int sk_config_set_successive_orders_legacy_interpolation(Config* config,
                                                         int enabled);
```

In `cpp/c_api/config.cpp` after the reduced-horizon definitions:

```cpp
int sk_config_get_successive_orders_legacy_interpolation(Config* config,
                                                         int* enabled) {
    if (config == nullptr || enabled == nullptr) {
        return -1;
    }
    *enabled = config->impl.successive_orders_legacy_interpolation() ? 1 : 0;
    return 0;
}

int sk_config_set_successive_orders_legacy_interpolation(Config* config,
                                                         int enabled) {
    if (config == nullptr) {
        return -1;
    }
    config->impl.set_successive_orders_legacy_interpolation(enabled != 0);
    return 0;
}
```

- [ ] **Step 5: Rust.** In `rust/sasktran2-sys/src/bindings.rs` after the reduced-horizon externs:

```rust
unsafe extern "C" {
    pub fn sk_config_get_successive_orders_legacy_interpolation(
        config: *mut Config,
        enabled: *mut ::std::os::raw::c_int,
    ) -> ::std::os::raw::c_int;
}
unsafe extern "C" {
    pub fn sk_config_set_successive_orders_legacy_interpolation(
        config: *mut Config,
        enabled: ::std::os::raw::c_int,
    ) -> ::std::os::raw::c_int;
}
```

In `rust/sasktran2-rs/src/bindings/config.rs` after `with_successive_orders_reduced_horizon_quadrature`:

```rust
    pub fn successive_orders_legacy_interpolation(&self) -> Result<bool> {
        let mut enabled = 0i32;
        let error_code = unsafe {
            ffi::sk_config_get_successive_orders_legacy_interpolation(self.config, &mut enabled)
        };
        if error_code != 0 {
            Err(anyhow!(
                "Error getting successive-orders legacy interpolation: error code {}",
                error_code
            ))
        } else {
            Ok(enabled != 0)
        }
    }

    pub fn with_successive_orders_legacy_interpolation(
        &mut self,
        enabled: bool,
    ) -> Result<&mut Self> {
        let error_code = unsafe {
            ffi::sk_config_set_successive_orders_legacy_interpolation(
                self.config,
                if enabled { 1 } else { 0 },
            )
        };
        if error_code != 0 {
            Err(anyhow!(
                "Error setting successive-orders legacy interpolation: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }
```

In `rust/sasktran2-py-ext/src/config.rs` after the reduced-horizon setter:

```rust
    #[getter]
    fn successive_orders_legacy_interpolation(&self) -> PyResult<bool> {
        self.config
            .successive_orders_legacy_interpolation()
            .into_pyresult()
    }

    #[setter]
    fn set_successive_orders_legacy_interpolation(&mut self, enabled: bool) -> PyResult<()> {
        self.config
            .with_successive_orders_legacy_interpolation(enabled)
            .into_pyresult()?;
        Ok(())
    }
```

- [ ] **Step 6: Python.** In `src/sasktran2/config.py` after the reduced-horizon property:

```python
    @property
    def successive_orders_legacy_interpolation(self) -> bool:
        """Restore the legacy successive-orders source interpolation.

        Disabled by default. Spherical successive-orders calculations then
        rotate each source column's angular grid into its local solar frame,
        so neighbouring columns sample identical local directions, and
        Geometry2D interpolates the observer line-of-sight source cubically in
        horizontal angle. Set this to ``True`` to use one globally oriented
        angular grid and bilinear line-of-sight interpolation.
        """
        return self._config.successive_orders_legacy_interpolation

    @successive_orders_legacy_interpolation.setter
    def successive_orders_legacy_interpolation(self, value: bool):
        self._config.successive_orders_legacy_interpolation = value
```

In `src/sasktran2/orbital.py::_structural_config_signature` add `"successive_orders_legacy_interpolation",` after `"successive_orders_reduced_horizon_quadrature",`.

- [ ] **Step 7: Settings mapping.** In `cpp/lib/successive_orders/geometry.h` `SourceGeometrySettings` add after `use_reduced_horizon_quadrature`:

```cpp
        /** Use one globally oriented angular grid and bilinear observer-LOS
         * interpolation instead of frame-aligned grids and cubic LOS. */
        bool legacy_interpolation = false;
```

In `cpp/lib/successive_orders/source.cpp::initialize_config` after the reduced-horizon assignment:

```cpp
            m_geometry_settings.legacy_interpolation =
                config.successive_orders_legacy_interpolation();
```

- [ ] **Step 8:** `pixi run build`, then rerun Step 2's command; expect PASS. Run `cargo clippy --all-targets --all-features`; expect no warnings.
- [ ] **Step 9: Commit** "Add successive_orders_legacy_interpolation config switch".

### Task 2: Vector scattering synthesis groups

**Files:**
- Modify: `cpp/lib/successive_orders/scattering.h:246-250,322` and `cpp/lib/successive_orders/scattering.cpp` (vector constructor near line 910, `apply` near line 1068)
- Modify: `cpp/lib/successive_orders/scattering_assembler.h` (vector private members), `cpp/lib/successive_orders/scattering_assembler.cpp:613-650,699-704`
- Test: `cpp/lib/tests/successive_orders/test_scattering.cpp`

**Interfaces:**
- Produces: `ScatteringOperator<3>(ScatteringBlockLayout, std::vector<std::shared_ptr<const VectorAngularBasis>> angular_bases, std::vector<int> synthesis_group_offsets = {})`. The offsets are CSR boundaries over atmospheric points: empty means no shared synthesis (per-point products); otherwise `front() == 0`, `back() == atmospheric_blocks()`, strictly increasing. Every point in `[offsets[g], offsets[g+1])` must have a basis whose synthesis equals that of `angular_bases[offsets[g]]`.
- Produces: `VectorScatteringAssembler` builds one basis per run of points sharing both sphere objects and one synthesis group per run of points sharing the outgoing sphere object.

- [ ] **Step 1: Failing test.** Append to `test_scattering.cpp`:

```cpp
namespace {
    class TestRotatedSphere final : public sasktran2::math::UnitSphere {
      public:
        TestRotatedSphere(int npoints, const Eigen::Matrix3d& rotation)
            : m_sphere(npoints), m_rotation(rotation) {}
        int num_points() const override { return m_sphere.num_points(); }
        Eigen::Vector3d get_quad_position(int index) const override {
            return m_rotation * m_sphere.get_quad_position(index);
        }
        double quadrature_weight(int index) const override {
            return m_sphere.quadrature_weight(index);
        }
        void interpolate(const Eigen::Vector3d& direction,
                         std::vector<std::pair<int, double>>& index_weights,
                         int& num_interp) const override {
            m_sphere.interpolate(m_rotation.transpose() * direction,
                                 index_weights, num_interp);
        }

      private:
        sasktran2::math::LebedevSphere m_sphere;
        Eigen::Matrix3d m_rotation;
    };
} // namespace

TEST_CASE("Vector successive-orders synthesis groups match per-point products",
          "[successive_orders][scattering][vector]") {
    using sasktran2::successive_orders::ScatteringBlockLayout;
    using sasktran2::successive_orders::ScatteringOperator;
    using sasktran2::successive_orders::VectorAngularBasis;
    constexpr int points = 4;
    constexpr int num_coefficients = 3;
    const sasktran2::math::LebedevSphere incoming(6);
    const sasktran2::math::LebedevSphere first_outgoing(14);
    const TestRotatedSphere second_outgoing(
        14, Eigen::AngleAxisd(0.3, Eigen::Vector3d::UnitY()).toRotationMatrix());
    std::vector<std::shared_ptr<const VectorAngularBasis>> bases;
    for (int point = 0; point < points; ++point) {
        const sasktran2::math::UnitSphere& outgoing =
            point < 2 ? static_cast<const sasktran2::math::UnitSphere&>(
                            first_outgoing)
                      : second_outgoing;
        bases.push_back(std::make_shared<const VectorAngularBasis>(
            incoming, outgoing, num_coefficients));
    }
    const auto make = [&](std::vector<int> offsets) {
        ScatteringOperator<3> scattering(
            ScatteringBlockLayout(points, 0, incoming.num_points(),
                                  first_outgoing.num_points(), 1, 2, 3),
            bases, std::move(offsets));
        Eigen::MatrixXd coefficients =
            Eigen::MatrixXd::Random(points, 4 * num_coefficients);
        coefficients.col(0).setOnes();
        scattering.set_atmospheric_coefficients(coefficients);
        return scattering;
    };
    std::srand(7);
    auto per_point = make({});
    std::srand(7);
    auto grouped = make({0, 2, 4});
    std::srand(7);
    auto single_group = make({0, 4});

    const Eigen::VectorXd incoming_values =
        Eigen::VectorXd::Random(per_point.input_size());
    Eigen::VectorXd expected(per_point.output_size());
    Eigen::VectorXd actual(per_point.output_size());
    Eigen::VectorXd wrong(per_point.output_size());
    auto workspace = per_point.make_workspace();
    per_point.apply(incoming_values, expected, workspace);
    grouped.apply(incoming_values, actual, workspace);
    single_group.apply(incoming_values, wrong, workspace);

    REQUIRE((actual - expected).norm() <= 1.0e-12 * expected.norm());
    REQUIRE((wrong - expected).norm() > 1.0e-6 * expected.norm());

    REQUIRE_THROWS_AS(make({0, 3}), std::invalid_argument);
    REQUIRE_THROWS_AS(make({0, 2, 2, 4}), std::invalid_argument);
}
```

- [ ] **Step 2:** Build C++ tests; expect a compile error (no vector-offset constructor).
- [ ] **Step 3: Operator.** In `scattering.h` replace the vector constructor's `bool point_bases_share_synthesis = false` with `std::vector<int> synthesis_group_offsets = {}` and the member `bool m_point_bases_share_synthesis = false;` with `std::vector<int> m_synthesis_group_offsets;`. In `scattering.cpp` change the vector constructor's signature and initializer accordingly (`m_synthesis_group_offsets(std::move(synthesis_group_offsets))`), and append this validation at the end of the constructor's existing checks:

```cpp
        if (!m_synthesis_group_offsets.empty()) {
            bool valid = m_synthesis_group_offsets.front() == 0 &&
                         m_synthesis_group_offsets.back() ==
                             m_layout.atmospheric_blocks();
            for (std::size_t group = 1;
                 valid && group < m_synthesis_group_offsets.size(); ++group) {
                valid = m_synthesis_group_offsets[group] >
                        m_synthesis_group_offsets[group - 1];
            }
            if (!valid) {
                throw std::invalid_argument(
                    "vector successive-orders synthesis groups must partition "
                    "the atmospheric points");
            }
        }
```

In `ScatteringOperator<3>::apply` replace `} else if (m_point_bases_share_synthesis) {` and its `m_basis->synthesize_active(...)` call with:

```cpp
            } else if (!m_synthesis_group_offsets.empty()) {
                const int active_modes =
                    m_active_coefficients * m_active_coefficients;
                for (auto& moments : workspace.m_angular.moments) {
                    moments.resize(m_layout.atmospheric_blocks(), active_modes);
                }
                for (int point = 0; point < m_layout.atmospheric_blocks();
                     ++point) {
                    m_point_bases[point]->analyze_active(
                        workspace.m_atmospheric_input.middleRows(point, 1),
                        m_active_coefficients, workspace.m_point_angular);
                    for (int moment = 0; moment < 6; ++moment) {
                        workspace.m_angular.moments[moment].row(point) =
                            workspace.m_point_angular.moments[moment].row(0);
                    }
                }
                m_basis->multiply_coefficients_active(
                    m_atmospheric_coefficients, m_active_coefficients,
                    workspace.m_angular);
                const int groups =
                    static_cast<int>(m_synthesis_group_offsets.size()) - 1;
                for (int group = 0; group < groups; ++group) {
                    const int begin = m_synthesis_group_offsets[group];
                    const int count =
                        m_synthesis_group_offsets[group + 1] - begin;
                    if (groups == 1) {
                        m_point_bases[begin]->synthesize_active(
                            workspace.m_atmospheric_output,
                            m_active_coefficients, workspace.m_angular);
                        continue;
                    }
                    for (int moment = 0; moment < 6; ++moment) {
                        workspace.m_point_angular.moments[moment] =
                            workspace.m_angular.moments[moment].middleRows(
                                begin, count);
                    }
                    m_point_bases[begin]->synthesize_active(
                        workspace.m_atmospheric_output.middleRows(begin, count),
                        m_active_coefficients, workspace.m_point_angular);
                }
                for (int point = 0; point < m_layout.atmospheric_blocks();
                     ++point) {
```

keeping the existing `add_frame_corrections_active` loop body that follows. A single group is the legacy shared path, unchanged: it synthesizes with `m_point_bases[0]`, which is `m_basis`.

- [ ] **Step 4: Assembler.** In `scattering_assembler.h` add `std::vector<int> m_synthesis_group_offsets;` to `VectorScatteringAssembler`'s private members. In `scattering_assembler.cpp` replace the `if (geometry.settings().use_reduced_horizon_quadrature) { ... }` block of the vector constructor with:

```cpp
        const auto& first_point = geometry.source_point(0);
        bool shares_grids = true;
        for (int point_index = 1;
             point_index < geometry.num_interior_points(); ++point_index) {
            const auto& point = geometry.source_point(point_index);
            shares_grids =
                shares_grids &&
                &point.incoming_sphere() == &first_point.incoming_sphere() &&
                &point.outgoing_sphere() == &first_point.outgoing_sphere();
        }
        if (!shares_grids) {
            // Q/U use global spin-harmonic frames, so points only share a
            // basis when they share both sphere objects, and synthesis only
            // when they share the outgoing sphere object.
            m_point_angular_bases.reserve(geometry.num_interior_points());
            m_synthesis_group_offsets.assign(1, 0);
            const sasktran2::math::UnitSphere* previous_incoming = nullptr;
            const sasktran2::math::UnitSphere* previous_outgoing = nullptr;
            std::shared_ptr<const VectorAngularBasis> basis;
            for (int point_index = 0;
                 point_index < geometry.num_interior_points(); ++point_index) {
                const auto& point = geometry.source_point(point_index);
                if (point_index != 0 &&
                    &point.outgoing_sphere() != previous_outgoing) {
                    m_synthesis_group_offsets.push_back(point_index);
                }
                if (point_index == 0) {
                    basis = m_angular_basis;
                } else if (&point.incoming_sphere() != previous_incoming ||
                           &point.outgoing_sphere() != previous_outgoing) {
                    basis = std::make_shared<const VectorAngularBasis>(
                        point.incoming_sphere(), point.outgoing_sphere(),
                        num_coefficients);
                }
                m_point_angular_bases.push_back(basis);
                previous_incoming = &point.incoming_sphere();
                previous_outgoing = &point.outgoing_sphere();
            }
            m_synthesis_group_offsets.push_back(
                geometry.num_interior_points());
        }
```

and in `create_operator` return `ScatteringOperator<3>(m_layout, m_point_angular_bases, m_synthesis_group_offsets);` for the point-basis case. With today's geometry this yields per-point bases and one group `{0, n}` for reduced horizon, identical to the previous `true`; plain Lebedev still uses the single basis.

- [ ] **Step 5:** Build C++ tests and run `./sasktran2_tests "[successive_orders]"`; expect all pass (124 cases).
- [ ] **Step 6:** `pixi run build`; run `pixi run python -m pytest tests/engine/test_successive_orders_2d.py tests/engine/test_successive_orders_cpp.py -q -p no:cacheprovider`; expect all pass (no behaviour change).
- [ ] **Step 7: Commit** "Partition vector successive-orders synthesis by outgoing grid".

### Task 3: Frame-aligned angular grids and class-shared scalar bases

**Files:**
- Modify: `cpp/lib/successive_orders/geometry.h` (`SourcePoint`)
- Modify: `cpp/lib/successive_orders/geometry.cpp` (`PoleAvoidingLebedevSphere`, `GroundUnitSphere`, `construct_source_points`)
- Modify: `cpp/lib/successive_orders/scattering_assembler.cpp:162-209` (scalar constructor)
- Test: `cpp/lib/tests/successive_orders/test_geometry.cpp`, `cpp/lib/tests/successive_orders/test_scattering_assembler.cpp`

**Interfaces:**
- Consumes: `SourceGeometrySettings::legacy_interpolation` (Task 1); vector synthesis groups (Task 2).
- Produces: `int SourcePoint::angular_class() const`. Interior points with equal class have scattering operators related by one rigid rotation applied to both of their grids. Aligned reduced-horizon gives class = altitude index; aligned or legacy Lebedev gives 0; legacy reduced-horizon gives the point index; ground points get -1.

- [ ] **Step 1: Failing geometry test.** Append to `test_geometry.cpp` inside the `#ifdef SKTRAN_RUST_SUPPORT` region used by the other Geometry2D tests:

```cpp
TEST_CASE("Successive-orders columns share frame-aligned outgoing grids",
          "[successive_orders][geometry][geometry2d]") {
    for (const bool reduced_horizon : {true, false}) {
        DYNAMIC_SECTION("reduced_horizon=" << reduced_horizon) {
            sasktran2::Geometry2D geometry(
                0.6, 0.4, 6372000.0, altitude_grid(), horizontal_angle_grid(),
                sasktran2::grids::interpolation::linear);
            sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
            const auto los = make_los_geometry(geometry, raytracer);
            sasktran2::successive_orders::SourceGeometrySettings settings;
            settings.num_incoming = 26;
            settings.num_outgoing = 26;
            settings.num_sza = 3;
            settings.num_threads = 1;
            settings.use_reduced_horizon_quadrature = reduced_horizon;
            sasktran2::successive_orders::SourceGeometry1D source(raytracer,
                                                                  geometry);
            source.initialize(los, settings);

            const int altitudes =
                static_cast<int>(source.source_altitudes_m().size());
            const auto frame = [&](const Eigen::Vector3d& position) {
                const Eigen::Vector3d up = position.normalized();
                const Eigen::Vector3d x =
                    sasktran2::successive_orders::solar_horizontal_reference(
                        up, geometry);
                Eigen::Matrix3d result;
                result << x, up.cross(x), up;
                return result;
            };
            const auto& reference = source.source_point(0);
            const Eigen::Matrix3d reference_frame =
                frame(reference.location().position);
            for (int index = 0; index < source.num_interior_points(); ++index) {
                const auto& point = source.source_point(index);
                const auto& column_first =
                    source.source_point(index - index % altitudes);
                REQUIRE(&point.outgoing_sphere() ==
                        &column_first.outgoing_sphere());
                REQUIRE(point.angular_class() ==
                        (reduced_horizon ? index % altitudes : 0));
                const Eigen::Matrix3d point_frame =
                    frame(point.location().position);
                for (int node = 0; node < point.num_outgoing(); ++node) {
                    REQUIRE((point_frame.transpose() *
                                 point.outgoing_sphere().get_quad_position(node) -
                             reference_frame.transpose() *
                                 reference.outgoing_sphere().get_quad_position(
                                     node))
                                .norm() < 1.0e-12);
                }
            }
            for (int ground = 0; ground < source.num_ground_points(); ++ground) {
                REQUIRE(source.source_point(source.num_interior_points() + ground)
                            .num_outgoing() ==
                        source.source_point(source.num_interior_points())
                            .num_outgoing());
            }
        }
    }
}

TEST_CASE("Successive-orders legacy interpolation keeps one global grid",
          "[successive_orders][geometry][geometry2d]") {
    sasktran2::Geometry2D geometry(0.6, 0.4, 6372000.0, altitude_grid(),
                                   horizontal_angle_grid(),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
    const auto los = make_los_geometry(geometry, raytracer);
    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 26;
    settings.num_outgoing = 26;
    settings.num_sza = 3;
    settings.num_threads = 1;
    settings.use_reduced_horizon_quadrature = true;
    settings.legacy_interpolation = true;
    sasktran2::successive_orders::SourceGeometry1D source(raytracer, geometry);
    source.initialize(los, settings);
    for (int index = 0; index < source.num_interior_points(); ++index) {
        REQUIRE(&source.source_point(index).outgoing_sphere() ==
                &source.source_point(0).outgoing_sphere());
        REQUIRE(source.source_point(index).angular_class() == index);
    }
}
```

- [ ] **Step 2: Failing assembler test.** Append to `test_scattering_assembler.cpp`, following that file's existing helpers for building a `SourceGeometry1D` from a `Geometry2D` (reuse its geometry factory if present, otherwise copy the setup from Step 1). For scalar (`ScalarScatteringAssembler`) and vector (`VectorScatteringAssembler`) in both quadrature modes:

```cpp
    auto assembled = assembler.create_operator();
    Eigen::MatrixXd coefficients = Eigen::MatrixXd::Random(
        assembled.atmospheric_coefficients().rows(),
        assembled.atmospheric_coefficients().cols());
    coefficients.col(0).setOnes();
    assembled.set_atmospheric_coefficients(coefficients);
    std::vector<std::shared_ptr<const BasisType>> exact;
    for (int point = 0; point < source.num_interior_points(); ++point) {
        exact.push_back(std::make_shared<const BasisType>(
            source.source_point(point).incoming_sphere(),
            source.source_point(point).outgoing_sphere(), num_coefficients));
    }
    ScatteringOperator<NSTOKES> reference(assembler.layout(), exact);
    reference.set_atmospheric_coefficients(coefficients);
    // Ground blocks are dense and identical in both; zero them.
    assembled.set_ground_values(Eigen::VectorXd::Zero(assembled.ground_value_size()));
    reference.set_ground_values(Eigen::VectorXd::Zero(reference.ground_value_size()));
    const Eigen::VectorXd input = Eigen::VectorXd::Random(assembled.input_size());
    Eigen::VectorXd expected(reference.output_size()), actual(assembled.output_size());
    auto workspace = reference.make_workspace();
    reference.apply(input, expected, workspace);
    assembled.apply(input, actual, workspace);
    REQUIRE((actual - expected).norm() <= 1.0e-11 * expected.norm());
```

The first argument of both operator constructors must be the layout the assembler exposes. For scalar also require that the number of distinct basis objects equals the number of source altitudes in reduced-horizon mode: count distinct `ScalarAngularBasis*` through `ScalarScatteringAssembler::storage_bytes()` decreasing relative to `legacy_interpolation = true`, i.e. `REQUIRE(aligned_assembler.storage_bytes() < legacy_assembler.storage_bytes())`.

- [ ] **Step 3:** Build C++ tests; expect compile failure (`angular_class` missing).
- [ ] **Step 4: SourcePoint.** In `geometry.h` add to `SourcePoint`'s public section `int angular_class() const { return m_angular_class; }` with the doc comment from Interfaces, and the private member `int m_angular_class = 0;`.
- [ ] **Step 5: Rotated spheres.** In `geometry.cpp`, replace class `PoleAvoidingLebedevSphere` with:

```cpp
        /** Rotation that moves Lebedev nodes off the coordinate poles, where
         * the meridian reference frame used by vector scattering is singular.
         * A rigid rotation preserves every weight and degree of exactness. */
        Eigen::Matrix3d pole_avoiding_rotation() {
            return (Eigen::AngleAxisd(0.01, Eigen::Vector3d::UnitZ()) *
                    Eigen::AngleAxisd(0.01, Eigen::Vector3d::UnitY()))
                .toRotationMatrix();
        }

        /** A Lebedev rule rigidly rotated by a fixed matrix. */
        class RotatedLebedevSphere final : public sasktran2::math::UnitSphere {
          public:
            RotatedLebedevSphere(int npoints, const Eigen::Matrix3d& rotation)
                : m_sphere(npoints), m_rotation(rotation) {}

            int num_points() const override { return m_sphere.num_points(); }

            Eigen::Vector3d get_quad_position(int index) const override {
                return m_rotation * m_sphere.get_quad_position(index);
            }

            double quadrature_weight(int index) const override {
                return m_sphere.quadrature_weight(index);
            }

            void interpolate(const Eigen::Vector3d& direction,
                             std::vector<std::pair<int, double>>& index_weights,
                             int& num_interp) const override {
                m_sphere.interpolate(m_rotation.transpose() * direction,
                                     index_weights, num_interp);
            }

          private:
            sasktran2::math::LebedevSphere m_sphere;
            Eigen::Matrix3d m_rotation;
        };
```

and, after `solar_horizontal_reference` is visible (it is declared in `interpolation.h`):

```cpp
        /** Rotation from canonical sphere axes to a point's local solar frame:
         * columns are the solar horizontal direction, up x that direction,
         * and the local vertical. Identity at a reference point whose solar
         * azimuth is zero. */
        Eigen::Matrix3d local_solar_frame(const Eigen::Vector3d& position,
                                          const sasktran2::Geometry& geometry) {
            const Eigen::Vector3d up = position.normalized();
            const Eigen::Vector3d x = solar_horizontal_reference(up, geometry);
            Eigen::Matrix3d frame;
            frame.col(0) = x;
            frame.col(1) = up.cross(x);
            frame.col(2) = up;
            return frame;
        }
```

In `GroundUnitSphere`'s constructor replace `if (m_full_sphere->get_quad_position(index).dot(location) > 0)` with:

```cpp
                // A frame-aligned rule places equator nodes exactly on the
                // horizon; reject them consistently rather than by roundoff.
                const double horizon_tolerance = 1.0e-12 * location.norm();
                ...
                    if (m_full_sphere->get_quad_position(index).dot(location) >
                        horizon_tolerance) {
```

(declare `horizon_tolerance` before the loop).

- [ ] **Step 6: Aligned construction.** In `construct_source_points` define after the ground-point loop that sets positions:

```cpp
        const bool aligned =
            !m_settings.legacy_interpolation &&
            m_geometry.coordinates().geometry_type() ==
                sasktran2::geometrytype::spherical;
        const int num_altitudes =
            static_cast<int>(m_source_altitudes_m.size());
        const auto column_frame = [&](int point_index) {
            return local_solar_frame(
                m_source_points[point_index].location().position, m_geometry);
        };
```

Make `volume_outgoing` in the reduced-horizon case `std::make_shared<const RotatedLebedevSphere>(m_settings.num_outgoing, pole_avoiding_rotation())` (bitwise identical to the old class). In the reduced-horizon interior loop replace `grid->outgoing = volume_outgoing;` with:

```cpp
                if (aligned) {
                    if (point_index % num_altitudes == 0) {
                        column_outgoing =
                            std::make_shared<const RotatedLebedevSphere>(
                                m_settings.num_outgoing,
                                column_frame(point_index) *
                                    pole_avoiding_rotation());
                    }
                    grid->outgoing = column_outgoing;
                    m_source_points[point_index].m_angular_class =
                        point_index % num_altitudes;
                } else {
                    grid->outgoing = volume_outgoing;
                    m_source_points[point_index].m_angular_class = point_index;
                }
```

declaring `std::shared_ptr<const sasktran2::math::UnitSphere> column_outgoing;` before the loop. In the Lebedev branch, when `aligned`, build one `AngularGridPair` per column, whose `incoming` and `outgoing` are `RotatedLebedevSphere(num, column_frame(first point of column))`, assign it to every point of the column and set `m_angular_class = 0`. Otherwise keep today's shared pair and set `m_angular_class = 0`. For ground points compute `const Eigen::Matrix3d ground_frame = aligned ? local_solar_frame(location, m_geometry) : Eigen::Matrix3d::Identity();`. Wrap `RotatedLebedevSphere(num, ground_frame)` instead of `LebedevSphere(num)` for the ground outgoing sphere, and for the ground incoming sphere in Lebedev mode, only when `aligned`. Keep the legacy `LebedevSphere` objects when not aligned, so legacy stays bitwise. Set `m_angular_class = -1` for ground points.

- [ ] **Step 7: Scalar assembler.** In `ScalarScatteringAssembler`'s constructor replace the reduced-horizon basis loop with:

```cpp
        if (geometry.settings().use_reduced_horizon_quadrature) {
            const auto& first_point = geometry.source_point(0);
            // Points of one angular class have grids related by one rigid
            // rotation, which leaves the scalar operator unchanged.
            std::unordered_map<int, std::shared_ptr<const ScalarAngularBasis>>
                class_bases;
            class_bases.emplace(first_point.angular_class(), m_angular_basis);
            m_point_angular_bases.reserve(geometry.num_interior_points());
            for (int point_index = 0;
                 point_index < geometry.num_interior_points(); ++point_index) {
                const auto& point = geometry.source_point(point_index);
                auto found = class_bases.find(point.angular_class());
                if (found == class_bases.end()) {
                    std::shared_ptr<const ScalarAngularBasis> basis;
                    if (&point.outgoing_sphere() ==
                        &first_point.outgoing_sphere()) {
                        basis = std::make_shared<const ScalarAngularBasis>(
                            point.incoming_sphere(), *m_angular_basis);
                    } else {
                        basis = std::make_shared<const ScalarAngularBasis>(
                            point.incoming_sphere(), point.outgoing_sphere(),
                            num_coefficients);
                    }
                    found = class_bases
                                .emplace(point.angular_class(), std::move(basis))
                                .first;
                }
                m_point_angular_bases.push_back(found->second);
            }
```

In the `SASKTRAN2_PROFILE_MEMORY` block, count `analysis_bytes` once per distinct basis pointer (use an `std::unordered_set<const ScalarAngularBasis*>`). Include `<unordered_map>`.

- [ ] **Step 8:** Build and run `./sasktran2_tests "[successive_orders]"`. New tests must pass. If existing tests pin legacy orientation (for example exact node directions or the interpolation of axial directions), set `settings.legacy_interpolation = true` in that test only when it explicitly tests legacy behaviour; otherwise update the expectation and explain why in a comment.
- [ ] **Step 9:** `pixi run build`, then run `pixi run python -m pytest tests/engine/test_successive_orders_2d.py tests/engine/test_successive_orders_cpp.py tests/engine/test_orbital_plane.py -q -p no:cacheprovider`. The adjoint and finite-difference tests (`test_2d_successive_orders_native_products_are_adjoint` for scalar and vector, both quadrature modes) must pass unchanged. For a value-pinned test that fails, regenerate its expected values only after confirming the change stays within the angular-discretization level (relative 1e-2 or less), and note the regeneration in the commit message.
- [ ] **Step 10: Commit** "Align successive-orders angular grids with each column's solar frame".

### Task 4: Cubic horizontal LOS interpolation

**Files:**
- Create: `cpp/lib/successive_orders/horizontal_interpolation.h`
- Modify: `cpp/lib/successive_orders/geometry.h` (member and accessor), `cpp/lib/successive_orders/geometry.cpp` (`AltitudeAngleSourceLocationInterpolator`, `initialize`, `compile_los_interpolation`)
- Test: `cpp/lib/tests/successive_orders/test_geometry.cpp`, `tests/engine/test_successive_orders_2d.py`

**Interfaces:**
- Produces: `enum class HorizontalInterpolation { linear, cubic };` and `void cubic_lagrange_weights(const Eigen::VectorXd& grid, double x, std::array<int, 4>& indices, std::array<double, 4>& weights)`. Precondition: `grid.size() >= 4` and `grid[0] < x < grid[grid.size() - 1]`.
- Produces: `SourceGeometry1D` compiles LOS rays with a cubic-mode interpolator when `!legacy_interpolation`, Geometry2D and at least 4 columns.

- [ ] **Step 1: Failing tests.** Append to `test_geometry.cpp`:

```cpp
#include "../../successive_orders/horizontal_interpolation.h"

TEST_CASE("Successive-orders cubic horizontal weights reproduce cubics",
          "[successive_orders][geometry]") {
    Eigen::VectorXd grid(6);
    grid << 0.0, 0.3, 0.7, 1.2, 2.0, 2.1;
    const auto cubic = [](double x) {
        return 1.0 + 2.0 * x - 0.5 * x * x + 0.3 * x * x * x;
    };
    std::array<int, 4> indices{};
    std::array<double, 4> weights{};
    for (double x = 0.01; x < 2.1; x += 0.0137) {
        sasktran2::successive_orders::cubic_lagrange_weights(grid, x, indices,
                                                             weights);
        double value = 0.0;
        double total = 0.0;
        for (int m = 0; m < 4; ++m) {
            REQUIRE(indices[m] == indices[0] + m);
            value += weights[m] * cubic(grid[indices[m]]);
            total += weights[m];
        }
        REQUIRE(indices[0] >= 0);
        REQUIRE(indices[3] <= 5);
        REQUIRE(value == Catch::Approx(cubic(x)).epsilon(1.0e-12));
        REQUIRE(total == Catch::Approx(1.0).epsilon(1.0e-13));
    }
    sasktran2::successive_orders::cubic_lagrange_weights(grid, 0.7, indices,
                                                         weights);
    REQUIRE(indices[0] == 1);
    REQUIRE(weights[1] == 1.0);
    REQUIRE(weights[0] == 0.0);
    REQUIRE(weights[2] == 0.0);
    REQUIRE(weights[3] == 0.0);
}

TEST_CASE("Successive-orders LOS uses a four-column stencil only on the LOS",
          "[successive_orders][geometry][geometry2d]") {
    Eigen::VectorXd horizontal(5);
    horizontal << -0.4, -0.2, 0.0, 0.2, 0.4;
    sasktran2::Geometry2D geometry(0.6, 0.0, 6372000.0, altitude_grid(),
                                   std::move(horizontal),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
    sasktran2::viewinggeometry::InternalViewingGeometry los;
    los.traced_rays.resize(1);
    sasktran2::viewinggeometry::ViewingRay ray;
    ray.observer.position =
        (geometry.coordinates().earth_radius() + 2500.0) *
        geometry.coordinates().unit_vector_from_angles(-0.33, 0.0);
    const Eigen::Vector3d target =
        (geometry.coordinates().earth_radius() + 2500.0) *
        geometry.coordinates().unit_vector_from_angles(0.33, 0.0);
    ray.look_away = (target - ray.observer.position).normalized();
    raytracer.trace_ray(ray, los.traced_rays.front());

    const auto column_counts = [&](bool legacy) {
        sasktran2::successive_orders::SourceGeometrySettings settings;
        settings.num_incoming = 14;
        settings.num_outgoing = 14;
        settings.num_sza = 5;
        settings.num_threads = 1;
        settings.legacy_interpolation = legacy;
        sasktran2::successive_orders::SourceGeometry1D source(raytracer,
                                                              geometry);
        source.initialize(los, settings);
        const int altitudes =
            static_cast<int>(source.source_altitudes_m().size());
        const auto& offsets = source.outgoing_point_offsets();
        const auto count = [&](const auto& interpolation, auto columns_for) {
            int maximum = 0;
            for (std::size_t ray_index = 0; ray_index < interpolation.size();
                 ++ray_index) {
                const auto columns = columns_for(ray_index).to_vector();
                const auto& compiled = interpolation[ray_index];
                for (std::size_t layer = 0; layer < compiled.layers.size();
                     ++layer) {
                    std::vector<int> touched;
                    for (const auto weight : compiled.source_for_layer(layer)) {
                        const int global = columns[weight.row_inner_index()];
                        const int point = static_cast<int>(
                            std::upper_bound(offsets.begin(), offsets.end(),
                                             global) -
                            offsets.begin() - 1);
                        if (point < source.num_interior_points()) {
                            touched.push_back(point / altitudes);
                        }
                    }
                    std::sort(touched.begin(), touched.end());
                    touched.erase(std::unique(touched.begin(), touched.end()),
                                  touched.end());
                    maximum = std::max(maximum,
                                       static_cast<int>(touched.size()));
                }
            }
            return maximum;
        };
        return std::make_pair(
            count(source.los_interpolation(),
                  [&](std::size_t r) {
                      return source.los_transport_columns_for_ray(r);
                  }),
            count(source.incoming_interpolation(), [&](std::size_t r) {
                return source.transport_columns_for_ray(r);
            }));
    };
    const auto [cubic_los, cubic_incoming] = column_counts(false);
    const auto [legacy_los, legacy_incoming] = column_counts(true);
    REQUIRE(cubic_los == 4);
    REQUIRE(legacy_los == 2);
    REQUIRE(cubic_incoming <= 2);
    REQUIRE(legacy_incoming <= 2);
}
```

If `source_for_layer` iteration does not yield objects with `row_inner_index()`, use the view's `visit` interface as `ray_transport.cpp` does. Append to `tests/engine/test_successive_orders_2d.py`:

```python
def _convergence_radiance(num_columns: int, *, legacy: bool) -> np.ndarray:
    altitudes = np.arange(0.0, 60_001.0, 2_000.0)
    horizontal = np.deg2rad(np.arange(-15.0, 15.01, 1.0))
    geometry = sk.Geometry2D(
        cos_sza=float(np.cos(np.deg2rad(60.0))),
        solar_azimuth=0.0,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=altitudes,
        horizontal_angle_grid_radians=horizontal,
    )
    config = successive_orders_config(single_scatter_source=sk.SingleScatterSource.NoSource)
    config.num_threads = 4
    config.num_successive_orders_incoming = 26
    config.num_successive_orders_outgoing = 26
    config.num_successive_orders_iterations = 30
    config.successive_orders_relative_tolerance = 1.0e-8
    config.successive_orders_altitude_grid_m = np.arange(1_000.0, 59_001.0, 4_000.0)
    config.num_sza = num_columns
    config.successive_orders_legacy_interpolation = legacy
    viewing = sk.ViewingGeometry()
    for tangent in (15_000.0, 25_000.0, 35_000.0):
        viewing.add_ray(
            sk.TangentAltitude(
                tangent_altitude_m=tangent,
                observer_altitude_m=600_000.0,
                horizontal_angle_radians=0.0,
                viewing_azimuth_radians=0.0,
            )
        )
    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=np.array([450.0]), calculate_derivatives=False)
    _, altitude = np.meshgrid(horizontal, altitudes, indexing="ij")
    atmosphere.storage.total_extinction[:, 0] = (1.2e-5 * np.exp(-altitude / 7_500.0)).ravel()
    atmosphere.storage.ssa[:, 0] = 0.999
    atmosphere.leg_coeff.a1[0] = 1.0
    atmosphere.leg_coeff.a1[2] = 0.5
    atmosphere.surface.albedo[:] = 0.3
    return sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere).radiance.values.ravel()


def test_2d_aligned_cubic_interpolation_needs_fewer_columns():
    aligned_reference = _convergence_radiance(31, legacy=False)
    legacy_reference = _convergence_radiance(31, legacy=True)
    aligned_error = np.max(np.abs(_convergence_radiance(7, legacy=False) / aligned_reference - 1))
    legacy_error = np.max(np.abs(_convergence_radiance(11, legacy=True) / legacy_reference - 1))
    assert aligned_error < legacy_error
```

Calibration: before finalizing, run this test's two errors once and record them in a comment. If the margin is under 2x, raise the comparison's strictness to the observed ratio halfway point; never loosen below `aligned_error < legacy_error`.

- [ ] **Step 2:** Build C++ tests; expect a compile failure (missing header).
- [ ] **Step 3: Weight function.** Create `cpp/lib/successive_orders/horizontal_interpolation.h`:

```cpp
#pragma once

#include <Eigen/Core>

#include <algorithm>
#include <array>

namespace sasktran2::successive_orders {
    /** Horizontal interpolation of the stored source between columns. */
    enum class HorizontalInterpolation { linear, cubic };

    /** Four-point Lagrange weights on a strictly increasing grid.
     *
     * Requires at least four nodes and grid[0] < x < grid[n-1]. The stencil
     * is centred on the interval containing x and shifted inward at the ends.
     * Weights reproduce cubic polynomials and are exactly one/zero at nodes.
     */
    inline void cubic_lagrange_weights(const Eigen::VectorXd& grid, double x,
                                       std::array<int, 4>& indices,
                                       std::array<double, 4>& weights) {
        const int size = static_cast<int>(grid.size());
        const int lower = static_cast<int>(
            std::upper_bound(grid.data(), grid.data() + size, x) - grid.data() -
            1);
        const int start = std::clamp(lower - 1, 0, size - 4);
        for (int m = 0; m < 4; ++m) {
            double weight = 1.0;
            for (int n = 0; n < 4; ++n) {
                if (n != m) {
                    weight *= (x - grid[start + n]) /
                              (grid[start + m] - grid[start + n]);
                }
            }
            indices[m] = start + m;
            weights[m] = weight;
        }
    }
} // namespace sasktran2::successive_orders
```

- [ ] **Step 4: Interpolator.** In `geometry.cpp` include `"horizontal_interpolation.h"`. Give `AltitudeAngleSourceLocationInterpolator` a final constructor parameter `HorizontalInterpolation horizontal_interpolation` stored as `const HorizontalInterpolation m_horizontal_interpolation;`, and add:

```cpp
            void horizontal_stencil(double angle, std::array<int, 4>& indices,
                                    std::array<double, 4>& weights,
                                    int& count) const {
                const auto& grid = m_horizontal_grid.grid();
                const auto size = grid.size();
                if (m_horizontal_interpolation ==
                        HorizontalInterpolation::cubic &&
                    size >= 4 && angle > grid[0] && angle < grid[size - 1]) {
                    cubic_lagrange_weights(grid, angle, indices, weights);
                    count = 4;
                    return;
                }
                std::array<int, 2> linear_indices{};
                std::array<double, 2> linear_weights{};
                m_horizontal_grid.calculate_interpolation_weights(
                    angle, linear_indices, linear_weights, count);
                for (int index = 0; index < count; ++index) {
                    indices[index] = linear_indices[index];
                    weights[index] = linear_weights[index];
                }
            }
```

Use it in `interior_interpolation_weights` and `ground_interpolation_weights` in place of `m_horizontal_grid.calculate_interpolation_weights`, widening their horizontal arrays to `std::array<int, 4>` / `std::array<double, 4>`.

- [ ] **Step 5: LOS interpolator.** In `geometry.h` add the private member `std::unique_ptr<sasktran2::grids::SourceLocationInterpolator> m_los_location_interpolator;`. Document it: "Observer-LOS interpolation, when it differs from diffuse-ray interpolation; null otherwise". In `initialize`, the 1D branch calls `m_los_location_interpolator.reset();`. The 2D branch passes `HorizontalInterpolation::linear` to the existing interpolator, then:

```cpp
            m_los_location_interpolator.reset();
            if (!m_settings.legacy_interpolation &&
                m_source_horizontal_angles_rad.size() >= 4) {
                m_los_location_interpolator =
                    std::make_unique<AltitudeAngleSourceLocationInterpolator>(
                        make_altitude_grid(), *m_geometry_2d,
                        m_settings.num_sza,
                        m_settings.horizontal_angle_grid_radians,
                        HorizontalInterpolation::cubic);
            }
```

In `compile_los_interpolation` use `auto& los_interpolator = m_los_location_interpolator ? *m_los_location_interpolator : *m_location_interpolator;` in place of `*m_location_interpolator`.

- [ ] **Step 6:** Build and run `./sasktran2_tests "[successive_orders]"`; expect all pass.
- [ ] **Step 7:** `pixi run build`; run `pixi run python -m pytest tests/engine/test_successive_orders_2d.py tests/engine/test_orbital_plane.py tests/engine/test_successive_orders_cpp.py -q -p no:cacheprovider`; expect all pass, including the new convergence test.
- [ ] **Step 8: Commit** "Interpolate the 2D successive-orders LOS source cubically".

### Task 5: Documentation and full verification

**Files:**
- Modify: `docs/sphinx/source/components/source_terms/successive_orders.md`

- [ ] **Step 1:** Add a "Source interpolation" subsection after "Structured 2D Geometry" covering four points:
  - each column's angular grid is aligned with its local solar frame, so neighbouring columns sample identical local directions;
  - Geometry2D interpolates the observer line-of-sight source cubically in horizontal angle when at least four columns exist, while diffuse rays stay bilinear;
  - this typically lets about 7 columns match the accuracy that previously needed 11;
  - `successive_orders_legacy_interpolation` restores the previous behaviour.

  Add `sasktran2.Config.successive_orders_legacy_interpolation` to the autosummary list.
- [ ] **Step 2:** Run `pixi run test` (full Python suite) and the full C++ suite `./sasktran2_tests`; expect all pass.
- [ ] **Step 3:** Run `pixi run pre-commit`; expect pass.
- [ ] **Step 4: Commit** "Document frame-aligned successive-orders interpolation".

### Task 6: Timing, memory and accuracy benchmarks

**Files:**
- Create: `tools/benchmarks/so_interpolation_benchmark.py`
- Modify: `tools/benchmarks/README.md`

- [ ] **Step 1:** Write `so_interpolation_benchmark.py`, a self-contained script with no external data. It sets up the standard Geometry2D limb scene:
  - tangents 10–50 km, observer 600 km, 350/525/750 nm;
  - Rayleigh, ozone and aerosol, albedo 0.3, 0.5° atmosphere grid over ±20°;
  - 25 source altitudes;
  - solar cases `fwd30, fwd60, back60, side60, oblique70, fwd75, back75, side80`.

  For each case and scheme (`legacy`, `aligned`) and column count (5, 7, 9, 11), it reports four numbers:
  - **(a) accuracy:** max relative MS error against that scheme's own 81-column reference at the same angular resolution, plus the legacy reference's difference from the aligned reference;
  - **(b) wall time:** engine construction plus radiance, with 8 threads;
  - **(c) wall time of radiance plus full Jacobian:** with `calculate_derivatives=True`;
  - **(d) peak RSS of a single-threaded subprocess run:** read from `resource.getrusage(RUSAGE_CHILDREN)`, plus the transport nonzeros and source-weight bytes from the `SASKTRAN2_PROFILE_MEMORY` records.

  It runs at 110 directions by default, writes JSON and a Markdown table, and exposes `--cases`, `--columns`, `--directions`, `--output`.
- [ ] **Step 2:** Run it at 110 directions. Also run one synthetic orbital-plane scene built like the fixtures in `tests/engine/test_orbital_plane.py`, with time, peak RSS and accuracy at 5/7/11 columns, legacy vs aligned.
- [ ] **Step 3:** Add a README section with the command lines and a summary of the results.
- [ ] **Step 4:** `pixi run pre-commit`, then commit "Add successive-orders interpolation benchmark".

## Self-Review

- Spec coverage:
  - frame-aligned grids → Task 3;
  - basis sharing → Tasks 2 and 3;
  - cubic LOS → Task 4;
  - config switch → Task 1;
  - docs → Task 5;
  - tests → each task;
  - benchmarks → Task 6;
  - plane-parallel unchanged → the `aligned` predicate in Task 3;
  - legacy bitwise → Task 3 keeps the legacy objects and Task 2 keeps the single-group path.
- Signature consistency: `SourcePoint::angular_class()`, `SourceGeometrySettings::legacy_interpolation`, `HorizontalInterpolation`, `cubic_lagrange_weights` and `ScatteringOperator<3>(layout, bases, offsets)` are used identically across tasks.
