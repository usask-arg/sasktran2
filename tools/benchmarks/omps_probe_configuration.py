"""Record candidate configuration from the copied OMPS benchmark process.

The runner passes overrides through the production ``model_kwargs`` interface.
This hook observes the completed forward-model construction; it does not alter
the processor, installed packages, engine configuration, or numerical arrays.
"""

from __future__ import annotations

import json
from pathlib import Path

FIELDS = (
    "successive_orders_transport_cache_wavelengths",
    "successive_orders_incoming_directions_by_altitude",
    "successive_orders_outgoing_directions_by_altitude",
    "successive_orders_altitude_grid_m",
    "successive_orders_horizontal_angle_grid_radians",
    "num_successive_orders_incoming",
    "num_successive_orders_outgoing",
    "num_successive_orders_iterations",
    "num_sza",
    "num_streams",
    "num_singlescatter_moments",
    "num_threads",
    "successive_orders_relative_tolerance",
    "successive_orders_absolute_tolerance",
    "successive_orders_reduced_horizon_quadrature",
)


def plain(value):
    if hasattr(value, "tolist"):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def install_configuration_recorder(destination: Path, overrides: dict) -> None:
    # Import the selected candidate only within the guarded child process.
    import sasktran2 as sk  # noqa: PLC0415
    from skretrieval.retrieval.processing import Retrieval  # noqa: PLC0415

    for name in overrides:
        # setattr on a Python Config could otherwise silently create an unused
        # instance attribute when the candidate has no native property.
        if not isinstance(getattr(sk.Config, name, None), property):
            msg = f"Candidate Config has no native property {name!r}"
            raise AttributeError(msg)
    original = Retrieval._construct_forward_model
    records = []

    def record(self):
        result = original(self)
        for name, model in result._forward_models.items():
            config = model._engine_config
            settings = {
                field: plain(getattr(config, field))
                for field in FIELDS
                if isinstance(getattr(type(config), field, None), property)
            }
            for field, requested in overrides.items():
                effective = settings[field]
                if field.endswith("_directions_by_altitude"):
                    effective = effective or []
                if effective != plain(requested):
                    msg = f"Effective {field} differs from requested override"
                    raise ValueError(msg)
            altitudes = settings.get("successive_orders_altitude_grid_m") or []
            for field in (
                "successive_orders_incoming_directions_by_altitude",
                "successive_orders_outgoing_directions_by_altitude",
            ):
                profile = settings.get(field) or []
                if profile and len(profile) != len(altitudes):
                    msg = f"{field} must match the actual source altitude grid"
                    raise ValueError(msg)
            engines = []
            for key, engine in model._engine.items():
                item = {"key": str(key), "class": type(engine).__name__}
                if hasattr(engine, "num_groups"):
                    item["occupied_time_groups"] = int(engine.num_groups)
                if hasattr(engine, "group_diagnostics"):
                    groups = engine.group_diagnostics
                    item["viewing_rays_by_group"] = [
                        {
                            "group_index": int(group["group_index"]),
                            "viewing_rays": len(group["observation_indices"]),
                        }
                        for group in groups
                    ]
                    item["viewing_rays"] = sum(
                        group["viewing_rays"] for group in item["viewing_rays_by_group"]
                    )
                engines.append(item)
            records.append(
                {
                    "forward_model": str(name),
                    "effective_config": settings,
                    "engines": engines,
                }
            )
        snapshot = {
            "requested_overrides": overrides,
            "records": records,
            "source_ray_count_evidence": (
                "Actual source/angular counts are recorded by the native "
                "source_angular_grid profile and copied into native_geometry.json."
            ),
        }
        (destination / "effective_native_config.json").write_text(
            json.dumps(snapshot, indent=2) + "\n"
        )
        parameters_path = destination / "parameters.json"
        parameters = json.loads(parameters_path.read_text())
        parameters["candidate_native_config"] = snapshot
        parameters_path.write_text(json.dumps(parameters, indent=2) + "\n")
        return result

    Retrieval._construct_forward_model = record
