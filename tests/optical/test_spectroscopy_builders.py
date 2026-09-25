"""Offline regression checks for acquisition of auditable spectroscopy data."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
from scipy.constants import c, e, epsilon_0, m_e


def _builder(name):
    path = Path(__file__).resolve().parents[2] / "tools/spectroscopy" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _zinc_row():
    return {
        "ID_i": "030001.000001",
        "ID_k": "030001.000003",
        "Aki(s^-1)": "3.8e+04",
        "fik": "1.62e-03",
        "J_i": "0",
        "J_k": "1",
        "ritz_wl_vac(nm)": "307.679063",
    }


def test_audited_zinc_rate_restored_but_new_conflicts_remain_quarantined():
    module = _builder("build_atomic_metals")
    row = _zinc_row()
    assert module.strength_reconciliation(row) is not None
    assert not module.inconsistent_decay_strength(row)
    row["Aki(s^-1)"] = "3.8e+05"
    assert module.strength_reconciliation(row) is None
    row["fik"] = "1.62e-01"
    assert module.inconsistent_decay_strength(row)


def test_builder_derives_zinc_absorption_from_verified_a(tmp_path, monkeypatch):
    module = _builder("build_atomic_metals")
    row = _zinc_row()
    row.update({"Ei(cm-1)": "0", "Ek(cm-1)": "32501.399", "conf_i": "4s2", "Type": ""})
    monkeypatch.setattr(
        module,
        "rows",
        lambda path: (
            [row]
            if path.stem.endswith("lines")
            else [{"Level (cm-1)": "0", "J": "0", "Level ID": "030001.000001"}]
        ),
    )
    monkeypatch.setattr(
        module, "acquire", lambda *_: {"retrieved_utc": "2026-09-24", "url": "test"}
    )
    record = module.build_species("Zn_I", tmp_path)
    assert record["status"] == "available"
    assert record["quarantined_decay_rate_count"] == 0
    assert record["reconciled_decay_rate_count"] == 1
    with xr.open_dataset(tmp_path / "Zn_I.nc") as result:
        # Independent SI conversion, not the builder's rounded cgs coefficient.
        expected = (
            m_e * epsilon_0 * c / (2 * np.pi * e**2) * (307.679063e-9) ** 2 * 3.8e4 * 3
        )
        np.testing.assert_allclose(result.oscillator_strength, expected, rtol=3e-6)
        np.testing.assert_allclose(result.nist_oscillator_strength, 1.62e-3)
        np.testing.assert_allclose(result.upper_total_a_s, 3.8e4)
        assert result.strength_reconciled.item() == 1


def test_molecular_reduction_does_not_overstate_source_temperature_validity(tmp_path):
    module = _builder("build_exomol_metals")
    full_path = tmp_path / "full.nc"
    xr.Dataset(attrs={"temperature_max_k": 500.0}).to_netcdf(full_path)
    with pytest.raises(ValueError, match="source dataset validity"):
        module.reduce_for_mesosphere(full_path, tmp_path, 5000, 501)


def test_molecular_serialization_failure_preserves_previous_cache(
    tmp_path, monkeypatch
):
    module = _builder("build_exomol_metals")
    output = tmp_path / "AlO.nc"
    output.write_bytes(b"previous complete cache")

    def fail_write(_self, path, **_kwargs):
        Path(path).write_bytes(b"incomplete replacement")
        msg = "interrupted serialization"
        raise OSError(msg)

    monkeypatch.setattr(xr.Dataset, "to_netcdf", fail_write)
    with pytest.raises(OSError, match="interrupted"):
        module.write_netcdf(xr.Dataset(), output)
    assert output.read_bytes() == b"previous complete cache"
    assert not output.with_suffix(".nc.part").exists()
