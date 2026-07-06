from __future__ import annotations

import importlib
import zipfile

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from sasktran2.database.o3 import O3BirkWagnerDatabase, O3SerdyuchenkoDatabase
from sasktran2.database.zenodo import download_zenodo_record


def test_o3_zenodo_optical_properties_are_exported():
    assert sk.optical.O3BirkWagner is not None
    assert sk.optical.O3Serdyuchenko is not None


def test_birk_wagner_database_conversion(tmp_path):
    raw_file = tmp_path.joinpath(
        "O3_UV_region_absorption_cross_section_database_13112018.zip"
    )

    with zipfile.ZipFile(raw_file, "w") as zf:
        for temperature, scale in [(193, 1.0), (293, 2.0)]:
            zf.writestr(
                f"ACS/O3_ACS_{temperature}K.asc",
                "\n".join(
                    [
                        f"O3_ACS_{temperature}K.asc",
                        "wavelength[nm] absorption cross section[cm^2/molecule] uncertainty[cm^2/molecule]",
                        f"300.0 {scale * 1e-18:.6e} 1e-20",
                        f"301.0 {scale * 2e-18:.6e} 1e-20",
                    ]
                ),
            )

    db = O3BirkWagnerDatabase(db_root=tmp_path, rel_path=None)
    output_file = db.path()

    with xr.open_dataset(output_file) as ds:
        np.testing.assert_allclose(ds["temperature_k"], [193.0, 293.0])
        np.testing.assert_allclose(ds["wavelength_nm"], [300.0, 301.0])
        np.testing.assert_allclose(
            ds["xs"].to_numpy(),
            [[1e-22, 2e-22], [2e-22, 4e-22]],
        )


def test_serdyuchenko_database_conversion(tmp_path):
    raw_file = tmp_path.joinpath("SerdyuchenkoGorshelev5digits_latest.dat")
    header = ["header"] * 45
    row_0 = "213.33 " + " ".join(f"{value:.1e}" for value in range(1, 12))
    row_1 = "213.34 " + " ".join(f"{value:.1e}" for value in range(2, 13))
    raw_file.write_text("\n".join([*header, row_0, row_1]))

    db = O3SerdyuchenkoDatabase(db_root=tmp_path, rel_path=None)
    output_file = db.path()

    with xr.open_dataset(output_file) as ds:
        np.testing.assert_allclose(
            ds["temperature_k"],
            [
                193.0,
                203.0,
                213.0,
                223.0,
                233.0,
                243.0,
                253.0,
                263.0,
                273.0,
                283.0,
                293.0,
            ],
        )
        np.testing.assert_allclose(ds["wavelength_nm"], [213.33, 213.34])
        np.testing.assert_allclose(ds["xs"].isel(temperature_k=0), [11e-4, 12e-4])
        np.testing.assert_allclose(ds["xs"].isel(temperature_k=-1), [1e-4, 2e-4])


def test_download_zenodo_record_missing_dependency(monkeypatch, tmp_path):
    real_import_module = importlib.import_module

    def fake_import_module(name, package=None):
        if name == "zenodo_get":
            msg = "missing test dependency"
            raise ImportError(msg)
        return real_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", fake_import_module)

    with pytest.raises(ImportError, match="zenodo_get is required"):
        download_zenodo_record("5793207", tmp_path)
