from __future__ import annotations

import importlib
import tarfile
import types

import pytest
from sasktran2.database import aer_line
from sasktran2.database.zenodo import (
    DOWNLOAD_RETRIES,
    DOWNLOAD_RETRY_PAUSE_S,
    download_zenodo_record,
)


def _fake_zenodo_get(monkeypatch, download):
    module = types.ModuleType("zenodo_get")
    module.download = download
    real_import_module = importlib.import_module

    def fake_import_module(name, package=None):
        if name == "zenodo_get":
            return module
        return real_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", fake_import_module)


@pytest.mark.parametrize(
    "error",
    [
        # zenodo-get 3: unreachable record metadata, and a file that failed
        # every attempt.
        ConnectionError("Timeout when fetching metadata for record 5120012"),
        RuntimeError("Download aborted for aer_v_3.8.1.tar.gz after 3 retries."),
    ],
)
def test_download_retries_then_raises_connection_error(monkeypatch, tmp_path, error):
    calls = []

    def download(
        record_or_doi, output_dir=".", retry_attempts=0, retry_pause=0.5, timeout=15.0
    ):
        calls.append((record_or_doi, output_dir, retry_attempts, retry_pause))
        raise error

    _fake_zenodo_get(monkeypatch, download)
    with pytest.raises(ConnectionError, match="Zenodo record 5120012") as caught:
        download_zenodo_record("5120012", tmp_path)
    assert caught.value.__cause__ is error
    assert calls == [
        ("5120012", tmp_path.as_posix(), DOWNLOAD_RETRIES, DOWNLOAD_RETRY_PAUSE_S)
    ]


def test_download_passes_only_supported_options(monkeypatch, tmp_path):
    received = {}

    def download(record_or_doi, output_dir="."):
        received.update(record=record_or_doi, output_dir=output_dir)

    _fake_zenodo_get(monkeypatch, download)
    download_zenodo_record("1485588", tmp_path)
    assert received == {"record": "1485588", "output_dir": tmp_path.as_posix()}


def test_missing_record_is_not_a_connection_error(monkeypatch, tmp_path):
    def download(record_or_doi, output_dir=".", **kwargs):
        msg = "HTTP error fetching metadata for record 0: 404 - NOT FOUND"
        raise ValueError(msg)

    _fake_zenodo_get(monkeypatch, download)
    with pytest.raises(ValueError, match="404"):
        download_zenodo_record("0", tmp_path)


def _aer_database(tmp_path, monkeypatch, write_tarball):
    pytest.importorskip("hapi")

    def fake_download(record_id, output_dir):
        write_tarball(output_dir.joinpath("aer_v_3.8.1.tar.gz"))

    monkeypatch.setattr(aer_line, "download_zenodo_record", fake_download)
    return aer_line.AERLineDatabase(db_root=tmp_path, rel_path=None)


def test_aer_database_is_extracted_into_place(monkeypatch, tmp_path):
    source = tmp_path.joinpath("source", "aer_v_3.8.1", "line_file")
    source.mkdir(parents=True)
    source.joinpath("aer_v_3.8.1_lines").write_text("lines\n")

    def write_tarball(path):
        with tarfile.open(path, "w:gz") as tar:
            tar.add(tmp_path.joinpath("source", "aer_v_3.8.1"), arcname="aer_v_3.8.1")

    db_root = tmp_path.joinpath("db")
    db = _aer_database(db_root, monkeypatch, write_tarball)
    path = db.path("")

    assert path == db_root.joinpath("aer_v_3.8.1")
    assert path.joinpath("line_file", "aer_v_3.8.1_lines").read_text() == "lines\n"
    # The tarball and the staging directory are gone.
    assert sorted(p.name for p in db_root.iterdir()) == ["aer_v_3.8.1"]


def test_aer_corrupt_download_leaves_no_partial_database(monkeypatch, tmp_path):
    def write_tarball(path):
        path.write_bytes(b"not a tarball")

    db = _aer_database(tmp_path, monkeypatch, write_tarball)
    with pytest.raises(tarfile.TarError):
        db.path("")
    # The next call downloads again instead of using a half-extracted database.
    assert not tmp_path.joinpath("aer_v_3.8.1").exists()
    assert not list(tmp_path.glob(".aer_extract_*"))
