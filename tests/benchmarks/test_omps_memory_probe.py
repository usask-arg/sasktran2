from __future__ import annotations

import hashlib
import runpy
import zipfile
from pathlib import Path

import pytest

PROBE = runpy.run_path(
    str(Path(__file__).resolve().parents[2] / "tools/benchmarks/omps_memory_probe.py")
)


@pytest.fixture()
def candidate_evidence(tmp_path):
    package = tmp_path / "package"
    files = {
        "sasktran2/__init__.py": b"version = 'candidate'\n",
        "sasktran2/_core_rust.abi3.so": b"candidate-native-binary",
    }
    wheel = tmp_path / "candidate.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        for name, data in files.items():
            path = package / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
            archive.writestr(name, data)
    identity = {
        "binary_sha256": hashlib.sha256(
            files["sasktran2/_core_rust.abi3.so"]
        ).hexdigest()
    }
    provenance = {
        "build_completed_utc": "2026-10-01T00:00:00+00:00",
        "binary_sha256": identity["binary_sha256"],
        "wheel_sha256": PROBE["sha256"](wheel),
    }
    return identity, package, provenance, wheel


def test_build_and_wheel_evidence_bind_to_imported_package(candidate_evidence):
    verified = PROBE["validate_candidate_provenance"](*candidate_evidence)
    assert all(value is True for value in verified.values())


@pytest.mark.parametrize("completed", [False, True])
def test_incomplete_or_wrong_build_is_rejected(candidate_evidence, completed):
    identity, package, provenance, wheel = candidate_evidence
    if completed:
        provenance["binary_sha256"] = "wrong-binary"
    else:
        provenance.pop("build_completed_utc")
    with pytest.raises(ValueError, match=r"completed build|runtime binary"):
        PROBE["validate_candidate_provenance"](identity, package, provenance, wheel)


def test_wrong_wheel_is_rejected(candidate_evidence):
    identity, package, provenance, wheel = candidate_evidence
    provenance["wheel_sha256"] = "wrong-wheel"
    with pytest.raises(ValueError, match="wheel does not match"):
        PROBE["validate_candidate_provenance"](identity, package, provenance, wheel)


def test_modified_python_package_is_rejected(candidate_evidence):
    identity, package, provenance, wheel = candidate_evidence
    (package / "sasktran2/__init__.py").write_text("version = 'different'\n")
    with pytest.raises(ValueError, match="package differs"):
        PROBE["validate_candidate_provenance"](identity, package, provenance, wheel)


def test_wheel_binary_is_checked_without_build_provenance(candidate_evidence):
    identity, package, _, wheel = candidate_evidence
    identity["binary_sha256"] = "different-runtime"
    with pytest.raises(ValueError, match="wheel does not contain"):
        PROBE["validate_candidate_provenance"](identity, package, None, wheel)


def test_reference_without_build_attestation_is_explicit():
    verified = PROBE["validate_candidate_provenance"]({}, None, None, None)
    assert all(value is None for value in verified.values())
