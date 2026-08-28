from __future__ import annotations

from types import SimpleNamespace

import sasktran2 as sk


def test_database_root_environment_overrides_user_config(monkeypatch, tmp_path):
    configured_root = tmp_path / "configured"
    environment_root = tmp_path / "environment"

    monkeypatch.setattr(
        sk.appconfig,
        "load_user_config",
        lambda: {"database_root": configured_root.as_posix()},
    )
    monkeypatch.setenv("SASKTRAN2_DATABASE_ROOT", environment_root.as_posix())

    assert sk.appconfig.database_root() == environment_root


def test_database_root_uses_user_config_without_environment(monkeypatch, tmp_path):
    configured_root = tmp_path / "configured"

    monkeypatch.setattr(
        sk.appconfig,
        "load_user_config",
        lambda: {"database_root": configured_root.as_posix()},
    )
    monkeypatch.delenv("SASKTRAN2_DATABASE_ROOT", raising=False)

    assert sk.appconfig.database_root() == configured_root


def test_empty_database_root_environment_uses_default(monkeypatch, tmp_path):
    default_data_root = tmp_path / "default-data"

    monkeypatch.setattr(sk.appconfig, "load_user_config", dict)
    monkeypatch.setattr(
        sk.appconfig, "APPDIRS", SimpleNamespace(user_data_dir=default_data_root)
    )
    monkeypatch.setenv("SASKTRAN2_DATABASE_ROOT", "")

    assert sk.appconfig.database_root() == default_data_root / "database"
