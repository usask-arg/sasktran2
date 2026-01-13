from __future__ import annotations

from sasktran2._core_rust import lto_test


def test_lto():
    lto_test()

    assert True
