import pytest
import sasktran2 as sk

try:
    import sasktran
except ImportError:
    sasktran = None


@pytest.mark.skipif(sasktran is None, reason="sasktran not installed")
def test_mie_simple_database():
    _ = sk.database.MieDatabase(
        sk.mie.distribution.LogNormalDistribution(),
        sk.mie.refractive.H2SO4(),
        [532, 1020],
        median_radius=[100, 200],
        mode_width=[1.5, 1.6, 1.7],
    )
