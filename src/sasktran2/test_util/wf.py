import numpy as np


def validate_wf(analytic, numerical, decimal=6):
    max_by_alt = np.abs(analytic).max(dim='altitude')

    max_by_alt.values[max_by_alt.values == 0] = 1e99

    rel_diff = (analytic - numerical) / max_by_alt

    np.testing.assert_array_almost_equal(rel_diff, 0, decimal=decimal)
