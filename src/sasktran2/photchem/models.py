from __future__ import annotations

from sasktran2._core_rust import PyYankovsky
from sasktran2.photchem import actinic_flux


class Yankovsky:
    _model: PyYankovsky

    def __init__(self):
        self._model = PyYankovsky()

    def run(self):
        flux = actinic_flux()

        return self._model.solve(flux)


if __name__ == "__main__":
    test = Yankovsky()

    state = test.run()

    import matplotlib.pyplot as plt

    plt.subplot(1, 4, 1)
    (state["O2(a)"] / 1e6).plot(y="altitude", label="O2(a)")
    (state["O2(a, v=1)"] / 1e6).plot(y="altitude", label="O2(a, v=1)")
    (state["O2(a, v=2)"] / 1e6).plot(y="altitude", label="O2(a, v=2)")
    (state["O2(a, v=3)"] / 1e6).plot(y="altitude", label="O2(a, v=3)")
    (state["O2(a, v=4)"] / 1e6).plot(y="altitude", label="O2(a, v=4)")

    plt.xscale("log")
    plt.xlim(1e-3, 1e11)

    plt.ylim(60000, 120000)

    plt.legend()
    plt.subplot(1, 4, 2)

    (state["O2(b)"] / 1e6).plot(y="altitude", label="O2(b)")
    (state["O2(b, v=1)"] / 1e6).plot(y="altitude", label="O2(b, v=1)")
    (state["O2(b, v=2)"] / 1e6).plot(y="altitude", label="O2(b, v=2)")

    plt.xscale("log")
    plt.xlim(1e-3, 1e11)

    plt.ylim(60000, 120000)

    plt.legend()

    plt.xlabel("Density (cm$^{-3}$)")
    plt.subplot(1, 4, 3)

    (state["O2(X, v=1)"] / 1e6).plot(y="altitude", label="O2(X, v=1)")
    (state["O2(X, v=2)"] / 1e6).plot(y="altitude", label="O2(X, v=2)")
    (state["O2(X, v=3)"] / 1e6).plot(y="altitude", label="O2(X, v=3)")
    (state["O2(X, v=4)"] / 1e6).plot(y="altitude", label="O2(X, v=4)")
    (state["O2(X, v=5)"] / 1e6).plot(y="altitude", label="O2(X, v=5)")

    plt.xscale("log")
    plt.xlim(1e-3, 1e11)

    plt.ylim(60000, 120000)

    plt.legend()

    plt.xlabel("Density (cm$^{-3}$)")

    plt.subplot(1, 4, 4)

    (state["O(1D)"] / 1e6).plot(y="altitude", label="O(1D)")

    plt.xscale("log")
    plt.xlim(1e-3, 1e5)

    plt.ylim(60000, 120000)

    plt.legend()

    plt.xlabel("Density (cm$^{-3}$)")

    plt.show()
