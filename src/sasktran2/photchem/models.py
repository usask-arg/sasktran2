from sasktran2._core_rust import PyYankovsky
import sasktran2 as sk
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
    state["O2(a)"].plot(y="altitude", label = "O2(a)")
    state["O2(a, v=1)"].plot(y="altitude", label = "O2(a, v=1)")
    state["O2(a, v=2)"].plot(y="altitude", label = "O2(a, v=2)")
    state["O2(a, v=3)"].plot(y="altitude", label = "O2(a, v=3)")
    state["O2(a, v=4)"].plot(y="altitude", label = "O2(a, v=4)")

    plt.xscale("log")
    plt.xlim(1e-3, 1e11)

    plt.ylim(60000, 120000)

    plt.legend()
    plt.show()

    state["O2(b)"].plot(y="altitude", label = "O2(b)")
    state["O2(b, v=1)"].plot(y="altitude", label = "O2(b, v=1)")
    state["O2(b, v=2)"].plot(y="altitude", label = "O2(b, v=2)")

    plt.xscale("log")
    plt.xlim(1e-3, 1e11)

    plt.ylim(60000, 120000)

    plt.legend()
    plt.show()

    pass