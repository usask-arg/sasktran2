from sasktran2._core_rust import PyYankovsky
import sasktran2 as sk
from sasktran2.photchem import actinic_flux


class Yankovsky:
    _model: PyYankovsky

    def __init__(self):
        self._model = PyYankovsky()

    def run(self):
        flux = actinic_flux()

        self._model.solve(flux)


if __name__ == "__main__":
    test = Yankovsky()

    test.run()

    pass