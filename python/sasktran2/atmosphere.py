import sasktran2 as sk


class Atmosphere:
    def __init__(self, numwavel: int, model_geometry: sk.Geometry1D, config: sk.Config, nstokes: int = 1):
        self._nstokes = nstokes

        if nstokes == 1:
            self._atmosphere = sk.AtmosphereStokes_1(numwavel, model_geometry, config)
        elif nstokes == 3:
            self._atmosphere = sk.AtmosphereStokes_3(numwavel, model_geometry, config)

    @property
    def storage(self):
        return self._atmosphere.storage

