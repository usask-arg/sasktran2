import sasktran2 as sk


class Engine(object):
    def __init__(self, config: sk.Config, model_geometry: sk.Geometry1D, viewing_geo: sk.ViewingGeometry, nstokes: int = 1
                 ):
        self._nstokes = nstokes

        if nstokes == 1:
            self._engine = sk.EngineStokes_1(config, model_geometry, viewing_geo)
        elif nstokes == 3:
            self._engine = sk.EngineStokes_3(config, model_geometry, viewing_geo)

    def calculate_radiance(self, atmosphere: sk.Atmosphere, output: sk.Output = None):
        if output is None:
            engine_output = sk.OutputIdeal(self._nstokes)
        else:
            engine_output = output

        self._engine.calculate_radiance(atmosphere._atmosphere, engine_output.internal_output())

        return engine_output.post_process()
