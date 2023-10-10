import sasktran2 as sk


class Engine(object):
    def __init__(self, config: sk.Config, model_geometry: sk.Geometry1D, viewing_geo: sk.ViewingGeometry
                 ):
        self._nstokes = config.num_stokes

        if self._nstokes == 1:
            self._engine = sk.EngineStokes_1(config, model_geometry, viewing_geo)
        elif self._nstokes == 3:
            self._engine = sk.EngineStokes_3(config, model_geometry, viewing_geo)

        self._model_geometry = model_geometry
        self._viewing_geometry = viewing_geo

    def calculate_radiance(self, atmosphere: sk.Atmosphere, output: sk.Output = None):
        engine_output = sk.OutputIdeal(self._nstokes) if output is None else output

        self._engine.calculate_radiance(atmosphere.internal_object(), engine_output.internal_output())

        return engine_output.post_process(atmosphere, self._model_geometry, self._viewing_geometry)
