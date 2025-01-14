from __future__ import annotations

import sasktran2 as sk
from sasktran2.viewinggeo.base import ViewingGeometryContainer


class Engine:
    def __init__(
        self,
        config: sk.Config,
        model_geometry: sk.Geometry1D,
        viewing_geo: sk.ViewingGeometry,
    ):
        """
        An Engine is the main class that handles the radiative transfer calculation.  The calculation takes
        place in two components.

        First, upon construction of the Engine, the majority of the geometry information is computed and
        cached.

        The main calculation takes place when calling :py:meth:`~calculate_radiance` with an
        :py:class:`sasktran2.Atmosphere` object where the actual radiative transfer calculation
        is performed.

        Parameters
        ----------
        config : sk.Config
            Configuration object
        model_geometry : sk.Geometry1D
            Geometry for the model
        viewing_geo : sk.ViewingGeometry
            Viewing geometry
        """
        self._nstokes = config.num_stokes

        if self._nstokes == 1:
            self._engine = sk.EngineStokes_1(config, model_geometry, viewing_geo)
        elif self._nstokes == 3:
            self._engine = sk.EngineStokes_3(config, model_geometry, viewing_geo)

        self._model_geometry = model_geometry
        self._viewing_geometry = viewing_geo

    def calculate_radiance(self, atmosphere: sk.Atmosphere, output: sk.Output = None):
        """
        Performs the radiative transfer calculation for a given atmosphere

        Parameters
        ----------
        atmosphere : sk.Atmosphere
            _description_
        output : sk.Output, optional
            Optional abstract output type.  Can be specified to change the
            output format. If set to None, the default :py:class:`sasktran2.OutputIdeal` class
            is used, by default None

        Returns
        -------
        varies
            Exact return type depends upon what output is set to.  For the default
            :py:class:`sasktran2.OutputIdeal`, the output format is an `xarray` Dataset.
            See for example, :py:meth:`sasktran2.OutputIdeal.post_process`
        """
        engine_output = (
            sk.OutputDerivMapped(self._nstokes) if output is None else output
        )

        self._engine.calculate_radiance(
            atmosphere.internal_object(), engine_output.internal_output()
        )

        result = engine_output.post_process(
            atmosphere, self._model_geometry, self._viewing_geometry
        )

        if isinstance(self._viewing_geometry, ViewingGeometryContainer):
            result = self._viewing_geometry.add_geometry_to_radiance(result)

        return result
