import abc
import sasktran2 as sk


class Output(abc.ABC):
    @abc.abstractmethod
    def internal_output(self):
        pass

    @abc.abstractmethod
    def post_process(self, atmo: sk.Atmosphere, geometry: sk.Geometry1D, viewing_geo: sk.ViewingGeometry):
        pass


class OutputIdeal(Output):
    def __init__(self, nstokes: int):
        self._nstokes = nstokes

        if nstokes == 1:
            self._output = sk.OutputIdealStokes_1()
        elif nstokes == 3:
            self._output = sk.OutputIdealStokes_3()
        else:
            msg = 'nstokes must be 1 or 3'
            raise ValueError(msg)

    def internal_output(self):
        return self._output

    def post_process(self, atmo: sk.Atmosphere, geometry: sk.Geometry1D, viewing_geo: sk.ViewingGeometry):
        return self._output
