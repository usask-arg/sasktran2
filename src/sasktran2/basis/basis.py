from sasktran2._core_rust import PyBasis

class Basis:
    def __init__(self, basis: PyBasis):
        self._basis = basis

    def lower_limit(self) -> float:
        return self._basis.lower_limit()
    
    def upper_limit(self) -> float:
        return self._basis.upper_limit()
    
    def evaluate(self, x: float) -> float:
        return self._basis.evaluate(x)
    
    def _internal_object(self) -> PyBasis:
        return self._basis


class Rectangle(Basis):
    def __init__(self, left: float, right: float):
        self._left = left
        self._right = right

        super().__init__(PyBasis.new_rectangle(left, right))

class Delta(Basis):
    def __init__(self, center: float):
        self._center = center

        super().__init__(PyBasis.new_delta(center))

class Gaussian(Basis):
    def __init__(self, center: float, stdev: float, max_stdev = 5):
        self._center = center
        self._stdev = stdev
        self._max_stdev = max_stdev

        super().__init__(PyBasis.new_gaussian(center, stdev, max_stdev))

class Triangle(Basis):
    def __init__(self, left: float, right: float, center: float):
        self._left = left
        self._right = right
        self._center = center

        super().__init__(PyBasis.new_triangle(left, right, center))
