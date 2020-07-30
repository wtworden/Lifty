




class AlgebraicNum:
    def __init__(self,sage_alg_number):

        self._algebraic = QQbar(sage_alg_number)
        self._interval = CIF(self._algebraic)
        self._prec = 53

    def prec(self):
        return self._prec

    def algebraic(self):
        return self._algebraic

    def alg(self):
        return self._algebraic

    def interval(self):
        return self._interval

    def cif(self):
        return self._interval

    def set_interval(self, prec):
        self._prec = prec
        CIF = ComplexIntervalField(self.prec())
        self._interval = self.alg().interval(CIF)

    def __repr__(self):
        return self.alg().__repr__()

    def __str__(self):
        return self.alg().__str__()

    def __eq__(self,other):
        return self.alg() == other.alg()

    def __ne__(self,other):
        return not self.alg() == other.alg()

    def __mul__(self,other):
        return AlgebraicNum(self.alg()*other.alg())

    def __add__(self,other):
        return AlgebraicNum(self.alg() + other.alg())

    def __sub__(self,other):
        return AlgebraicNum(self.alg() - other.alg())

    def __pow__(self,integer):
        return AlgebraicNum(self.alg()**integer)

    def __neg__(self):
        return AlgebraicNum(-self.alg())

    def overlaps(self,other):
        return self.cif().overlaps(other.cif())


