from dolfin import *


class DisplacementSpace:
    def __init__(self, mesh):
        self.V = VectorFunctionSpace(mesh, "CG", 1)
        self.test = TestFunction(self.V)
        self.trial = TrialFunction(self.V)

    def get_test(self):
        return self.test
    
    def get_trial(self):
        return self.trial

class PhiSpace:
    def __init__(self, mesh):
        self.V = FunctionSpace(mesh, "CG", 1)
        self.test = TestFunction(self.V)
        self.trial = TrialFunction(self.V)

    def get_test(self):
        return self.test
    
    def get_trial(self):
        return self.trial