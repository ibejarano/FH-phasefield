from dolfin import Function, VectorFunctionSpace, TrialFunction, TestFunction

class DisplacementField:
    def __init__(self, mesh, V=None):
        self.V = V if V is not None else VectorFunctionSpace(mesh, "CG", 1)
        self.new = Function(self.V, name="displacement")
        self.old = Function(self.V)
        self.temp = Function(self.V)

    def update(self):
        self.old.assign(self.new)

    def get_trialfunction(self):
        return TrialFunction(self.V)

    def get_testfunction(self):
        return TestFunction(self.V)

    def get(self):
        return self.new
    
