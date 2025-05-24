from dolfin import Function, FunctionSpace, TrialFunction, TestFunction, errornorm

class PhaseField:
    def __init__(self, mesh):
        self.mesh = mesh
        self.V = FunctionSpace(mesh, "CG", 1)
        self.new = Function(self.V, name="phi")
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
    
    def get_old(self):
        return self.old
    
    def get_error(self):
        return errornorm(self.new, self.old, norm_type='l2', mesh=self.mesh)