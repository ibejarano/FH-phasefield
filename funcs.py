from ufl import Identity, tr, inner, dev, conditional, sym, grad, lt


def lame_lambda(E, nu):
    return E*nu / ((1+nu) * (1-2*nu))

def lame_mu(E, nu):
    return E / (2*(1+nu))


#lmbda = lame_lambda(E, nu)
#mu = lame_mu(E, nu)

#Gc = Gc_param(subdomains, Gc1, Gc2)
lmbda = lame_lambda(E, nu)
mu = lame_mu(E, nu)

def epsilon(u):
  return sym(grad(u))

def sigma_0(u):
  return 2.0*mu*epsilon(u)+lmbda*tr(epsilon(u))*Identity(len(u))

def sigma(u, phi):
  return (1-phi)**2 * sigma_0(u)

# Densidad de energia elastica eq (12)
def psi(u):
  return 0.5*(lmbda+mu)*(0.5*(tr(epsilon(u)) + abs(tr(epsilon(u)))))**2 + mu*inner(dev(epsilon(u)), dev(epsilon(u)))

# Parámetro de trayectoria eq (11)
def H(uold, unew, Hold):
  return conditional(lt(psi(uold), psi(unew)), psi(unew), Hold)