from dolfin import *
from mpi4py import MPI
import numpy as np
import json
import logging

logger = logging.getLogger(__name__)

def export_phi_to_csv(phi, mesh, output_dir):
    coordenadas = mesh.coordinates()
    valores_phi = phi.compute_vertex_values(mesh)

    # 1 Obtenemos los indices donde phi es igual a 1
    indices = np.where(np.greater_equal(valores_phi, 0.5))

    coordenadas_filtradas = coordenadas[indices]

    nombre_archivo = output_dir / "output_coords.csv"
    csv_header = "x,y" if mesh.geometry().dim() == 2 else "x,y,z"

    np.savetxt(
        nombre_archivo,
        coordenadas_filtradas,
        delimiter=",",
        header=csv_header,
        comments=''
    )

    logger.info(f"Coordenadas maximas X:{coordenadas_filtradas[:,0].max():.1e} Y:{coordenadas_filtradas[:,1].max():.1e}")


def fracture_length(phi, x1=-1, x2=1, y=0.0, npoints=5000, cutoff=0.6):
    """
    Aproxima el largo de la fractura evaluando phi sobre la línea y sumando los tramos donde phi < cutoff.
    Args:
        phi: Function de FEniCS (campo de fase)
        x1, x2: extremos de la línea (float)
        y: coordenada y (float, default 0.0)
        npoints: número de puntos a muestrear
        cutoff: umbral de phi para considerar fractura
    Returns:
        float: largo de fractura aproximado
    """
    xs = np.linspace(x1, x2, npoints)
    points = [(x, y) for x in xs]
    phi_vals = np.array([phi(point) for point in points])
    dx = abs(x2 - x1) / (npoints - 1)
    mask = phi_vals > cutoff
    length = np.sum(mask) * dx
    return length

def compute_opening_grad(uvec, phivec, grad_phivec, lelem, lfracmax, tol=0.01):
    #grad_phi.assign(project(grad(phit), WW))
    w_overx = list()
    TOL_PHI = tol

    xsn = np.arange(-lelem, -lfracmax/2, -lelem)
    xsp = np.arange(0, lfracmax/2, lelem)
    xs = np.concatenate([np.flip(xsn), xsp])


    for x in xs:
        ys = np.arange(0, 5e-2, lelem)
        wx = 0.0
        for y in ys:
            u_y = uvec(x, y)
            phi_y = phivec(x, y)
            grad_phi_y = grad_phivec(x, y)

            if phi_y < TOL_PHI:
                break
            wx += -u_y[1] * grad_phi_y[1] * lelem
            
        w_overx.append(wx)
    return xs, np.array(w_overx)

def compute_opening_cutoff(uvec, phivec, lelem, lfracmax, phi_val=0.5):
    xsn = np.arange(-lelem, -lfracmax/2, -lelem)
    xsp = np.arange(0, lfracmax/2, lelem)
    xs = np.concatenate([np.flip(xsn), xsp])
    w_overx = list()

    for x in xs:
        ys = np.flip(np.arange(0, 5e-2, lelem))
        wx_value = 0.0
        for y in ys:
            u_y = uvec(x, y)
            phi_y = phivec(x, y)

            if phi_y >= phi_val:
                wx_value = u_y[1]
                break

        w_overx.append(wx_value)
    return xs, np.array(w_overx)

def compute_opening_overtime(uvec, phivec, lelem, x= 0.0, phi_val=0.5):
    """
    Calcula la apertura de la fractura a lo largo del tiempo.
    Args:
        uvec: Función vectorial de desplazamiento
        phivec: Función de fase
        lelem: Longitud del elemento
        lfracmax: Longitud máxima de la fractura
        phi_val: Valor de fase para considerar apertura (default 0.5)
    Returns:
        w_overx: Apertura a lo largo de la fractura
    """
    y_plus = np.flip(np.arange(lelem, 10e-2, lelem))
    y_minus = np.arange(-10e-2, -lelem, lelem)

    w_plus_value = 0.0
    w_minus_value = 0.0
    for y in y_plus:
        u_y = uvec(x, y)
        phi_y = phivec(x, y)
        if phi_y >= phi_val:
            w_plus_value = u_y[1]
            break
            
    for y in y_minus:
        u_y = uvec(x, y)
        phi_y = phivec(x, y)
        if phi_y >= phi_val:
            w_minus_value = u_y[1]
            break

    return w_plus_value, w_minus_value

def read_data(fname, overrrides=None):
    with open(fname, "r") as f:
        data = json.load(f)


    if overrrides is not None:
        overrides = parse_overrides(overrrides)
        for k, v in overrides.items():
            if k in data:
                data[k] = v
            else:
                print(f"Warning: Override key '{k}' not found in data. Skipping.")
    return data

def parse_overrides(args):
    overrides = {}
    for arg in args:
        if "=" in arg:
            k, v = arg.split("=", 1)
            # Intenta convertir a float, si falla deja como string
            try:
                v_eval = eval(v, {}, {})
            except Exception:
                v_eval = v
            overrides[k] = v_eval
    return overrides

# ("stress_0"+"stress_4")/2 + sqrt((("stress_0"-"stress_4")/2)^2 + "stress_1"^2)

if __name__ == "__main__":
    # Ejemplo de uso
    import matplotlib.pyplot as plt
    u = lambda x, y: np.array([0, 1/(y**4)])
    phi = lambda x, y: 0.0 if abs(y) > 0.005 else 1


    yp, yn = compute_opening_overtime(u, phi, lelem=0.0001, phi_val=0.5)

    #plt.plot(yp,label="Apertura positiva")
    #plt.plot(yn, label="Apertura negativa")
    plt.plot(yp, 1/abs(yp), label="Apertura positiva")
    plt.plot(yn, 1/abs(yn), label="Apertura negativa")
    plt.ylim(-1, 100)
    plt.show()