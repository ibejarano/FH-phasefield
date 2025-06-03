import numpy as np

def read_from_json(file_path):
    """
    Lee un archivo JSON y devuelve su contenido.
    Args:
        file_path: Ruta al archivo JSON.
    Returns:
        dict: Contenido del archivo JSON.
    """
    import json
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data

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