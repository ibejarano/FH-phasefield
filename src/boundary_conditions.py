from dolfin import DirichletBC, SubDomain, Constant, Point # etc.

def create_crack_domain(center, l0, w0):
    class CrackDomain(SubDomain):
        def inside(self, x, on_boundary):
            # Usa center, l0, w0 pasados como argumentos
            return abs(x[0] - center[0]) <= l0 and abs(x[1] - center[1]) <= w0
    return CrackDomain()

def setup_boundary_conditions(V, W, boundary_markers, data):
    bcs_u = []
    bcs_phi = []
    bc_config = data.get("boundary_conditions", {}) # Leer config del JSON

    # Ejemplo: BCs de desplazamiento
    for marker_id, bc_data in bc_config.get("displacement", {}).items():
        value = bc_data.get("value", [0.0, 0.0]) # Valor por defecto
        print(marker_id, value)
        bc = DirichletBC(W, Constant(value), boundary_markers, int(marker_id))
        bcs_u.append(bc)

    # Ejemplo: BC de phi (grieta inicial)
    crack_config = bc_config.get("initial_crack", {})
    if crack_config:
        center = crack_config.get("center", [0.0, 0.0])
        l0 = crack_config.get("l0", data.get("linit", 0.0)) # Usar linit si no está especificado aquí
        w0 = crack_config.get("w0", data.get("h", 0.0))     # Usar h si no está especificado aquí
        crack_subdomain = create_crack_domain(center, l0, w0)
        bc_phi_crack = DirichletBC(V, Constant(1.0), crack_subdomain)
        bcs_phi.append(bc_phi_crack)

    # Añadir otras BCs de phi si se definen en el JSON...

    return bcs_u, bcs_phi

# --- En main.py ---
# from src.boundary_conditions import setup_boundary_conditions
# ...
# bc_u, bc_phi = setup_boundary_conditions(W, V, boundary_markers, data)
# ...
