import numpy as np

def cm_to_u(r_cm, L_unit_cm=5.11e12):
    """
    da cm a unità interne.
    """
    return r_cm / L_unit_cm


def u_to_cm(r_u, L_unit_cm=5.11e12):
    """
    da unità interne a cm.
    """
    return r_u * L_unit_cm


# =========================================================
# Condizioni iniziali cartesiane al pericentro
# =========================================================

def two_body_ic_pericenter(m1, m2, a, e, G=1.0):
    """
    - moto nel piano xy
    - z = vz = 0
    - il vettore relativo iniziale è lungo +x
    - la velocità relativa iniziale è lungo +y
    - frame del centro di massa

    """
    if not (0.0 <= e < 1.0):
        raise ValueError("Per orbita bound serve 0 <= e < 1.")

    M = m1 + m2

    # Distanza al pericentro
    r_p = a * (1.0 - e)

    # Velocità al pericentro
    v_p = np.sqrt(G * M * (1.0 + e) / (a * (1.0 - e)))

    # Vettori relativi cartesiani
    r_rel = np.array([r_p, 0.0, 0.0], dtype=float)
    v_rel = np.array([0.0, v_p, 0.0], dtype=float)

    # Passaggio al frame del centro di massa
    r1 = -(m2 / M) * r_rel
    r2 = +(m1 / M) * r_rel

    v1 = -(m2 / M) * v_rel
    v2 = +(m1 / M) * v_rel

    # Check utili
    eps = 0.5 * np.dot(v_rel, v_rel) - G * M / np.linalg.norm(r_rel)
    eps_theory = -G * M / (2.0 * a)
    T = 2.0 * np.pi * np.sqrt(a**3 / (G * M))

    return {
        "m1": m1,
        "m2": m2,
        "M": M,
        "a": a,
        "e": e,
        "r_p": r_p,
        "v_p": v_p,
        "r1": r1,
        "r2": r2,
        "v1": v1,
        "v2": v2,
        "r_rel": r_rel,
        "v_rel": v_rel,
        "eps": eps,
        "eps_theory": eps_theory,
        "T": T,
    }



def print_ic(ic):
    print("=== PARAMETRI ===")
    print(f"m1 = {ic['m1']}")
    print(f"m2 = {ic['m2']}")
    print(f"M  = {ic['M']}")
    print(f"a  = {ic['a']}")
    print(f"e  = {ic['e']}")
    print()

    print("=== MOTO RELATIVO ===")
    print(f"r_rel = {ic['r_rel']}")
    print(f"v_rel = {ic['v_rel']}")
    print(f"r_p   = {ic['r_p']}")
    print(f"v_p   = {ic['v_p']}")
    print()

    print("=== PARTICELLA 1 ===")
    print(f"r1 = {ic['r1']}")
    print(f"v1 = {ic['v1']}")
    print()

    print("=== PARTICELLA 2 ===")
    print(f"r2 = {ic['r2']}")
    print(f"v2 = {ic['v2']}")
    print()

    print("=== CHECK ENERGIA ===")
    print(f"eps numerica = {ic['eps']}")
    print(f"eps teorica  = {ic['eps_theory']}")
    print(f"T            = {ic['T']}")



def write_ic_file(filename, ic, t0=0.0):
    
    x1, y1, z1 = ic["r1"]
    vx1, vy1, vz1 = ic["v1"]

    x2, y2, z2 = ic["r2"]
    vx2, vy2, vz2 = ic["v2"]

    with open(filename, "w") as f:
        f.write("2\n")
        f.write(f"{t0:.16e}\n")
        f.write(
            f"{ic['m1']:.16e} {x1:.16e} {y1:.16e} {z1:.16e} "
            f"{vx1:.16e} {vy1:.16e} {vz1:.16e}\n"
        )
        f.write(
            f"{ic['m2']:.16e} {x2:.16e} {y2:.16e} {z2:.16e} "
            f"{vx2:.16e} {vy2:.16e} {vz2:.16e}\n"
        )



if __name__ == "__main__":
    # -------------------------
    # INPUT
    # -------------------------
    m1 = 1.0e-6
    m2 = 1.0
    e = 0.0167


    a_cm = 1.5e13
    a = cm_to_u(a_cm)  

    G = 1.0
    t0 = 0.0


    ic = two_body_ic_pericenter(m1=m1, m2=m2, a=a, e=e, G=G)

    
    print_ic(ic)


    outname = "IC2b.txt"
    write_ic_file(outname, ic, t0=t0)

    print()
    print(f"File '{outname}' scritto correttamente.")
    print (2*np.pi*np.sqrt(a**3/m2))