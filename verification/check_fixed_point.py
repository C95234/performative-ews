"""
Verification INDEPENDANTE (reimplementation propre, pas une copie du script
fourni) de l'equation de point fixe du document performative_ews.tex :

    lambda_eff = lambda(t) + kappa * A( sigma^2 / (2*|lambda_eff|) )
    A(V) = 1 / (1 + exp(-k*(V - Vc)))

Objectifs :
1. Retrouver exactement les tableaux du document (kappa=-0.1, kappa=+0.05).
2. Scanner finement l'existence de racines pour kappa>0 pres du seuil de
   disparition annonce, pour verifier qu'il n'y a pas de bistabilite cachee
   (deuxieme racine que la bissection simple aurait pu manquer).
"""
import numpy as np

SIGMA = 1.0
K = 10.0
VC = 2.0


def alarm(V, k=K, vc=VC):
    return 1.0 / (1.0 + np.exp(-k * (V - vc)))


def residual(lam_eff, lam_t, kappa, sigma=SIGMA, k=K, vc=VC):
    V = sigma**2 / (2 * abs(lam_eff))
    return lam_eff - lam_t - kappa * alarm(V, k, vc)


def find_all_roots(lam_t, kappa, lo=-20.0, hi=-1e-6, n_scan=200000):
    """Balaye finement residual() sur [lo, hi] et renvoie TOUS les
    changements de signe (racines), pas seulement la premiere trouvee --
    detecte une bistabilite que la bissection simple manquerait."""
    grid = np.linspace(lo, hi, n_scan)
    vals = np.array([residual(x, lam_t, kappa) for x in grid])
    roots = []
    sign = np.sign(vals)
    for i in range(len(grid) - 1):
        if sign[i] == 0:
            roots.append(grid[i])
        elif sign[i] != sign[i + 1] and sign[i] != 0 and sign[i + 1] != 0:
            # bissection fine entre grid[i] et grid[i+1]
            a, b = grid[i], grid[i + 1]
            for _ in range(60):
                m = (a + b) / 2
                if np.sign(residual(a, lam_t, kappa)) == np.sign(residual(m, lam_t, kappa)):
                    a = m
                else:
                    b = m
            roots.append((a + b) / 2)
    return roots


def variance_at(lam_eff):
    return SIGMA**2 / (2 * abs(lam_eff))


print("=== Reproduction independante du tableau kappa=-0.1 (auto-invalidant) ===")
for lam_t in [-1.000, -0.100, -0.010, -0.001]:
    roots = find_all_roots(lam_t, kappa=-0.1)
    print(f"  lambda(t)={lam_t:+.3f} : {len(roots)} racine(s) -> {[f'{r:.4f}' for r in roots]}"
          f"  Var={[f'{variance_at(r):.2f}' for r in roots]}")

print()
print("=== Reproduction independante du tableau kappa=+0.05 (auto-realisateur) ===")
for lam_t in [-0.500, -0.100, -0.050, -0.020]:
    roots = find_all_roots(lam_t, kappa=0.05)
    label = [f"{r:.4f}" for r in roots] if roots else "AUCUNE SOLUTION"
    print(f"  lambda(t)={lam_t:+.3f} : {len(roots)} racine(s) -> {label}")

print()
print("=== Scan fin autour du seuil de disparition annonce (kappa=+0.05, entre -0.06 et -0.05) ===")
for lam_t in np.linspace(-0.06, -0.05, 21):
    roots = find_all_roots(lam_t, kappa=0.05)
    print(f"  lambda(t)={lam_t:+.5f} : {len(roots)} racine(s)")

print()
print("=== Balayage plus large de kappa>0 : lambda_c(kappa) (premier lambda(t) sans solution) ===")
for kappa in [0.01, 0.02, 0.05, 0.08, 0.10, 0.15]:
    lam_c = None
    for lam_t in np.linspace(-2.0, -0.0001, 4000):
        if not find_all_roots(lam_t, kappa, n_scan=4000):
            lam_c = lam_t
            break
    print(f"  kappa={kappa:+.2f} : lambda_c approx = {lam_c}")
