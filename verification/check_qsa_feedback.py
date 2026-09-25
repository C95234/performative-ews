"""
Verification de l'AQS DANS LE CAS AVEC RETROACTION (kappa != 0) : compare
la variance d'ENSEMBLE mesuree sur la vraie dynamique stochastique complete
(chaque trajectoire reagit a SA PROPRE variance glissante, comme dans le
modele reel) a la prediction du point fixe quasi-stationnaire
lambda_eff(lambda_bare(t), kappa) resolu independamment
(check_fixed_point.py) -- c'est l'AQS telle qu'utilisee dans l'analyse du
point fixe du document, pas la version simplifiee sans retroaction.
"""
import numpy as np
import time
import sys
sys.path.insert(0, ".")
from check_fixed_point import find_all_roots  # noqa: E402


def alarm(V, k=10.0, vc=2.0):
    return 1.0 / (1.0 + np.exp(-k * (V - vc)))


def simulate_ensemble_with_feedback(kappa, n_ens, seed0, lam0=-1.0, rate=0.002,
                                     sigma=1.0, dt=0.01, t_max=490.0,
                                     window_time=30.0, snapshot_every=2000):
    rng = np.random.default_rng(seed0)
    window_n = int(window_time / dt)
    n_steps = int(t_max / dt)
    y = np.zeros(n_ens)
    buf = np.zeros((n_ens, window_n))
    s1 = np.zeros(n_ens)
    s2 = np.zeros(n_ens)
    pos, filled = 0, 0
    snapshots_t, snapshots_var = [], []

    for step in range(n_steps):
        t = step * dt
        lam_bare = lam0 + rate * t
        if filled >= window_n:
            n = window_n
            V = np.clip(s2 / n - (s1 / n) ** 2, 0, None)
        else:
            V = np.zeros(n_ens)
        lam_eff = lam_bare + kappa * alarm(V)

        noise = rng.standard_normal(n_ens)
        y = y + dt * (lam_eff * y) + sigma * np.sqrt(dt) * noise

        old = buf[:, pos]
        s1 += y - old
        s2 += y * y - old * old
        buf[:, pos] = y
        pos = (pos + 1) % window_n
        filled += 1

        if step % snapshot_every == 0:
            snapshots_t.append(t)
            snapshots_var.append(np.var(y))  # variance d'ENSEMBLE (verite terrain)

    return np.array(snapshots_t), np.array(snapshots_var)


if __name__ == "__main__":
    t0 = time.time()
    N_ENS = 20000
    lam0, rate = -1.0, 0.002

    for kappa in [-0.1, 0.05]:
        print(f"\n=== kappa={kappa:+.2f} : variance d'ensemble (dynamique reelle) vs AQS (point fixe) ===")
        ts, var_emp = simulate_ensemble_with_feedback(kappa, N_ENS, seed0=777 + int(kappa * 1000))
        lams_bare = lam0 + rate * ts

        print(f"{'t':>8} {'lambda_bare':>12} {'Var empirique':>15} {'Var AQS (pt fixe)':>18} {'ecart relatif':>14}")
        for t, lam_bare, ve in zip(ts, lams_bare, var_emp):
            roots = find_all_roots(lam_bare, kappa, n_scan=4000)
            if not roots:
                print(f"{t:8.1f} {lam_bare:12.4f} {ve:15.4f} {'PAS DE POINT FIXE':>18}")
                continue
            lam_eff_theory = roots[0]
            v_theory = 1.0 / (2 * abs(lam_eff_theory))
            rel = abs(ve - v_theory) / v_theory
            flag = "  <-- ecart>20%" if rel > 0.20 else ""
            print(f"{t:8.1f} {lam_bare:12.4f} {ve:15.4f} {v_theory:18.4f} {rel:13.1%}{flag}")

    print(f"\nTemps ecoule: {time.time()-t0:.1f}s")
