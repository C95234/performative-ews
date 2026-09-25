"""
Verification de l'approximation quasi-stationnaire (AQS) : le document
suppose V_t ~ sigma^2/(2|lambda_eff|) a chaque instant (separation
d'echelles de temps entre la relaxation rapide de la variance et la derive
lente de lambda(t)). Teste ici directement, independamment du reste :

Methode -- variance d'ENSEMBLE (pas la fenetre glissante d'une seule
trajectoire, qui melange biais de mesure et validite de l'AQS) : simule
N_ENS realisations INDEPENDANTES en parallele, prend la variance a travers
l'ensemble a chaque instant t (estimateur non biaise de Var(y_t) au sens
statistique usuel), et compare a la prediction quasi-stationnaire
sigma^2/(2|lambda(t)|) (kappa=0, cas de base sans retroaction -- teste la
validite de l'AQS elle-meme, independamment de la retroaction performative).
"""
import numpy as np
import time


def simulate_ensemble_no_tipping(n_ens, seed0, lam0=-1.0, rate=0.002, sigma=1.0,
                                  dt=0.01, t_max=490.0, snapshot_every=2000):
    """N_ENS trajectoires independantes, SANS seuil de bascule (on veut la
    variance d'ensemble loin avant toute bascule, jusqu'a t_max=490 < 500
    -- le temps de bascule moyen documente -- pour rester dans un regime
    ou quasiment aucune realisation n'a encore bascule)."""
    rng = np.random.default_rng(seed0)
    n_steps = int(t_max / dt)
    y = np.zeros(n_ens)
    snapshots_t = []
    snapshots_var = []
    for step in range(n_steps):
        t = step * dt
        lam = lam0 + rate * t
        noise = rng.standard_normal(n_ens)
        y = y + dt * (lam * y) + sigma * np.sqrt(dt) * noise
        if step % snapshot_every == 0:
            snapshots_t.append(t)
            snapshots_var.append(np.var(y))
    return np.array(snapshots_t), np.array(snapshots_var)


if __name__ == "__main__":
    t0 = time.time()
    N_ENS = 20000
    print(f"=== Validite de l'AQS : variance d'ensemble vs prediction quasi-stationnaire (kappa=0, N={N_ENS}) ===")
    ts, var_emp = simulate_ensemble_no_tipping(N_ENS, seed0=555)
    lam0, rate = -1.0, 0.002
    lams = lam0 + rate * ts
    var_theory = 1.0 / (2 * np.abs(lams))

    print(f"{'t':>8} {'lambda(t)':>10} {'Var empirique':>15} {'Var AQS':>10} {'ecart relatif':>14}")
    for t, lam, ve, vt in zip(ts, lams, var_emp, var_theory):
        rel = abs(ve - vt) / vt
        flag = "  <-- ecart>20%" if rel > 0.20 else ""
        print(f"{t:8.1f} {lam:10.4f} {ve:15.4f} {vt:10.4f} {rel:13.1%}{flag}")

    # Premier instant ou l'ecart relatif depasse durablement 20%
    rel_all = np.abs(var_emp - var_theory) / var_theory
    breakdown_idx = None
    for i in range(len(rel_all)):
        if all(rel_all[i:i+3] > 0.20) if i + 3 <= len(rel_all) else rel_all[i] > 0.20:
            breakdown_idx = i
            break
    if breakdown_idx is not None:
        print(f"\nL'AQS decroche (ecart >20%, de facon soutenue) a partir de t={ts[breakdown_idx]:.1f}, "
              f"lambda(t)={lams[breakdown_idx]:.4f}")
    else:
        print("\nL'AQS reste dans les 20% d'ecart sur toute la fenetre testee (jusqu'a t=490).")

    print(f"\nTemps ecoule: {time.time()-t0:.1f}s")
