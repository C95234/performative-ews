"""
Reimplementation INDEPENDANTE (ecrite from scratch a partir des equations du
document, sans copier le script fourni) de la simulation stochastique
complete du modele EWS performatif -- contre-verification : si cette
implementation independante retrouve le meme phenomene qualitatif (taux de
bascule monotone en kappa), c'est une preuve plus solide que le mecanisme
est reel plutot qu'un artefact d'un script particulier.

Difference deliberee avec le script fourni : la variance glissante est mise
a jour a CHAQUE pas (pas tous les 20 pas), et les graines utilisees
(seed0=10000) ne recoupent jamais celles des scripts fournis (seed=42,
range(500) a partir de 42).

dy = lambda_eff(t) * y * dt + sigma * dW,  lambda_eff = lambda_bare + kappa*A(V)
V = variance glissante de y sur une fenetre de window_time unites de temps.
"""
import numpy as np
import time


def alarm(V, k=10.0, vc=2.0):
    return 1.0 / (1.0 + np.exp(-k * (V - vc)))


def simulate_batch(kappa, n_trials, seed0, lam0=-1.0, rate=0.002, sigma=1.0,
                    dt=0.01, t_max=1200.0, window_time=30.0, y_thresh=8.0,
                    k=10.0, vc=2.0):
    rng = np.random.default_rng(seed0)
    window_n = int(window_time / dt)
    n_steps = int(t_max / dt)

    y = np.zeros(n_trials)
    buf = np.zeros((n_trials, window_n))
    s1 = np.zeros(n_trials)  # somme
    s2 = np.zeros(n_trials)  # somme des carres
    pos = 0
    filled = 0
    alive = np.ones(n_trials, dtype=bool)
    tipped = np.zeros(n_trials, dtype=bool)
    tipped_at = np.full(n_trials, np.nan)

    for step in range(n_steps):
        t = step * dt
        lam_bare = lam0 + rate * t
        if lam_bare > 0.05:
            break

        # variance glissante O(1) mise a jour a CHAQUE pas (different du
        # script fourni, qui ne recalcule que tous les 20 pas)
        if filled >= window_n:
            n = window_n
            V = np.clip(s2 / n - (s1 / n) ** 2, 0, None)
        else:
            V = np.zeros(n_trials)
        lam_eff = lam_bare + kappa * alarm(V, k, vc)

        noise = rng.standard_normal(n_trials)
        y_new = y + dt * (lam_eff * y) + sigma * np.sqrt(dt) * noise
        y = np.where(alive, y_new, y)

        old = buf[:, pos]
        s1 += np.where(alive, y - old, 0.0)
        s2 += np.where(alive, y * y - old * old, 0.0)
        buf[:, pos] = np.where(alive, y, buf[:, pos])
        pos = (pos + 1) % window_n
        filled += 1

        newly_tipped = alive & (np.abs(y) > y_thresh)
        tipped_at = np.where(newly_tipped, t, tipped_at)
        tipped |= newly_tipped
        alive &= ~newly_tipped
        if not alive.any():
            break

    return tipped, tipped_at


if __name__ == "__main__":
    t0 = time.time()
    n_trials = 500
    SEED0 = 10_000  # jamais utilisee par les scripts fournis (seed=42 / range depuis 42)
    print(f"=== Contre-verification independante ({n_trials} graines, seed0={SEED0}, jamais utilisees par les scripts fournis) ===")
    print(f"{'kappa':>10} {'bascules':>12} {'taux':>8} {'t_moy':>10} {'t_sd':>8}")
    rates = {}
    for kappa in [0.0, -0.05, -0.1, -0.2, 0.02, 0.05, 0.08]:
        tipped, tipped_at = simulate_batch(kappa, n_trials, seed0=SEED0 + int(round(kappa * 1000)))
        n_tip = int(tipped.sum())
        rate_tip = n_tip / n_trials
        valid = tipped_at[tipped]
        t_mean = np.nanmean(valid) if len(valid) else float("nan")
        t_sd = np.nanstd(valid) if len(valid) else float("nan")
        rates[kappa] = rate_tip
        label = "0 (classique)" if kappa == 0 else f"{kappa:+.2f}"
        print(f"{label:>10} {n_tip:5d}/{n_trials:<5d} {rate_tip:7.1%} {t_mean:10.1f} {t_sd:8.1f}")
    print(f"\nTemps ecoule: {time.time()-t0:.1f}s")

    print()
    seq = [rates[k] for k in [-0.2, -0.1, -0.05, 0.0, 0.02, 0.05, 0.08]]
    strict = all(seq[i] < seq[i + 1] for i in range(len(seq) - 1))
    print(f"Monotonie stricte sur cette implementation INDEPENDANTE : {strict}")
    print(f"  sequence : {[f'{r:.1%}' for r in seq]}")
