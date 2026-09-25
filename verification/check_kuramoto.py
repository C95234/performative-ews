"""
Reimplementation INDEPENDANTE de l'extension Kuramoto du modele EWS
performatif -- ecrite from scratch a partir de la description du document
(pas une copie du script fourni), graines disjointes.

dtheta_i = (omega_i + K_eff(t) * r * sin(psi - theta_i)) dt + bruit
K_eff(t) = K0(t) + kappa * A(r_t)
r*e^{i*psi} = (1/N) * sum_j e^{i*theta_j}   (parametre d'ordre de Kuramoto)
"""
import numpy as np
import time


def alarm(r, k=10.0, rc=0.5):
    return 1.0 / (1.0 + np.exp(-k * (r - rc)))


def simulate_batch(kappa, n_trials, seed0, N=60, K0_start=0.5, K0_rate=0.003,
                    dt=0.05, t_max=500.0, r_tip=0.7, sigma_omega=0.5, sigma_noise=0.05):
    rng = np.random.default_rng(seed0)
    n_steps = int(t_max / dt)

    omega = rng.normal(0, sigma_omega, size=(n_trials, N))
    theta = rng.uniform(0, 2 * np.pi, size=(n_trials, N))
    tipped = np.zeros(n_trials, dtype=bool)
    tipped_at = np.full(n_trials, np.nan)
    alive = np.ones(n_trials, dtype=bool)

    for step in range(n_steps):
        t = step * dt
        K0 = K0_start + K0_rate * t

        c, s = np.cos(theta), np.sin(theta)
        rx, ry = c.mean(axis=1), s.mean(axis=1)
        r = np.hypot(rx, ry)
        psi = np.arctan2(ry, rx)

        K_eff = K0 + kappa * alarm(r)

        dtheta = omega + (K_eff * r)[:, None] * np.sin(psi[:, None] - theta)
        noise = sigma_noise * np.sqrt(dt) * rng.standard_normal((n_trials, N))
        theta_new = theta + dt * dtheta + noise
        theta = np.where(alive[:, None], theta_new, theta)

        newly_tipped = alive & (r > r_tip)
        tipped_at = np.where(newly_tipped, t, tipped_at)
        tipped |= newly_tipped
        alive &= ~newly_tipped
        if not alive.any():
            break

    return tipped, tipped_at


if __name__ == "__main__":
    t0 = time.time()
    n_trials = 500
    SEED0 = 20_000  # jamais utilisee par le script fourni (seed=42)
    print(f"=== Contre-verification independante Kuramoto ({n_trials} graines, seed0={SEED0}) ===")
    print(f"{'kappa':>10} {'bascules':>12} {'taux':>8} {'t_moy':>10}")
    for kappa in [0.0, -0.5, -1.0, -2.0, 0.5, 1.0, 2.0]:
        tipped, tipped_at = simulate_batch(kappa, n_trials, seed0=SEED0 + int(round(kappa * 100)))
        n_tip = int(tipped.sum())
        rate_tip = n_tip / n_trials
        valid = tipped_at[tipped]
        t_mean = np.nanmean(valid) if len(valid) else float("nan")
        label = "0 (classique)" if kappa == 0 else f"{kappa:+.2f}"
        print(f"{label:>10} {n_tip:5d}/{n_trials:<5d} {rate_tip:7.1%} {t_mean:10.1f}")
    print(f"\nTemps ecoule: {time.time()-t0:.1f}s")

    print()
    print("=== Sensibilite au seuil r_tip (l'asymetrie documentee dans le §0bis) ===")
    for r_tip in [0.6, 0.7, 0.8]:
        print(f"--- r_tip={r_tip} ---")
        for kappa in [-2.0, -1.0, 0.0, 1.0, 2.0]:
            tipped, tipped_at = simulate_batch(kappa, 200, seed0=SEED0 + 500, r_tip=r_tip)
            n_tip = int(tipped.sum())
            print(f"    kappa={kappa:+.1f} : {n_tip}/200 ({n_tip/200:.1%})")
