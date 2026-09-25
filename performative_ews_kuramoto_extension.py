"""
Extension du modele performatif au cas Kuramoto (H4) :
K_eff(t) = K0(t) + kappa * A(r_t)
kappa < 0 : la coalition se desynchronise quand l'alerte monte (auto-invalidant)
kappa > 0 : la coalition se synchronise davantage quand l'alerte monte (panique/auto-realisateur)
"""
import numpy as np
import time

def alarm(r, k=10.0, rc=0.5):
    return 1.0/(1.0+np.exp(-k*(r-rc)))

def simulate_kuramoto_batch(kappa, n_trials, N=60, K0_start=0.5, K0_rate=0.003,
                             dt=0.05, T=500, seed=0, r_tip=0.7, sigma_omega=0.5):
    rng = np.random.default_rng(seed)
    steps = int(T/dt)

    omega = rng.normal(0, sigma_omega, size=(n_trials, N))
    theta = rng.uniform(0, 2*np.pi, size=(n_trials, N))
    tipped = np.zeros(n_trials, dtype=bool)
    tipped_at = np.full(n_trials, np.nan)
    alive = np.ones(n_trials, dtype=bool)
    max_r = np.zeros(n_trials)

    for step in range(steps):
        t = step*dt
        K0 = K0_start + K0_rate*t

        cos_t, sin_t = np.cos(theta), np.sin(theta)
        rx, ry = cos_t.mean(axis=1), sin_t.mean(axis=1)
        r = np.sqrt(rx**2 + ry**2)
        psi = np.arctan2(ry, rx)
        max_r = np.maximum(max_r, r)

        a = alarm(r, k=10.0, rc=0.5)
        K_eff = K0 + kappa*a

        dtheta = omega + K_eff[:, None]*r[:, None]*np.sin(psi[:, None] - theta)
        theta_new = theta + dt*dtheta + 0.05*np.sqrt(dt)*rng.standard_normal((n_trials, N))
        theta = np.where(alive[:, None], theta_new, theta)

        newly_tipped = alive & (r > r_tip)
        tipped_at = np.where(newly_tipped, t, tipped_at)
        tipped |= newly_tipped
        alive &= ~newly_tipped
        if not alive.any():
            break

    return tipped, tipped_at, max_r

if __name__ == "__main__":
    t0 = time.time()
    n_trials = 500
    print(f"=== Extension Kuramoto : {n_trials} graines par kappa ===")
    print(f"{'kappa':>10} {'bascules':>12} {'taux':>8} {'t_moy':>10}")
    for kappa in [0.0, -0.5, -1.0, -2.0, 0.5, 1.0, 2.0]:
        tipped, tipped_at, max_r = simulate_kuramoto_batch(kappa, n_trials, seed=42)
        n_tip = tipped.sum()
        rate_tip = n_tip/n_trials
        valid_times = tipped_at[tipped]
        t_mean = np.nanmean(valid_times) if len(valid_times) else float('nan')
        label = "0 (classique)" if kappa==0 else f"{kappa:+.2f}"
        print(f"{label:>10} {n_tip:5d}/{n_trials:<5d} {rate_tip:7.1%} {t_mean:10.1f}")
    print(f"\nTemps ecoule: {time.time()-t0:.1f}s")
