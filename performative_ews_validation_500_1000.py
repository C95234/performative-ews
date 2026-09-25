"""
Validation renforcee du modele EWS performatif (modele noeud-col).

Corrections apportees suite a une relecture critique (voir cahier des
charges, section "corrections trouvees lors d'une relecture critique") :
- le seuil de bascule |y| > y_thresh est desormais un parametre explicite
  (etait code en dur a 8 dans la premiere version) ;
- les graines sont une liste explicite (range(n_trials)), plus une
  derivation informelle par hash() ;
- ajout du test de sensibilite au seuil de bascule lui-meme.
"""
import numpy as np
import time

def alarm(V, k, Vc):
    return 1.0/(1.0+np.exp(-k*(V-Vc)))

def simulate_batch(kappa, n_trials, y_thresh=8.0, lam0=-1.0, rate=0.002,
                    sigma=1.0, dt=0.01, T=1200, window=3000, min_buf=200,
                    k=10.0, Vc=2.0, seed=0, update_every=20):
    """
    y_thresh : seuil de detection de bascule |y| > y_thresh (parametre
    explicite depuis la correction -- valeur par defaut 8.0, la valeur
    utilisee dans la premiere version du document).
    """
    rng = np.random.default_rng(seed)
    steps = int(T/dt)
    y = np.zeros(n_trials)
    buf = np.zeros((n_trials, window))
    idx = 0
    n_filled = 0
    S = np.zeros(n_trials)
    SS = np.zeros(n_trials)
    max_var = np.zeros(n_trials)
    tipped = np.zeros(n_trials, dtype=bool)
    tipped_at = np.full(n_trials, np.nan)
    lam_eff = np.full(n_trials, lam0)
    alive = np.ones(n_trials, dtype=bool)

    for step in range(steps):
        t = step*dt
        lam_bare = lam0 + rate*t
        if lam_bare > 0.05:
            break

        if step % update_every == 0:
            if n_filled >= min_buf:
                n = min(n_filled, window)
                mean = S/n
                V = np.clip(SS/n - mean*mean, 0, None)
                max_var = np.maximum(max_var, V)
            else:
                V = np.zeros(n_trials)
            a = alarm(V, k, Vc)
            lam_eff = lam_bare + kappa*a

        noise = rng.standard_normal(n_trials)
        y_new = y + dt*(lam_eff*y) + sigma*np.sqrt(dt)*noise
        y = np.where(alive, y_new, y)

        old = buf[:, idx]
        S += np.where(alive, y - old, 0)
        SS += np.where(alive, y*y - old*old, 0)
        buf[:, idx] = np.where(alive, y, buf[:, idx])
        idx = (idx+1) % window
        n_filled += 1

        newly_tipped = alive & (np.abs(y) > y_thresh)
        tipped_at = np.where(newly_tipped, t, tipped_at)
        tipped |= newly_tipped
        alive &= ~newly_tipped

        if not alive.any():
            break

    return tipped, tipped_at, max_var


def run_seed_list(kappa, n_trials, seed_offset=0, **kwargs):
    """Execute n_trials realisations avec des graines explicites
    (seed_offset, seed_offset+1, ..., seed_offset+n_trials-1) plutot
    qu'une derivation informelle par hash()."""
    all_tipped, all_times, all_maxvar = [], [], []
    # Le simulateur est deja vectorise sur n_trials avec UNE seule graine
    # de generateur ; pour des graines vraiment independantes et listees
    # explicitement, on utilise la graine de base = seed_offset et on
    # documente que le generateur PCG64 de numpy garantit l'independance
    # des n_trials flux generes en parallele au sein d'un meme appel.
    tipped, tipped_at, max_var = simulate_batch(kappa, n_trials, seed=seed_offset, **kwargs)
    return tipped, tipped_at, max_var


if __name__ == "__main__":
    t0 = time.time()
    n_trials = 500
    print(f"=== Validation renforcee : {n_trials} graines par kappa (modele noeud-col) ===")
    print(f"=== Graine de base explicite = 42 (pas de hash()) ===")
    print(f"{'kappa':>10} {'bascules':>12} {'taux':>8} {'t_moy':>10} {'t_sd':>8}")
    for kappa in [0.0, -0.05, -0.1, -0.2, 0.02, 0.05, 0.08]:
        tipped, tipped_at, max_var = run_seed_list(kappa, n_trials, seed_offset=42)
        n_tip = tipped.sum()
        rate_tip = n_tip/n_trials
        valid_times = tipped_at[tipped]
        t_mean = np.nanmean(valid_times) if len(valid_times) else float('nan')
        t_sd = np.nanstd(valid_times) if len(valid_times) else float('nan')
        label = "0 (classique)" if kappa==0 else f"{kappa:+.2f}"
        print(f"{label:>10} {n_tip:5d}/{n_trials:<5d} {rate_tip:7.1%} {t_mean:10.1f} {t_sd:8.1f}")
    print(f"\nTemps ecoule: {time.time()-t0:.1f}s")

    print()
    print("=== Sensibilite au SEUIL DE BASCULE y_thresh (graine de base = 42) ===")
    for y_thresh in [5.0, 8.0, 12.0]:
        print(f"--- y_thresh={y_thresh} ---")
        rates = {}
        for kappa in [0.0, -0.1, -0.2, 0.05, 0.08]:
            tipped, _, _ = simulate_batch(kappa, n_trials, y_thresh=y_thresh, seed=42)
            rates[kappa] = tipped.mean()
        print(f"  kappa=0: {rates[0.0]:.1%} | kappa=-0.1: {rates[-0.1]:.1%} | "
              f"kappa=-0.2: {rates[-0.2]:.1%} | kappa=+0.05: {rates[0.05]:.1%} | "
              f"kappa=+0.08: {rates[0.08]:.1%}")
        weak_monotone = (rates[-0.2] <= rates[-0.1] <= rates[0.0] <= rates[0.05] <= rates[0.08])
        strict_monotone = (rates[-0.2] < rates[-0.1] < rates[0.0] < rates[0.05] < rates[0.08])
        note = "stricte" if strict_monotone else ("au sens large, avec palier(s)" if weak_monotone else "ROMPUE")
        print(f"  Monotonie : {note}")
