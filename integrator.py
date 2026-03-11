"""
integrator.py — Intégration numérique (PyTorch)

Ce module est responsable uniquement des schémas d'intégration :
  - Intégration scalaire (quadrature) pour valider des intégrales spectrales
  - Intégration de l'équation de transfert radiatif (RTE)

Il est complètement indépendant de la physique : il reçoit des tenseurs
ou des tableaux en entrée et retourne des valeurs numériques.

Équation de transfert radiatif (1D, plan-parallèle, propagation verticale) :
  dI/dz = −alpha(z) * [I(z) − B(nu, T(z))]

Schéma d'intégration par couche homogène :
  tau   = alpha(z_mid) * dz_cm         [sans dimension]
  I_out = I_in * exp(−tau) + B(T_mid) * (1 − exp(−tau))

Ce schéma est exact pour une couche d'absorption constante. Il est
inconditionnellement stable et garantit physiquement 0 ≤ I_out ≤ max(I_in, B).

Toutes les opérations spectrales (sur l'axe nu) sont vectorisées avec torch.
La seule boucle Python restante est la boucle sur les couches verticales z,
ce qui permet de garder un empreinte mémoire bornée (pas de tensor 3D nu×z).
"""

import torch


class Integrator:
    """
    Intégrateurs scalaires (quadrature 1D sur tableaux torch).

    Toutes les méthodes opèrent sur des tenseurs X, Y représentant une fonction
    échantillonnée. Les calculs se font sur le device des tenseurs d'entrée.
    """

    def integrateSquare(self, X: torch.Tensor, Y: torch.Tensor) -> float:
        """Méthode des rectangles (ordre 1). Borne gauche."""
        dX = torch.diff(X)
        return float(torch.sum(dX * Y[:-1]).item())

    def integrateTrapezoidal(self, X: torch.Tensor, Y: torch.Tensor) -> float:
        """Méthode des trapèzes (ordre 2)."""
        dX  = torch.diff(X)
        midY = 0.5 * (Y[:-1] + Y[1:])
        return float(torch.sum(dX * midY).item())

    def integrateModifiedSimpson(self, X: torch.Tensor, Y: torch.Tensor) -> float:
        """Règle de Simpson par intervalles avec point milieu interpolé (ordre 2)."""
        dX   = torch.diff(X)
        midY = 0.5 * (Y[:-1] + Y[1:])
        return float(torch.sum(dX / 6.0 * (Y[:-1] + 4.0 * midY + Y[1:])).item())

    def integrateMonteCarlo(self, X: torch.Tensor, Y: torch.Tensor,
                             niter: int = 10_000) -> float:
        """
        Intégration Monte-Carlo par interpolation linéaire sur [X[0], X[-1]].

        Les tirages aléatoires sont effectués sur le même device que X.
        """
        x_rand = torch.empty(niter, dtype=X.dtype, device=X.device).uniform_(
            float(X[0]), float(X[-1])
        )
        # torch.searchsorted pour une interpolation linéaire vectorisée
        idx  = torch.searchsorted(X.contiguous(), x_rand.contiguous()).clamp(1, len(X) - 1)
        x0   = X[idx - 1];  x1 = X[idx]
        y0   = Y[idx - 1];  y1 = Y[idx]
        t    = (x_rand - x0) / (x1 - x0)
        y_rand = y0 + t * (y1 - y0)
        return float(((X[-1] - X[0]) * y_rand.mean()).item())


class RTEIntegrator:
    """
    Intégrateur de l'équation de transfert radiatif (RTE) — PyTorch.

    Propage un spectre I(nu) depuis z=0 jusqu'à z_max en appliquant
    le schéma exact couche homogène à chaque pas dz.

    Seule la boucle sur les couches verticales reste en Python.
    Toutes les opérations spectrales (sur l'axe nu) sont entièrement
    vectorisées via torch — aucun transfert CPU↔GPU.

    Paramètres
    ----------
    solver : Solver
        Fournit solver.alpha(nu, z) [cm-1] et solver.planck(nu, T) [W/m²/sr/cm-1].
        Le device et le dtype sont hérités de solver.device / solver.dtype.
    """

    def __init__(self, solver):
        self.solver = solver

    def propagate(self, nu: torch.Tensor, dz_m: float, z_max: float,
                  verbose: bool = True):
        """
        Propage le rayonnement de z=0 à z_max.

        Paramètres
        ----------
        nu      : torch.Tensor [cm-1]   grille spectrale — sur solver.device
        dz_m    : float [m]             pas vertical
        z_max   : float [m]             altitude maximale
        verbose : bool                  affiche la progression

        Retour
        ------
        I_surface : torch.Tensor [W/m²/sr/cm-1]  spectre corps noir à z=0
        I_top     : torch.Tensor [W/m²/sr/cm-1]  spectre au sommet z_max
        z_steps   : torch.Tensor [m]              altitudes des couches
        """
        dz_cm   = dz_m * 100.0
        n_steps = int(z_max / dz_m)

        # Grille z sur le device — pas de numpy
        z_steps = torch.arange(n_steps,
                               dtype=self.solver.dtype,
                               device=self.solver.device) * dz_m

        # Condition initiale : corps noir à la surface
        T_surf    = self.solver.surface_temperature()
        I_surface = self.solver.planck(nu, T_surf)
        I         = I_surface.clone()

        progress_interval = max(1, n_steps // 20)

        # -------------------------------------------------------
        # Boucle verticale — seule boucle Python autorisée
        # Toutes les opérations internes sont sur torch (nu vectorisé)
        # -------------------------------------------------------
        for step in range(n_steps):
            z_mid = (step + 0.5) * dz_m

            # Propriétés scalaires depuis le modèle atmo (Python pur, OK)
            T_mid, _, _, _ = self.solver.atm.getAtmProperties(z_mid)

            # alpha : torch.Tensor (N_nu,) — sur device
            alpha = self.solver.alpha(nu, z_mid)

            # Épaisseur optique (N_nu,)
            tau   = alpha * dz_cm

            # Luminance d'équilibre local (N_nu,)
            B     = self.solver.planck(nu, T_mid)

            # Schéma exact couche homogène — tout vectorisé
            e_tau = torch.exp(-tau)
            I     = I * e_tau + B * (1.0 - e_tau)

            if verbose and step % progress_interval == 0:
                pct = 100.0 * step / n_steps
                print(f"    z = {step * dz_m:7.1f} m  ({pct:5.1f}%)")

        if verbose:
            print(f"    z = {z_max:7.1f} m  (100.0%)")

        return I_surface, I, z_steps