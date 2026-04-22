"""
controlador.py - Diseño del controlador (Parte B, Ejercicio 1)

Contiene funciones para:
- Diseño de compensadores lead/lag por root locus
- Evaluación de especificaciones de desempeño
- Comparación Bode compensado vs no compensado
- Seguimiento de rampa y adición de integrador

Curso: IOR442 Sistemas de Control - Dr. Mariano Scaramal
"""

import control as ctrl
import numpy as np
import matplotlib.pyplot as plt
from planta import crear_planta, crear_planta_perturbacion, Kt, J, B, L, R, Ke


# =============================================================================
# Especificaciones de diseño (Parte B.a)
# =============================================================================
SPECS = {
    'ts_max': 7.5,       # Tiempo de establecimiento máximo [s] (2%)
    'Mp_max': 15.0,      # Sobrepico máximo [%]
    'error_pert_max': 5.0,  # Error máximo ante perturbación 20% de ref [%]
    'amp_pert': 0.2,     # Amplitud perturbación (20% de escalón unitario)
}

# De las especificaciones derivamos restricciones en el plano s:
#   Mp ≤ 15% → ζ ≥ 0.5169
#   ts ≤ 7.5s (2%) → σ = 4/ts ≥ 0.5333 rad/s
ZETA_MIN = -np.log(SPECS['Mp_max'] / 100) / np.sqrt(np.pi**2 + np.log(SPECS['Mp_max'] / 100)**2)
SIGMA_MIN = 4.0 / SPECS['ts_max']


def graficar_root_locus_con_region(G, C=None, titulo="Root Locus"):
    """
    Grafica el root locus del sistema con la región admisible sombreada.

    La región admisible es la intersección de:
    - Re(s) ≤ -σ_min  (condición de ts)
    - Dentro del cono de amortiguamiento ζ ≥ ζ_min

    Args:
        G: planta
        C: controlador (si se quiere ver el RL del sistema compensado)
        titulo: título del gráfico
    """
    if C is not None:
        L_sys = C * G
    else:
        L_sys = G

    fig, ax = plt.subplots(figsize=(10, 8))

    # Dibujar root locus
    rlist, klist = ctrl.root_locus(L_sys, plot=False)

    for i in range(rlist.shape[1]):
        ax.plot(np.real(rlist[:, i]), np.imag(rlist[:, i]), 'b-', linewidth=1)

    # Marcar polos (x) y ceros (o) de lazo abierto
    polos = ctrl.poles(L_sys)
    ceros = ctrl.zeros(L_sys)
    ax.plot(np.real(polos), np.imag(polos), 'rx', markersize=10, markeredgewidth=2, label='Polos LA')
    if len(ceros) > 0:
        ax.plot(np.real(ceros), np.imag(ceros), 'go', markersize=10, markeredgewidth=2, label='Ceros LA')

    # Región de σ_min (línea vertical)
    y_lim = max(abs(np.imag(rlist).max()), 15)
    ax.axvline(-SIGMA_MIN, color='r', linestyle='--', alpha=0.5, label=f'σ_min = {SIGMA_MIN:.3f}')

    # Región de ζ_min (líneas desde el origen)
    theta = np.arccos(ZETA_MIN)
    r_max = y_lim * 1.5
    ax.plot([0, -r_max * np.cos(theta)], [0, r_max * np.sin(theta)],
            'm--', alpha=0.5, label=f'ζ_min = {ZETA_MIN:.3f}')
    ax.plot([0, -r_max * np.cos(theta)], [0, -r_max * np.sin(theta)],
            'm--', alpha=0.5)

    # Sombrear región admisible
    from matplotlib.patches import Polygon
    # Construir la región como un polígono (simplificada)
    sigma_lim = -SIGMA_MIN
    y_top = abs(sigma_lim) * np.tan(theta)

    region_x = [sigma_lim, -r_max, -r_max, sigma_lim]
    region_y_top = [y_top, r_max * np.sin(theta), 0, 0]
    region_y_bot = [-y_top, -r_max * np.sin(theta), 0, 0]

    # Región simplificada: rectángulo + conos
    ax.fill_betweenx(
        np.linspace(-y_lim, y_lim, 100),
        -r_max, sigma_lim,
        alpha=0.05, color='green'
    )

    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)
    ax.set_xlabel('Re(s)')
    ax.set_ylabel('Im(s)')
    ax.set_title(titulo)
    ax.legend(loc='best')
    ax.set_xlim([-20, 5])
    ax.set_ylim([-y_lim, y_lim])
    plt.tight_layout()
    return fig


def disenar_lead_gridsearch(G, verbose=True):
    """
    Diseña un compensador lead genuino (|polo| > |cero|) mediante grid search
    sobre pares (z_c, s_d):
      - z_c  : posición del cero  (cero en s = -z_c, z_c > 0)
      - s_d  : polo dominante deseado en la región admisible (σ, ζ)

    Para cada par, calcula p_c y K por contribución angular + condición de
    magnitud del root locus, acepta solo soluciones con p_c > z_c (lead).
    Devuelve el mejor compensador que cumple ts y Mp (y minimiza e_pert).

    Returns:
        mejor : dict con claves z_c, p_c, K, C, ts, Mp, e_pert, C0
        todos : lista completa de soluciones factibles
    """
    from analisis import calcular_metricas_temporales

    mejores = []

    z_vals     = np.arange(0.2, 8.0, 0.4)
    sigma_vals = np.arange(SIGMA_MIN + 0.2, 9.0, 0.4)
    zeta_vals  = np.arange(ZETA_MIN + 0.02, 0.92, 0.08)

    for z_c in z_vals:
        for sigma_d in sigma_vals:
            for zeta_d in zeta_vals:
                omega_d = sigma_d * np.tan(np.arccos(zeta_d))
                sd = complex(-sigma_d, omega_d)

                # --- Contribución angular de la planta ---
                G_val    = ctrl.evalfr(G, sd)
                angle_G  = np.degrees(np.angle(G_val))

                # Fase que debe aportar el compensador para root locus (lead: > 0)
                phi_c = -180.0 - angle_G
                # Normalizar al primer positivo en (-360°, 360°)
                phi_c = phi_c % 360
                if phi_c > 180:
                    phi_c -= 360
                if phi_c <= 0:
                    continue          # necesita fase negativa → lag, descartar

                # --- Contribución del cero ---
                phi_z = np.degrees(np.angle(sd + z_c))
                if not (0 < phi_z < 180):
                    continue

                # --- Ángulo requerido del polo ---
                phi_p = phi_z - phi_c
                if not (0 < phi_p < 180):
                    continue

                # --- Posición del polo ---
                p_c = omega_d / np.tan(np.radians(phi_p)) + sigma_d
                if p_c <= z_c + 0.05:   # condición lead: |polo| > |cero|
                    continue

                # --- Ganancia por condición de magnitud ---
                C_base  = ctrl.tf([1, z_c], [1, p_c])
                CG_val  = ctrl.evalfr(C_base * G, sd)
                if abs(CG_val) < 1e-12:
                    continue
                K = 1.0 / abs(CG_val)
                C = K * C_base

                # --- Estabilidad del lazo cerrado ---
                T      = ctrl.feedback(C * G, 1)
                p_lc   = ctrl.poles(T)
                if any(np.real(p_lc) > 1e-6):
                    continue

                # --- Evaluación de specs ---
                try:
                    met = calcular_metricas_temporales(G, C, t_final=40)
                    ts_v = met['ts']
                    Mp_v = met['Mp']
                    if ts_v > SPECS['ts_max'] or Mp_v > SPECS['Mp_max']:
                        continue

                    C0      = float(abs(ctrl.evalfr(C, 1e-10)))
                    e_pert  = 64.0 / C0 if C0 > 0 else np.inf

                    mejores.append({
                        'z_c': round(z_c, 2), 'p_c': round(p_c, 4),
                        'K':   round(K, 4),   'C':   C,
                        'C0':  round(C0, 4),  'ts':  round(ts_v, 3),
                        'Mp':  round(Mp_v, 3),'e_pert': round(e_pert, 2),
                        'sd':  sd, 'sigma_d': sigma_d, 'zeta_d': zeta_d,
                    })
                except Exception:
                    continue

    if not mejores:
        raise RuntimeError("Grid search no encontró ningún compensador lead factible.")

    # Prioridad: menor e_pert, desempate por ts
    mejores.sort(key=lambda x: (x['e_pert'], x['ts']))
    best = mejores[0]

    if verbose:
        sep = '=' * 56
        print(f"\n{sep}")
        print(f"  GRID SEARCH — Compensador Lead (p_c > z_c)")
        print(sep)
        print(f"  Cero : s = -{best['z_c']:.2f}")
        print(f"  Polo : s = -{best['p_c']:.4f}   (ratio p/z = {best['p_c']/best['z_c']:.2f})")
        print(f"  K    = {best['K']:.4f}")
        print(f"  C(s) = {best['K']:.4f} · (s + {best['z_c']}) / (s + {best['p_c']:.4f})")
        print(f"  C(0) = {best['C0']:.4f}")
        print(f"  ts   = {best['ts']:.2f} s    (máx {SPECS['ts_max']} s)  ✓")
        print(f"  Mp   = {best['Mp']:.2f} %   (máx {SPECS['Mp_max']} %)  ✓")
        print(f"  e_pert ≈ {best['e_pert']:.1f} %  (máx {SPECS['error_pert_max']} %)  ✗")
        print(sep)

    return best, mejores


def disenar_lead_lag_gridsearch(G, C0_min=13.5, verbose=True):
    """
    Diseña un compensador lead-lag para cumplir simultáneamente ts<=7.5, Mp<=15%, e_pert<5%.

    La estructura es:
        C(s) = Kc * (s+z_L)/(s+p_L) * (s+z_lag)/(s+p_lag)
      - Lead (p_L > z_L): aporta fase en banda media → controla ts y Mp.
      - Lag  (z_lag > p_lag, con z_lag, p_lag << wcp): sube ganancia DC sin
        degradar fase en cruce → permite C(0) > 12.8 y así e_pert < 5%.

    Kc se calcula por condición de magnitud del root locus en sd (el polo dominante
    deseado), lo que garantiza que el RL pase exactamente por sd.
    e_pert se calcula con la fórmula analítica exacta: e_pert[%] = 64 / C(0).

    Returns:
        best : dict con z_L, p_L, z_lag, p_lag, Kc, C, C0, ts, Mp, e_pert
        todos: lista de todos los factibles
    """
    from analisis import calcular_metricas_temporales

    mejores = []

    z_vals      = np.arange(0.2, 8.0, 0.4)
    sigma_vals  = np.arange(SIGMA_MIN + 0.2, 9.0, 0.4)
    zeta_vals   = np.arange(ZETA_MIN + 0.02, 0.92, 0.08)
    p_lag_vals  = [0.005, 0.01, 0.02]
    ratio_lags  = [15, 20, 25, 30]

    for z_L in z_vals:
        for sigma_d in sigma_vals:
            for zeta_d in zeta_vals:
                omega_d = sigma_d * np.tan(np.arccos(zeta_d))
                sd = complex(-sigma_d, omega_d)

                # --- Contribución angular lead ---
                G_val   = ctrl.evalfr(G, sd)
                angle_G = np.degrees(np.angle(G_val))

                phi_c = (-180.0 - angle_G) % 360
                if phi_c > 180:
                    phi_c -= 360
                if phi_c <= 0:
                    continue   # necesita fase negativa → solo lead no sirve aquí

                phi_z = np.degrees(np.angle(sd + z_L))
                if not (0 < phi_z < 180):
                    continue

                phi_p = phi_z - phi_c
                if not (0 < phi_p < 180):
                    continue

                p_L = omega_d / np.tan(np.radians(phi_p)) + sigma_d
                if p_L <= z_L + 0.05:
                    continue

                C_lead_base = ctrl.tf([1, z_L], [1, p_L])

                # Estimar wcp con lead solo (para filtro z_lag << wcp)
                try:
                    K_lead_est = 1.0 / abs(ctrl.evalfr(C_lead_base * G, sd))
                    _, _, _, wcp_est = ctrl.margin(K_lead_est * C_lead_base * G)
                    if not (np.isfinite(wcp_est) and wcp_est > 0):
                        wcp_est = 1.0
                except Exception:
                    wcp_est = 1.0

                # --- Agregar lag y verificar ---
                for p_lag in p_lag_vals:
                    for ratio_lag in ratio_lags:
                        z_lag = p_lag * ratio_lag

                        # Condición: z_lag << wcp (máximo wcp/5)
                        if z_lag > wcp_est / 5:
                            continue

                        C_ll_base = C_lead_base * ctrl.tf([1, z_lag], [1, p_lag])

                        CllG_val = ctrl.evalfr(C_ll_base * G, sd)
                        if abs(CllG_val) < 1e-12:
                            continue
                        Kc = 1.0 / abs(CllG_val)
                        C_full = Kc * C_ll_base

                        # Ganancia DC y e_pert analíticos (exactos por FVT)
                        C0 = float(abs(ctrl.evalfr(C_full, 1e-10)))
                        if C0 < C0_min:
                            continue
                        e_pert = 64.0 / C0
                        if e_pert >= SPECS['error_pert_max']:
                            continue

                        # Estabilidad
                        try:
                            T = ctrl.feedback(C_full * G, 1)
                            p_lc = ctrl.poles(T)
                            if any(np.real(p_lc) > 1e-6):
                                continue
                        except Exception:
                            continue

                        # Simulación transitoria (solo candidatos que pasan lo anterior)
                        try:
                            met = calcular_metricas_temporales(G, C_full, t_final=40)
                            ts_v = met['ts']
                            Mp_v = met['Mp']
                            if ts_v > SPECS['ts_max'] or Mp_v > SPECS['Mp_max']:
                                continue
                        except Exception:
                            continue

                        mejores.append({
                            'z_L':   round(z_L,    4),
                            'p_L':   round(p_L,    4),
                            'z_lag': round(z_lag,  6),
                            'p_lag': round(p_lag,  6),
                            'Kc':    round(Kc,     4),
                            'C':     C_full,
                            'C0':    round(C0,     4),
                            'ts':    round(ts_v,   3),
                            'Mp':    round(Mp_v,   3),
                            'e_pert':round(e_pert, 2),
                            'sd':    sd,
                            'sigma_d': sigma_d,
                            'zeta_d':  zeta_d,
                        })

    if not mejores:
        raise RuntimeError(
            "Grid search lead-lag no encontró ningún compensador factible. "
            "Intentá reducir C0_min o ampliar la grilla."
        )

    # Prioridad: menor Mp (transitorio más suave), desempate por e_pert
    mejores.sort(key=lambda x: (x['Mp'], x['e_pert']))
    best = mejores[0]

    if verbose:
        sep = '=' * 62
        print(f"\n{sep}")
        print(f"  GRID SEARCH — Compensador Lead-Lag")
        print(sep)
        print(f"  Lead: cero = -{best['z_L']:.4f},  polo = -{best['p_L']:.4f}")
        print(f"  Lag:  cero = -{best['z_lag']:.5f}, polo = -{best['p_lag']:.5f}")
        print(f"        ratio lag = {best['z_lag']/best['p_lag']:.1f}")
        print(f"  Kc   = {best['Kc']:.4f}")
        print(f"  C(0) = {best['C0']:.4f}  (mín requerido: {C0_min})")
        print(f"  ts   = {best['ts']:.2f} s    (máx {SPECS['ts_max']} s)  ✓")
        print(f"  Mp   = {best['Mp']:.2f} %   (máx {SPECS['Mp_max']} %)  ✓")
        print(f"  e_pert ≈ {best['e_pert']:.1f} %  (máx {SPECS['error_pert_max']} %)  ✓")
        print(sep)

    return best, mejores


def disenar_pi_gridsearch(G, verbose=True):
    """
    Diseña un compensador PI puro C(s) = Kc * (s + z_i) / s para cumplir
    SIMULTÁNEAMENTE ts<=7.5, Mp<=15, e_pert<5%.

    Justificación:
      - El sistema sin compensar ya cumple ts y Mp (G es Tipo 1 con PM=72°,
        Mp≈0), pero falla e_pert (= 64%).
      - Agregar un integrador eleva el tipo del sistema a 2 => e_pert = 0.
      - El cero z_i se ubica << wcp para que el aporte de fase negativa del
        integrador quede casi totalmente cancelado por el cero en la banda
        de cruce, preservando ts y Mp.

    Kc se recalcula por condición de magnitud del root locus en el polo
    dominante deseado sd.

    Returns:
        best : dict con z_i, Kc, C, ts, Mp, sd
        todos: lista de todos los factibles
    """
    from analisis import calcular_metricas_temporales

    mejores = []
    z_i_vals   = np.concatenate([np.arange(0.02, 0.5, 0.02),
                                  np.arange(0.5, 3.0, 0.1)])
    sigma_vals = np.arange(SIGMA_MIN + 0.2, 9.0, 0.3)
    zeta_vals  = np.arange(ZETA_MIN + 0.02, 0.98, 0.05)

    for sigma_d in sigma_vals:
        for zeta_d in zeta_vals:
            omega_d = sigma_d * np.tan(np.arccos(zeta_d))
            sd = complex(-sigma_d, omega_d)
            for z_i in z_i_vals:
                C_base = ctrl.tf([1, z_i], [1, 0])
                CG = ctrl.evalfr(C_base * G, sd)
                ang = np.degrees(np.angle(CG)) % 360
                if ang > 180:
                    ang -= 360
                if abs(abs(ang) - 180) > 3:
                    continue
                Kc = 1.0 / abs(CG)
                C_full = Kc * C_base
                try:
                    T = ctrl.feedback(C_full * G, 1)
                    if any(np.real(ctrl.poles(T)) > 1e-6):
                        continue
                    met = calcular_metricas_temporales(G, C_full, t_final=60)
                    if met['ts'] > SPECS['ts_max'] or met['Mp'] > SPECS['Mp_max']:
                        continue
                    mejores.append({
                        'z_i': round(z_i, 4),
                        'Kc':  round(Kc, 4),
                        'C':   C_full,
                        'ts':  round(met['ts'], 3),
                        'Mp':  round(met['Mp'], 3),
                        'sd':  sd, 'sigma_d': sigma_d, 'zeta_d': zeta_d,
                    })
                except Exception:
                    continue

    if not mejores:
        raise RuntimeError("Grid search PI no encontró compensador factible.")

    mejores.sort(key=lambda x: (x['Mp'], x['ts']))
    best = mejores[0]

    if verbose:
        sep = '=' * 60
        print(f"\n{sep}")
        print(f"  GRID SEARCH — Compensador PI (Tipo 2)")
        print(sep)
        print(f"  Cero: s = -{best['z_i']:.4f}   Polo: s = 0 (integrador)")
        print(f"  Kc   = {best['Kc']:.4f}")
        print(f"  C(s) = {best['Kc']:.4f} · (s + {best['z_i']}) / s")
        print(f"  ts   = {best['ts']:.2f} s    (máx {SPECS['ts_max']} s)  ✓")
        print(f"  Mp   = {best['Mp']:.2f} %   (máx {SPECS['Mp_max']} %)  ✓")
        print(f"  e_pert = 0   (Tipo 2)  ✓")
        print(f"  e_rampa = 0  (Tipo 2)  ✓")
        print(sep)

    return best, mejores


def disenar_lead_pi_gridsearch(G, verbose=True):
    """
    Diseña un compensador lead + PI para cumplir SIMULTÁNEAMENTE:
      ts <= 7.5 s, Mp <= 15%, e_pert < 5%, e_rampa = 0.

    Estructura:
        C(s) = Kc * (s+z_L)/(s+p_L) * (s+z_i)/s
      - Lead (p_L > z_L): aporta fase en banda media → controla ts y Mp.
      - PI   ((s+z_i)/s, con z_i << wcp): eleva el tipo del sistema a 2
        (=> e_pert = 0, e_rampa = 0) sin degradar la fase en cruce si z_i
        se ubica bien por debajo de wcp.

    Kc se recalcula por condición de magnitud en sd del compensador combinado.

    Returns:
        best : dict con z_L, p_L, z_i, Kc, C, ts, Mp, e_pert, e_rampa
    """
    from analisis import calcular_metricas_temporales

    mejores = []

    z_vals      = np.arange(0.2, 8.0, 0.4)
    sigma_vals  = np.arange(SIGMA_MIN + 0.2, 9.0, 0.4)
    zeta_vals   = np.arange(ZETA_MIN + 0.02, 0.92, 0.08)
    z_i_vals    = [0.02, 0.05, 0.1, 0.2, 0.3, 0.5]

    for z_L in z_vals:
        for sigma_d in sigma_vals:
            for zeta_d in zeta_vals:
                omega_d = sigma_d * np.tan(np.arccos(zeta_d))
                sd = complex(-sigma_d, omega_d)

                G_val   = ctrl.evalfr(G, sd)
                angle_G = np.degrees(np.angle(G_val))
                phi_c = (-180.0 - angle_G) % 360
                if phi_c > 180:
                    phi_c -= 360
                if phi_c <= 0:
                    continue

                phi_z = np.degrees(np.angle(sd + z_L))
                if not (0 < phi_z < 180):
                    continue

                phi_p = phi_z - phi_c
                if not (0 < phi_p < 180):
                    continue

                p_L = omega_d / np.tan(np.radians(phi_p)) + sigma_d
                if p_L <= z_L + 0.05:
                    continue

                C_lead_base = ctrl.tf([1, z_L], [1, p_L])

                try:
                    K_lead_est = 1.0 / abs(ctrl.evalfr(C_lead_base * G, sd))
                    _, _, _, wcp_est = ctrl.margin(K_lead_est * C_lead_base * G)
                    if not (np.isfinite(wcp_est) and wcp_est > 0):
                        wcp_est = 1.0
                except Exception:
                    wcp_est = 1.0

                for z_i in z_i_vals:
                    if z_i > wcp_est / 5:
                        continue

                    PI = ctrl.tf([1, z_i], [1, 0])
                    C_base = C_lead_base * PI

                    CG_val = ctrl.evalfr(C_base * G, sd)
                    if abs(CG_val) < 1e-12:
                        continue
                    Kc = 1.0 / abs(CG_val)
                    C_full = Kc * C_base

                    try:
                        T = ctrl.feedback(C_full * G, 1)
                        p_lc = ctrl.poles(T)
                        if any(np.real(p_lc) > 1e-6):
                            continue
                    except Exception:
                        continue

                    try:
                        met = calcular_metricas_temporales(G, C_full, t_final=40)
                        ts_v = met['ts']
                        Mp_v = met['Mp']
                        if ts_v > SPECS['ts_max'] or Mp_v > SPECS['Mp_max']:
                            continue
                    except Exception:
                        continue

                    mejores.append({
                        'z_L':   round(z_L, 4),
                        'p_L':   round(p_L, 4),
                        'z_i':   round(z_i, 4),
                        'Kc':    round(Kc, 4),
                        'C':     C_full,
                        'ts':    round(ts_v, 3),
                        'Mp':    round(Mp_v, 3),
                        'sd':    sd,
                    })

    if not mejores:
        raise RuntimeError("Grid search lead+PI no encontró compensador factible.")

    # Prioridad: menor Mp, desempate por ts
    mejores.sort(key=lambda x: (x['Mp'], x['ts']))
    best = mejores[0]

    if verbose:
        sep = '=' * 62
        print(f"\n{sep}")
        print(f"  GRID SEARCH — Compensador Lead + PI (Tipo 2)")
        print(sep)
        print(f"  Lead: cero = -{best['z_L']:.4f},  polo = -{best['p_L']:.4f}")
        print(f"  PI:   cero = -{best['z_i']:.4f},  polo = 0 (integrador)")
        print(f"  Kc   = {best['Kc']:.4f}")
        print(f"  ts   = {best['ts']:.2f} s   (máx {SPECS['ts_max']} s)  ✓")
        print(f"  Mp   = {best['Mp']:.2f} %  (máx {SPECS['Mp_max']} %)  ✓")
        print(f"  e_pert = 0  (sistema Tipo 2)  ✓")
        print(f"  e_rampa = 0  (sistema Tipo 2)  ✓")
        print(sep)

    return best, mejores


def evaluar_specs(G, C, verbose=True, t_final=40):
    """
    Evalúa si el sistema compensado cumple las especificaciones.

    Args:
        G: planta
        C: controlador
        verbose: imprimir resultados

    Returns:
        dict con resultados y booleano 'cumple'
    """
    from analisis import calcular_metricas_temporales, analizar_error_estacionario

    T = ctrl.feedback(C * G, 1)
    metricas = calcular_metricas_temporales(G, C, t_final=t_final)

    # Error ante perturbación de 0.2
    G_n = crear_planta_perturbacion()
    T_pert = G_n * ctrl.feedback(ctrl.tf(1, 1), C * G)
    t = np.linspace(0, t_final, 5000)
    _, y_pert = ctrl.step_response(T_pert, t)
    y_pert = SPECS['amp_pert'] * np.array(y_pert).flatten()
    error_pert_ss = np.abs(y_pert[-1])  # error en régimen permanente por perturbación
    error_pert_pct = error_pert_ss * 100  # como % de la referencia unitaria

    cumple_ts = metricas['ts'] <= SPECS['ts_max']
    cumple_Mp = metricas['Mp'] <= SPECS['Mp_max']
    cumple_pert = error_pert_pct <= SPECS['error_pert_max']
    cumple = cumple_ts and cumple_Mp and cumple_pert

    if verbose:
        print(f"\n{'='*50}")
        print(f"  EVALUACIÓN DE ESPECIFICACIONES")
        print(f"{'='*50}")
        print(f"  ts  = {metricas['ts']:.2f} s    (máx {SPECS['ts_max']} s)    {'✓' if cumple_ts else '✗'}")
        print(f"  Mp  = {metricas['Mp']:.2f}%   (máx {SPECS['Mp_max']}%)   {'✓' if cumple_Mp else '✗'}")
        print(f"  e_pert = {error_pert_pct:.2f}% (máx {SPECS['error_pert_max']}%)   {'✓' if cumple_pert else '✗'}")
        print(f"  RESULTADO: {'CUMPLE ✓' if cumple else 'NO CUMPLE ✗'}")
        print(f"{'='*50}")

    return {
        'ts': metricas['ts'],
        'Mp': metricas['Mp'],
        'error_pert_ss': error_pert_ss,
        'error_pert_pct': error_pert_pct,
        'cumple': cumple,
        'cumple_ts': cumple_ts,
        'cumple_Mp': cumple_Mp,
        'cumple_pert': cumple_pert,
    }


# =============================================================================
# Parte B.c - Comparación Bode compensado vs no compensado
# =============================================================================
def comparar_bode(G, C, titulo="Comparación Bode: Compensado vs No Compensado"):
    """
    Grafica Bode del sistema compensado y no compensado en el mismo gráfico.

    Args:
        G: planta
        C: controlador
        titulo: título del gráfico
    """
    omega = np.logspace(-2, 3, 1000)

    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(10, 7))

    for sys, label, color, ls in [(G, 'Sin compensar (C=1)', 'b', '--'),
                                    (C * G, 'Compensado', 'r', '-')]:
        resp = ctrl.frequency_response(sys, omega)
        mag_dB = 20 * np.log10(np.abs(resp.fresp[0, 0, :]))
        phase_deg = np.degrees(np.unwrap(np.angle(resp.fresp[0, 0, :])))

        ax_mag.semilogx(omega, mag_dB, color=color, linestyle=ls, linewidth=1.5, label=label)
        ax_phase.semilogx(omega, phase_deg, color=color, linestyle=ls, linewidth=1.5, label=label)

    ax_mag.axhline(0, color='k', linestyle=':', linewidth=0.8)
    ax_mag.set_ylabel('Magnitud [dB]')
    ax_mag.set_title(titulo)
    ax_mag.legend(loc='best')

    ax_phase.axhline(-180, color='k', linestyle=':', linewidth=0.8)
    ax_phase.set_ylabel('Fase [°]')
    ax_phase.set_xlabel('Frecuencia [rad/s]')
    ax_phase.legend(loc='best')

    # Anotar márgenes
    for sys, label, color in [(G, 'C=1', 'b'), (C * G, 'Comp.', 'r')]:
        gm, pm, wcg, wcp = ctrl.margin(sys)
        gm_dB = 20 * np.log10(gm) if np.isfinite(gm) else float('inf')
        print(f"  {label}: GM = {gm_dB:.1f} dB, PM = {pm:.1f}°, ωcg = {wcg:.2f}, ωcp = {wcp:.2f}")

    plt.tight_layout()
    return fig


# =============================================================================
# Parte B.d/e - Seguimiento de rampa y adición de integrador
# =============================================================================
def agregar_integrador(C):
    """
    Agrega un integrador al controlador para lograr seguimiento de rampa:
        C_new(s) = C(s) · (1/s)

    Esto incrementa el tipo del sistema en 1:
    - Sistema Tipo 1 → Tipo 2: error nulo ante rampa.

    Riesgo: el integrador adicional reduce el margen de fase.
    Hay que verificar que el sistema siga siendo estable.

    Args:
        C: controlador original

    Returns:
        C_new: controlador con integrador
    """
    integrador = ctrl.tf([1], [1, 0])
    C_new = C * integrador

    print("--- Integrador agregado ---")
    print(f"  C_new(s) = C(s) / s")
    print(f"  Nuevo tipo del sistema: 2 (seguimiento perfecto de rampa)")

    return C_new


def graficar_respuesta_rampa_compensada(G, C, C_int=None, t_final=30,
                                         titulo="Seguimiento de rampa"):
    """
    Compara la respuesta ante rampa del sistema con y sin integrador.

    Args:
        G: planta
        C: controlador sin integrador
        C_int: controlador con integrador (None → no se grafica)
        t_final: tiempo de simulación
    """
    t = np.linspace(0, t_final, 3000)
    ref = t  # rampa unitaria

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t, ref, 'k--', linewidth=1, label='Referencia (rampa)')

    for Ci, label, color in [(C, 'Sin integrador', 'b'),
                              (C_int, 'Con integrador', 'r')]:
        if Ci is None:
            continue
        T = ctrl.feedback(Ci * G, 1)
        _, y = ctrl.forced_response(T, t, ref)
        y = np.array(y).flatten()
        ax.plot(t, y, color=color, linewidth=1.5, label=label)

    ax.set_xlabel('Tiempo [s]')
    ax.set_ylabel('Posición angular θ [rad]')
    ax.set_title(titulo)
    ax.legend(loc='best')
    ax.set_xlim([0, t_final])
    plt.tight_layout()
    return fig
