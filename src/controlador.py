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


def disenar_compensador_lead(G, polo_comp, cero_comp, K_comp=1.0):
    """
    Crea un compensador lead de la forma:
        C(s) = K_comp · (s - cero_comp) / (s - polo_comp)

    donde |polo_comp| > |cero_comp| (polo más a la izquierda que el cero).

    El compensador lead aporta fase positiva → mejora transitorio.

    Args:
        G: planta (para verificación)
        polo_comp: polo del compensador (valor negativo real)
        cero_comp: cero del compensador (valor negativo real)
        K_comp: ganancia del compensador

    Returns:
        C: función de transferencia del compensador
    """
    C = K_comp * ctrl.tf([1, -cero_comp], [1, -polo_comp])

    print(f"--- Compensador Lead ---")
    print(f"  C(s) = {K_comp} · (s + {-cero_comp:.4f}) / (s + {-polo_comp:.4f})")
    print(f"  Cero en s = {cero_comp:.4f}")
    print(f"  Polo en s = {polo_comp:.4f}")

    return C


def ajustar_ganancia_rlocus(G, C, polos_deseados):
    """
    Ajusta la ganancia K del controlador para que el root locus
    pase por los polos deseados (o lo más cerca posible).

    Usa la condición de magnitud del root locus:
        |C(s_d)·G(s_d)| = 1  →  K = 1 / |C(s_d)·G(s_d)|

    Args:
        G: planta
        C: compensador (sin ganancia final)
        polos_deseados: polo dominante deseado (complejo)

    Returns:
        K: ganancia necesaria
        C_final: compensador con ganancia ajustada
    """
    s_d = polos_deseados

    # Evaluar CG en el polo deseado
    CG = C * G
    CG_eval = ctrl.evalfr(CG, s_d)
    K = 1.0 / np.abs(CG_eval)

    C_final = K * C

    print(f"  Polo deseado: s_d = {s_d}")
    print(f"  |CG(s_d)| = {np.abs(CG_eval):.6f}")
    print(f"  Ganancia ajustada: K = {K:.4f}")

    return K, C_final


def evaluar_specs(G, C, verbose=True):
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
    metricas = calcular_metricas_temporales(G, C, t_final=40)

    # Error ante perturbación de 0.2
    G_n = crear_planta_perturbacion()
    T_pert = G_n * ctrl.feedback(ctrl.tf(1, 1), C * G)
    t = np.linspace(0, 40, 5000)
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
