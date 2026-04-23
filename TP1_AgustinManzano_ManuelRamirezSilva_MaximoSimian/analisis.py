"""
analisis.py - Funciones de análisis del sistema (Parte A, Ejercicio 1)

Contiene funciones para:
- Cálculo de márgenes de ganancia y fase
- Análisis de perturbaciones (escalón y rampa)
- Simulación del sistema en lazo cerrado
- Gráficos de Bode, respuesta temporal, etc.

Curso: IOR442 Sistemas de Control - Dr. Mariano Scaramal
"""

import control as ctrl
import numpy as np
import matplotlib.pyplot as plt
from planta import crear_planta, crear_planta_perturbacion, crear_planta_ss

# Configuración global de matplotlib para gráficos consistentes
plt.rcParams.update({
    'figure.figsize': (10, 6),
    'font.size': 11,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'lines.linewidth': 1.5,
})


# =============================================================================
# Parte A.b - Márgenes de estabilidad
# =============================================================================
def calcular_margenes(G, titulo="Planta G(s)"):
    """
    Calcula y muestra los márgenes de ganancia y fase de un sistema.

    Args:
        G: función de transferencia (lazo abierto)
        titulo: nombre para mostrar en el resultado

    Returns:
        dict con gm (dB), pm (grados), wcg (rad/s), wcp (rad/s)
    """
    gm, pm, wcg, wcp = ctrl.margin(G)
    gm_dB = 20 * np.log10(gm) if gm and np.isfinite(gm) else float('inf')

    print(f"--- Márgenes de estabilidad: {titulo} ---")
    print(f"  Margen de ganancia: {gm_dB:.2f} dB  (a ω = {wcg:.4f} rad/s)")
    print(f"  Margen de fase:     {pm:.2f}°       (a ω = {wcp:.4f} rad/s)")

    return {"gm": gm, "gm_dB": gm_dB, "pm": pm, "wcg": wcg, "wcp": wcp}


def graficar_bode_con_margenes(G, titulo="Diagrama de Bode"):
    """
    Grafica el diagrama de Bode con indicación de márgenes.

    Args:
        G: función de transferencia
        titulo: título del gráfico
    """
    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(10, 7))

    # Calcular Bode manualmente para tener control del gráfico
    omega = np.logspace(-2, 3, 1000)
    mag, phase, omega_out = ctrl.frequency_response(G, omega)

    mag_dB = 20 * np.log10(np.abs(mag))
    phase_deg = np.degrees(np.angle(mag * np.exp(1j * np.radians(phase)) / np.abs(mag)))
    # Usar la fase directamente del frequency_response
    phase_deg = np.degrees(np.unwrap(np.angle(mag * np.exp(1j * 0))))

    # Recalcular correctamente
    resp = ctrl.frequency_response(G, omega)
    mag_dB = 20 * np.log10(np.abs(resp.fresp[0, 0, :]))
    phase_deg = np.degrees(np.unwrap(np.angle(resp.fresp[0, 0, :])))

    # Magnitud
    ax_mag.semilogx(omega_out, mag_dB, 'b-')
    ax_mag.axhline(0, color='k', linestyle='--', linewidth=0.8)
    ax_mag.set_ylabel('Magnitud [dB]')
    ax_mag.set_title(titulo)

    # Fase
    ax_phase.semilogx(omega_out, phase_deg, 'r-')
    ax_phase.axhline(-180, color='k', linestyle='--', linewidth=0.8)
    ax_phase.set_ylabel('Fase [°]')
    ax_phase.set_xlabel('Frecuencia [rad/s]')

    # Marcar márgenes
    gm, pm, wcg, wcp = ctrl.margin(G)
    if np.isfinite(wcp) and wcp > 0:
        ax_mag.axvline(wcp, color='g', linestyle=':', alpha=0.7, label=f'ωcp = {wcp:.2f} rad/s')
        ax_phase.axvline(wcp, color='g', linestyle=':', alpha=0.7)
        # Buscar fase en wcp
        idx_cp = np.argmin(np.abs(omega_out - wcp))
        ax_phase.annotate(f'PM = {pm:.1f}°',
                          xy=(wcp, phase_deg[idx_cp]),
                          xytext=(wcp * 3, phase_deg[idx_cp] + 20),
                          arrowprops=dict(arrowstyle='->', color='green'),
                          color='green', fontsize=10)

    if np.isfinite(wcg) and wcg > 0:
        gm_dB = 20 * np.log10(gm)
        ax_mag.axvline(wcg, color='m', linestyle=':', alpha=0.7, label=f'ωcg = {wcg:.2f} rad/s')
        ax_phase.axvline(wcg, color='m', linestyle=':', alpha=0.7)
        idx_cg = np.argmin(np.abs(omega_out - wcg))
        ax_mag.annotate(f'GM = {gm_dB:.1f} dB',
                        xy=(wcg, mag_dB[idx_cg]),
                        xytext=(wcg * 3, mag_dB[idx_cg] + 10),
                        arrowprops=dict(arrowstyle='->', color='purple'),
                        color='purple', fontsize=10)

    ax_mag.legend(loc='best')
    plt.tight_layout()
    return fig


# =============================================================================
# Parte A.c - Respuesta ante escalón con perturbación
# =============================================================================
def simular_escalon_con_perturbacion(G, C=None, amp_ref=1.0, amp_pert=0.1,
                                      t_final=30, n_puntos=2000):
    """
    Simula el sistema en lazo cerrado ante referencia escalón unitario
    y perturbación de carga constante (simultáneas desde t=0).

    El diagrama de bloques es:
        R(s) →(+)→ C(s) → G_motor → (+) → 1/s → Θ(s)
               -↑                     ↑ N(s)        |
                └─────────────────────────────────────┘

    En lazo cerrado:
        Θ(s) = [C(s)·G(s)/(1+C(s)·G(s))]·R(s) + [G_n(s)/(1+C(s)·G(s))]·N(s)

    Args:
        G: planta G(s)
        C: controlador C(s), None → C(s) = 1
        amp_ref: amplitud de la referencia escalón
        amp_pert: amplitud de la perturbación
        t_final: tiempo de simulación [s]
        n_puntos: cantidad de puntos temporales

    Returns:
        t: vector de tiempo
        y_total: salida total (referencia + perturbación)
        y_ref: componente debida a la referencia
        y_pert: componente debida a la perturbación
    """
    if C is None:
        C = ctrl.tf([1], [1])

    G_n_ol = crear_planta_perturbacion()

    # Transferencias en lazo cerrado
    CG = C * G
    T_ref = ctrl.feedback(CG, 1)                    # Θ/R = CG/(1+CG)
    T_pert = G_n_ol * ctrl.feedback(ctrl.tf(1, 1), CG)  # Θ/N = G_n/(1+CG)

    t = np.linspace(0, t_final, n_puntos)

    # Respuesta a cada entrada por separado (superposición)
    t_r, y_ref = ctrl.step_response(T_ref, t)
    t_n, y_pert = ctrl.step_response(T_pert, t)

    y_ref = amp_ref * np.array(y_ref).flatten()
    y_pert = amp_pert * np.array(y_pert).flatten()
    y_total = y_ref + y_pert

    return t, y_total, y_ref, y_pert


def graficar_respuesta_escalon_pert(t, y_total, y_ref, y_pert, amp_ref=1.0,
                                     titulo="Respuesta ante escalón + perturbación"):
    """
    Grafica la respuesta temporal mostrando referencia, salida total
    y las componentes individuales.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(t, np.ones_like(t) * amp_ref, 'k--', linewidth=1, label='Referencia')
    ax.plot(t, y_total, 'b-', linewidth=2, label='Salida total θ(t)')
    ax.plot(t, y_ref, 'g--', linewidth=1, alpha=0.7, label='Componente referencia')
    ax.plot(t, y_pert, 'r--', linewidth=1, alpha=0.7, label='Componente perturbación')

    ax.set_xlabel('Tiempo [s]')
    ax.set_ylabel('Posición angular θ [rad]')
    ax.set_title(titulo)
    ax.legend(loc='best')
    ax.set_xlim([0, t[-1]])
    plt.tight_layout()
    return fig


# =============================================================================
# Parte A.d - Respuesta ante rampa con perturbación
# =============================================================================
def simular_rampa_con_perturbacion(G, C=None, amp_pert=0.1,
                                    t_final=30, n_puntos=2000):
    """
    Simula el sistema en lazo cerrado ante referencia rampa unitaria r(t)=t
    y perturbación de carga constante.

    Para la rampa usamos simulación con forced_response ya que
    step_response solo sirve para escalón.

    Args:
        G: planta G(s)
        C: controlador C(s), None → C(s) = 1
        amp_pert: amplitud de la perturbación
        t_final: tiempo de simulación [s]
        n_puntos: cantidad de puntos temporales

    Returns:
        t: vector de tiempo
        y_total: salida total
        y_rampa: componente debida a la rampa
        y_pert: componente debida a la perturbación
        ref: señal de referencia rampa
    """
    if C is None:
        C = ctrl.tf([1], [1])

    G_n_ol = crear_planta_perturbacion()

    CG = C * G
    T_ref = ctrl.feedback(CG, 1)
    T_pert = G_n_ol * ctrl.feedback(ctrl.tf(1, 1), CG)

    t = np.linspace(0, t_final, n_puntos)
    ref = t  # rampa unitaria r(t) = t

    # Respuesta a la rampa (forced_response)
    t_r, y_rampa = ctrl.forced_response(T_ref, t, ref)
    y_rampa = np.array(y_rampa).flatten()

    # Respuesta a la perturbación (escalón de amplitud amp_pert)
    t_n, y_pert = ctrl.step_response(T_pert, t)
    y_pert = amp_pert * np.array(y_pert).flatten()

    y_total = y_rampa + y_pert

    return t, y_total, y_rampa, y_pert, ref


def graficar_respuesta_rampa_pert(t, y_total, y_rampa, y_pert, ref,
                                   titulo="Respuesta ante rampa + perturbación"):
    """
    Grafica la respuesta ante rampa mostrando referencia, salida y error.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), height_ratios=[3, 1])

    # Respuesta temporal
    ax1.plot(t, ref, 'k--', linewidth=1, label='Referencia (rampa)')
    ax1.plot(t, y_total, 'b-', linewidth=2, label='Salida total θ(t)')
    ax1.plot(t, y_rampa, 'g--', linewidth=1, alpha=0.7, label='Componente rampa')
    ax1.plot(t, y_pert, 'r--', linewidth=1, alpha=0.7, label='Componente perturbación')

    ax1.set_ylabel('Posición angular θ [rad]')
    ax1.set_title(titulo)
    ax1.legend(loc='best')
    ax1.set_xlim([0, t[-1]])

    # Error
    error = ref - y_total
    ax2.plot(t, error, 'm-', linewidth=1.5)
    ax2.set_xlabel('Tiempo [s]')
    ax2.set_ylabel('Error [rad]')
    ax2.set_title('Error de seguimiento e(t) = r(t) - θ(t)')
    ax2.set_xlim([0, t[-1]])

    plt.tight_layout()
    return fig


# =============================================================================
# Parte A.e - Análisis para proponer controlador
# =============================================================================
def analizar_error_estacionario(G, C=None):
    """
    Calcula el error en régimen permanente para distintas entradas.

    Para un sistema con lazo abierto L(s) = C(s)·G(s):
    - Escalón: e_ss = 1 / (1 + Kp), donde Kp = lim_{s→0} L(s)
    - Rampa:   e_ss = 1 / Kv, donde Kv = lim_{s→0} s·L(s)

    Args:
        G: planta
        C: controlador (None → C=1)

    Returns:
        dict con Kp, Kv, error_escalon, error_rampa
    """
    if C is None:
        C = ctrl.tf([1], [1])

    L = C * G

    # Kp = lim_{s→0} L(s) → para sistema Tipo 1, Kp = ∞
    # Kv = lim_{s→0} s·L(s)
    s = ctrl.tf([1, 0], [1])
    sL = ctrl.minreal(s * L)

    # Evaluar en s = 0 (o muy cerca)
    Kp = np.abs(ctrl.evalfr(L, 1e-10))
    Kv = np.abs(ctrl.evalfr(sL, 1e-10))

    error_escalon = 1 / (1 + Kp) if np.isfinite(Kp) else 0
    error_rampa = 1 / Kv if Kv > 0 else float('inf')

    print("--- Errores en régimen permanente ---")
    print(f"  Kp = {Kp:.4f} → Error ante escalón: {error_escalon:.6f}")
    print(f"  Kv = {Kv:.4f} → Error ante rampa:   {error_rampa:.6f}")

    return {"Kp": Kp, "Kv": Kv, "error_escalon": error_escalon, "error_rampa": error_rampa}


def calcular_metricas_temporales(G, C=None, t_final=30):
    """
    Calcula sobrepico (Mp), tiempo de establecimiento (ts) y valor final
    de la respuesta al escalón en lazo cerrado.

    Args:
        G: planta
        C: controlador (None → C=1)
        t_final: tiempo de simulación

    Returns:
        dict con Mp (%), ts (s), y_final
    """
    if C is None:
        C = ctrl.tf([1], [1])

    T = ctrl.feedback(C * G, 1)
    t = np.linspace(0, t_final, 5000)
    t_out, y = ctrl.step_response(T, t)
    y = np.array(y).flatten()

    y_final = y[-1]

    # Sobrepico
    y_max = np.max(y)
    if y_final > 0:
        Mp = (y_max - y_final) / y_final * 100
    else:
        Mp = 0

    # Tiempo de establecimiento (2%)
    tolerancia = 0.02
    banda_sup = y_final * (1 + tolerancia)
    banda_inf = y_final * (1 - tolerancia)

    # Buscar el último instante donde sale de la banda
    fuera_banda = np.where((y > banda_sup) | (y < banda_inf))[0]
    if len(fuera_banda) > 0:
        ts = t_out[fuera_banda[-1]]
    else:
        ts = 0

    print(f"--- Métricas temporales (escalón unitario) ---")
    print(f"  Valor final:               y_∞ = {y_final:.4f}")
    print(f"  Sobrepico:                 Mp  = {Mp:.2f}%")
    print(f"  Tiempo de establecimiento: ts  = {ts:.2f} s (criterio 2%)")

    return {"Mp": Mp, "ts": ts, "y_final": y_final, "t": t_out, "y": y}
