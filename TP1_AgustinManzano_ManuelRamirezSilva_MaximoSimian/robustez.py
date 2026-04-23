"""
robustez.py - Análisis de robustez y diseño de compensadores (Ejercicio 2)

Contiene funciones para:
- Diagrama de Nyquist y análisis de estabilidad
- Simulación de sistemas LTV (cambio de B en t=15s)
- Root locus comparativo (nominal vs falla)
- Diseño de compensadores lead normalizados (ganancia DC = 1)
- Comparación de Bode para múltiples compensadores

Curso: IOR442 Sistemas de Control - Dr. Mariano Scaramal
"""

import control as ctrl
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from planta2 import crear_planta2, crear_planta2_ss, K, J, B_NOM, B_FALLA

plt.rcParams.update({
    'figure.figsize': (10, 6),
    'font.size': 11,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'lines.linewidth': 1.5,
})


# =============================================================================
# Parte A.a - Diagrama de Nyquist
# =============================================================================
def graficar_nyquist(G, titulo="Diagrama de Nyquist", ax=None):
    """
    Grafica el diagrama de Nyquist de G(s) con el punto crítico (-1,0).

    Para sistemas con polo en el origen, python-control indenta el contorno
    automáticamente. Usamos indent_direction='right' para evitar rodear el polo.

    Args:
        G: función de transferencia (lazo abierto)
        titulo: título del gráfico
        ax: eje existente (None → crear nuevo)

    Returns:
        fig: figura de matplotlib (o None si se pasó ax)
        count: número de encirculamientos de -1
    """
    fig_created = ax is None
    if fig_created:
        fig, ax = plt.subplots(figsize=(8, 8))
    else:
        fig = None

    # Nyquist usando python-control
    count = ctrl.nyquist_plot(G, omega_limits=[1e-4, 1e3],
                               plot=True, ax=ax,
                               indent_direction='right')

    ax.plot(-1, 0, 'rx', markersize=12, markeredgewidth=2, label='Punto crítico (-1,0)')
    ax.set_title(titulo)
    ax.legend(loc='best')
    ax.set_xlabel('Re')
    ax.set_ylabel('Im')

    if fig_created:
        plt.tight_layout()

    return fig, count


def analizar_estabilidad_nyquist(G, B_val):
    """
    Analiza la estabilidad usando el criterio de Nyquist.

    Criterio: N = Z - P
    - P = número de polos de G(s) en el semiplano derecho (SPD)
    - N = número de encirculamientos de (-1,0) en sentido horario
    - Z = polos de lazo cerrado en SPD = N + P

    Para estabilidad en LC: Z = 0 → N = -P

    Args:
        G: función de transferencia de lazo abierto
        B_val: valor de B (para contexto)

    Returns:
        dict con P, N, Z y si es estable
    """
    polos_la = ctrl.poles(G)
    P = np.sum(np.real(polos_la) > 0)

    # Encirculamientos (calculamos usando la función de python-control)
    count = ctrl.nyquist_plot(G, omega_limits=[1e-4, 1e3], plot=False,
                               indent_direction='right')

    # count retorna el número de encirculamientos CW de -1
    # En python-control, count = N (encirculamientos horarios)
    N = count if isinstance(count, (int, float)) else 0
    Z = int(N) + P

    estable = (Z == 0)

    print(f"--- Criterio de Nyquist (B = {B_val}) ---")
    print(f"  Polos de LA en SPD (P): {P}")
    print(f"  Encirculamientos CW (N): {N}")
    print(f"  Polos de LC en SPD (Z = N + P): {Z}")
    print(f"  Estable en LC: {'SÍ' if estable else 'NO'}")

    return {"P": P, "N": int(N), "Z": int(Z), "estable": estable}


# =============================================================================
# Parte A.b - Simulación LTV (falla en t=15s)
# =============================================================================
def simular_falla(C_tf=None, t_falla=15, t_final=40, B_pre=B_NOM, B_post=B_FALLA):
    """
    Simula el sistema en lazo cerrado con realimentación unitaria ante escalón,
    donde B cambia de B_pre a B_post en t = t_falla.

    Como el sistema es variante en el tiempo (LTV), no podemos usar
    funciones de transferencia directamente. Usamos espacio de estados
    con scipy.integrate.solve_ivp.

    El lazo cerrado con C(s) y realimentación unitaria se resuelve
    en espacio de estados aumentado [planta + controlador].

    Para C(s) = 1 (sin compensar):
        e = r - θ,  u = e
        dθ/dt = ω
        dω/dt = -(B(t)/J)·ω + (K/J)·u = -(B(t)/J)·ω + (K/J)·(r - θ)

    Para C(s) = Kc·(s+z)/(s+p) (lead):
        Necesitamos realizar el controlador en espacio de estados
        y combinar con la planta.

    Args:
        C_tf: controlador como TransferFunction (None → C=1)
        t_falla: instante de la falla [s]
        t_final: tiempo total de simulación [s]
        B_pre: valor de B antes de la falla
        B_post: valor de B después de la falla

    Returns:
        t: vector de tiempo
        theta: salida θ(t)
        omega: velocidad ω(t)
    """
    r = 1.0  # referencia escalón unitario
    n_puntos = 5000
    t_span = (0, t_final)
    t_eval = np.linspace(0, t_final, n_puntos)

    if C_tf is None or (len(ctrl.poles(C_tf)) == 0):
        # Controlador proporcional (C=K constante)
        # Extraer ganancia DC
        Kc = float(np.abs(ctrl.evalfr(C_tf, 1e-10))) if C_tf is not None else 1.0

        def dxdt(t, x):
            theta, omega = x
            B_t = B_pre if t < t_falla else B_post
            e = r - theta
            u = Kc * e
            dtheta = omega
            domega = -(B_t / J) * omega + (K / J) * u
            return [dtheta, domega]

        x0 = [0, 0]
        sol = solve_ivp(dxdt, t_span, x0, t_eval=t_eval, method='RK45',
                        max_step=0.01)
        return sol.t, sol.y[0], sol.y[1]

    else:
        # Controlador dinámico: C(s) = num/den
        # Realizamos C(s) en espacio de estados y aumentamos con la planta
        C_ss = ctrl.tf2ss(C_tf)
        Ac = np.array(C_ss.A)
        Bc = np.array(C_ss.B)
        Cc = np.array(C_ss.C)
        Dc = np.array(C_ss.D)
        nc = Ac.shape[0]  # orden del controlador

        # Estado aumentado: x = [θ, ω, xc1, xc2, ...]
        # donde xc son los estados del controlador
        def dxdt(t, x):
            theta = x[0]
            omega = x[1]
            xc = x[2:2+nc].reshape(-1, 1)

            B_t = B_pre if t < t_falla else B_post
            e = r - theta

            # Salida del controlador: u = Cc·xc + Dc·e
            u = (Cc @ xc + Dc * e).item()

            # Derivadas de la planta
            dtheta = omega
            domega = -(B_t / J) * omega + (K / J) * u

            # Derivadas del controlador: dxc = Ac·xc + Bc·e
            dxc = (Ac @ xc + Bc * e).flatten()

            return [dtheta, domega] + list(dxc)

        x0 = [0, 0] + [0] * nc
        sol = solve_ivp(dxdt, t_span, x0, t_eval=t_eval, method='RK45',
                        max_step=0.01)
        return sol.t, sol.y[0], sol.y[1]


def graficar_respuesta_falla(t, theta, titulo="Respuesta con falla", ax=None,
                              label=None, color=None, t_falla=15):
    """Grafica la respuesta temporal marcando el instante de falla."""
    fig_created = ax is None
    if fig_created:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = None

    kwargs = {}
    if label: kwargs['label'] = label
    if color: kwargs['color'] = color

    ax.plot(t, theta, linewidth=1.5, **kwargs)
    ax.axhline(1.0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axvline(t_falla, color='r', linestyle=':', linewidth=1, alpha=0.7,
               label=f'Falla (t={t_falla}s)' if fig_created else '')

    ax.set_xlabel('Tiempo [s]')
    ax.set_ylabel('$\\theta$ [rad]')
    ax.set_title(titulo)
    if fig_created:
        ax.legend()
    ax.set_xlim([0, t[-1]])

    if fig_created:
        plt.tight_layout()
    return fig


# =============================================================================
# Parte A.c - Root Locus comparativo
# =============================================================================
def graficar_rlocus_comparativo(B_vals, labels=None, colores=None):
    """
    Grafica el root locus para varios valores de B en el mismo gráfico.

    Args:
        B_vals: lista de valores de B
        labels: lista de etiquetas
        colores: lista de colores
    """
    if labels is None:
        labels = [f'B = {b}' for b in B_vals]
    if colores is None:
        colores = ['b', 'r', 'g', 'm']

    fig, ax = plt.subplots(figsize=(10, 8))

    for B_val, label, color in zip(B_vals, labels, colores):
        G = crear_planta2(B_val)
        rlist, klist = ctrl.root_locus(G, plot=False)

        for i in range(rlist.shape[1]):
            lbl = label if i == 0 else None
            ax.plot(np.real(rlist[:, i]), np.imag(rlist[:, i]),
                    color=color, linewidth=1.2, label=lbl)

        # Marcar polos de lazo abierto
        polos = ctrl.poles(G)
        ax.plot(np.real(polos), np.imag(polos), 'x', color=color,
                markersize=10, markeredgewidth=2)

        # Marcar polos de lazo cerrado (K=1, realimentación unitaria)
        T = ctrl.feedback(G, 1)
        polos_lc = ctrl.poles(T)
        ax.plot(np.real(polos_lc), np.imag(polos_lc), '*', color=color,
                markersize=12, markeredgewidth=1.5)

    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)
    ax.set_xlabel('Re(s)')
    ax.set_ylabel('Im(s)')
    ax.set_title('Root Locus comparativo')
    ax.legend(loc='best')
    plt.tight_layout()
    return fig


# =============================================================================
# Parte B.a - Diseño de compensadores lead normalizados
# =============================================================================
def disenar_lead_por_margen_fase(G, PM_deseado, B_val=B_NOM):
    """
    Diseña un compensador lead para lograr un margen de fase deseado.

    Compensador lead:
        C(s) = Kc · (s + z) / (s + p),   p > z > 0

    Normalización: ganancia en DC = 1 → C(0) = Kc·z/p = 1 → Kc = p/z

    Esto es posible porque un lead tiene |C(0)| = z/p < 1, y al multiplicar
    por Kc = p/z se normaliza a 1. La ventaja es que no altera el error
    de posición del sistema.

    Procedimiento de diseño por respuesta en frecuencia:
    1. Calcular PM actual del sistema
    2. Fase adicional necesaria: φ_max = PM_deseado - PM_actual + margen_seguridad
    3. α = (1 - sin(φ_max)) / (1 + sin(φ_max))
    4. Frecuencia de máxima fase: ωm donde |G(jωm)| = -10·log10(α) dB
       (esto asegura que ωm sea la nueva frecuencia de cruce)
    5. z = ωm · √α,  p = ωm / √α
    6. Kc = p / z = 1/α  (para normalizar DC)

    Args:
        G: planta
        PM_deseado: margen de fase objetivo [grados]
        B_val: valor de B (para info)

    Returns:
        C: compensador lead (TransferFunction)
        info: dict con parámetros del diseño
    """
    # PM actual
    gm, pm_actual, wcg, wcp = ctrl.margin(G)

    omega = np.logspace(-3, 3, 10000)
    resp  = ctrl.frequency_response(G, omega)
    mag   = np.abs(resp.fresp[0, 0, :])

    # Iteración: aumenta el margen de seguridad hasta que PM_obtenido >= PM_deseado
    C = pm_c = alpha = wm = z = p = Kc = None
    for extra in range(5, 80, 3):
        phi_extra   = PM_deseado - pm_actual + extra
        phi_max_rad = np.radians(phi_extra)
        sin_phi     = np.sin(phi_max_rad)
        alpha_i     = (1 - sin_phi) / (1 + sin_phi)
        target_mag  = np.sqrt(alpha_i)

        idx  = np.argmin(np.abs(mag - target_mag))
        wm_i = omega[idx]
        z_i  = wm_i * np.sqrt(alpha_i)
        p_i  = wm_i / np.sqrt(alpha_i)
        Kc_i = p_i / z_i

        C_i  = Kc_i * ctrl.tf([1, z_i], [1, p_i])
        _, pm_c_i, _, _ = ctrl.margin(C_i * G)

        if pm_c_i >= PM_deseado:
            C, pm_c, alpha, wm, z, p, Kc = C_i, pm_c_i, alpha_i, wm_i, z_i, p_i, Kc_i
            break

    # Fallback: usar el resultado de la última iteración si ninguna alcanzó el target
    if C is None:
        C, pm_c, alpha, wm, z, p, Kc = C_i, pm_c_i, alpha_i, wm_i, z_i, p_i, Kc_i

    info = {
        'PM_deseado': PM_deseado,
        'PM_actual':  pm_actual,
        'PM_obtenido': pm_c,
        'alpha': alpha,
        'wm': wm,
        'z': z,
        'p': p,
        'Kc': Kc,
        'ganancia_DC': float(np.abs(ctrl.evalfr(C, 1e-10))),
    }

    print(f"--- Compensador Lead para PM = {PM_deseado}° (B = {B_val}) ---")
    print(f"  PM sin compensar: {pm_actual:.1f}°")
    print(f"  α = {alpha:.4f}")
    print(f"  ωm = {wm:.4f} rad/s")
    print(f"  Cero: z = {z:.4f},  Polo: p = {p:.4f}")
    print(f"  Kc = {Kc:.4f} (ganancia DC = {info['ganancia_DC']:.4f})")
    print(f"  PM obtenido: {pm_c:.1f}°")
    print(f"  C(s) = {Kc:.4f} · (s + {z:.4f}) / (s + {p:.4f})")

    return C, info


def comparar_bode_multiples(G, compensadores, labels, colores=None,
                             titulo="Bode comparativo"):
    """
    Grafica Bode de G(s) compensada con múltiples controladores.

    Args:
        G: planta
        compensadores: lista de C(s)
        labels: lista de etiquetas
        colores: lista de colores
    """
    if colores is None:
        colores = ['b', 'r', 'g', 'm', 'c']

    omega = np.logspace(-3, 3, 2000)
    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(10, 7))

    for C, label, color in zip(compensadores, labels, colores):
        CG = C * G
        resp = ctrl.frequency_response(CG, omega)
        mag_dB = 20 * np.log10(np.abs(resp.fresp[0, 0, :]))
        phase_deg = np.degrees(np.unwrap(np.angle(resp.fresp[0, 0, :])))

        ax_mag.semilogx(omega, mag_dB, color=color, linewidth=1.5, label=label)
        ax_phase.semilogx(omega, phase_deg, color=color, linewidth=1.5, label=label)

        # Marcar PM en fase
        gm, pm, wcg, wcp = ctrl.margin(CG)
        if np.isfinite(wcp) and wcp > 0:
            ax_phase.axvline(wcp, color=color, linestyle=':', alpha=0.3)

    ax_mag.axhline(0, color='k', linestyle=':', linewidth=0.8)
    ax_mag.set_ylabel('Magnitud [dB]')
    ax_mag.set_title(titulo)
    ax_mag.legend(loc='best')

    ax_phase.axhline(-180, color='k', linestyle=':', linewidth=0.8)
    ax_phase.set_ylabel('Fase [°]')
    ax_phase.set_xlabel('Frecuencia [rad/s]')
    ax_phase.legend(loc='best')

    plt.tight_layout()
    return fig


# =============================================================================
# Parte B.d - Nyquist compensado
# =============================================================================
def graficar_nyquist_comparativo(G, compensadores, labels, colores=None,
                                  titulo="Nyquist comparativo"):
    """
    Grafica Nyquist del sistema compensado para múltiples controladores.
    """
    if colores is None:
        colores = ['b', 'r', 'g', 'm']

    fig, ax = plt.subplots(figsize=(8, 8))

    # Usar solo frecuencias > 0.05 rad/s para evitar el blow-up del integrador
    omega = np.logspace(-1.3, 2, 3000)

    for C, label, color in zip(compensadores, labels, colores):
        CG = C * G
        resp = ctrl.frequency_response(CG, omega)
        H = resp.frdata[0, 0, :]

        ax.plot(np.real(H), np.imag(H), color=color, linewidth=1.5, label=label)
        ax.plot(np.real(H), -np.imag(H), color=color, linewidth=0.8, alpha=0.5)

    ax.plot(-1, 0, 'rx', markersize=12, markeredgewidth=2, label='(-1, 0)')
    ax.set_xlabel('Re')
    ax.set_ylabel('Im')
    ax.set_title(titulo)
    ax.legend(loc='best')
    # Fija los límites para mantener el punto crítico visible
    ax.set_xlim([-3, 1])
    ax.set_ylim([-3, 3])
    ax.set_aspect('equal')
    plt.tight_layout()
    return fig
