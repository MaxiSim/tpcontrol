"""
generar_figuras.py - Genera todas las figuras del informe en formato PDF/PNG

Ejecutar desde la carpeta src/:
    python3 generar_figuras.py

Las figuras se guardan en ../informe/figuras/
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # backend sin GUI para generar archivos
import matplotlib.pyplot as plt
import control as ctrl

from planta import (crear_planta, crear_planta_ss, crear_planta_perturbacion,
                    info_sistema, L, R, Ke, J, B, Kt)
from analisis import (calcular_margenes, simular_escalon_con_perturbacion,
                      simular_rampa_con_perturbacion, analizar_error_estacionario,
                      calcular_metricas_temporales)
from controlador import (SPECS, ZETA_MIN, SIGMA_MIN,
                         graficar_root_locus_con_region, evaluar_specs,
                         comparar_bode, graficar_respuesta_rampa_compensada)

# Configuración
plt.rcParams.update({
    'figure.figsize': (10, 6),
    'font.size': 11,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'lines.linewidth': 1.5,
    'savefig.dpi': 150,
    'savefig.bbox': 'tight',
})

FIG_DIR = os.path.join(os.path.dirname(__file__), '..', 'informe', 'figuras')
os.makedirs(FIG_DIR, exist_ok=True)

def guardar(fig, nombre):
    path = os.path.join(FIG_DIR, nombre)
    fig.savefig(path)
    plt.close(fig)
    print(f"  ✓ {nombre}")


def main():
    print("Generando figuras del informe...\n")

    G = crear_planta()

    # =========================================================================
    # PARTE A
    # =========================================================================

    # --- A.b) Bode con márgenes ---
    omega = np.logspace(-2, 3, 1000)
    resp = ctrl.frequency_response(G, omega)
    mag_dB = 20 * np.log10(np.abs(resp.fresp[0, 0, :]))
    phase_deg = np.degrees(np.unwrap(np.angle(resp.fresp[0, 0, :])))
    omega_out = resp.omega

    gm, pm, wcg, wcp = ctrl.margin(G)
    gm_dB = 20 * np.log10(gm) if np.isfinite(gm) else float('inf')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7))
    ax1.semilogx(omega_out, mag_dB, 'b-')
    ax1.axhline(0, color='k', linestyle='--', linewidth=0.8)
    if np.isfinite(wcg) and wcg > 0:
        ax1.axvline(wcg, color='m', linestyle=':', alpha=0.7, label=f'$\\omega_{{cg}}$ = {wcg:.2f} rad/s')
        idx = np.argmin(np.abs(omega_out - wcg))
        ax1.annotate(f'GM = {gm_dB:.1f} dB', xy=(wcg, mag_dB[idx]),
                     xytext=(wcg*3, mag_dB[idx]+10),
                     arrowprops=dict(arrowstyle='->', color='purple'), color='purple')
    if np.isfinite(wcp) and wcp > 0:
        ax1.axvline(wcp, color='g', linestyle=':', alpha=0.7, label=f'$\\omega_{{cp}}$ = {wcp:.2f} rad/s')
    ax1.set_ylabel('Magnitud [dB]')
    ax1.set_title('Diagrama de Bode — G(s) a lazo abierto')
    ax1.legend()

    ax2.semilogx(omega_out, phase_deg, 'r-')
    ax2.axhline(-180, color='k', linestyle='--', linewidth=0.8)
    if np.isfinite(wcp) and wcp > 0:
        ax2.axvline(wcp, color='g', linestyle=':', alpha=0.7)
        idx = np.argmin(np.abs(omega_out - wcp))
        ax2.annotate(f'PM = {pm:.1f}°', xy=(wcp, phase_deg[idx]),
                     xytext=(wcp*3, phase_deg[idx]+15),
                     arrowprops=dict(arrowstyle='->', color='green'), color='green')
    if np.isfinite(wcg) and wcg > 0:
        ax2.axvline(wcg, color='m', linestyle=':', alpha=0.7)
    ax2.set_ylabel('Fase [°]')
    ax2.set_xlabel('Frecuencia [rad/s]')
    plt.tight_layout()
    guardar(fig, 'bode_planta.pdf')

    # --- A.c) Escalón + perturbación ---
    t, y_total, y_ref, y_pert = simular_escalon_con_perturbacion(
        G, C=None, amp_ref=1.0, amp_pert=0.1, t_final=30
    )
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t, np.ones_like(t), 'k--', linewidth=1, label='Referencia')
    ax.plot(t, y_total, 'b-', linewidth=2, label='Salida total $\\theta(t)$')
    ax.plot(t, y_ref, 'g--', linewidth=1, alpha=0.7, label='Componente referencia')
    ax.plot(t, y_pert, 'r--', linewidth=1, alpha=0.7, label='Componente perturbación')
    ax.set_xlabel('Tiempo [s]')
    ax.set_ylabel('Posición angular $\\theta$ [rad]')
    ax.set_title('Respuesta ante escalón unitario + perturbación de carga (0.1)')
    ax.legend()
    ax.set_xlim([0, 30])
    plt.tight_layout()
    guardar(fig, 'escalon_perturbacion.pdf')

    # --- A.d) Rampa + perturbación ---
    t, y_total, y_rampa, y_pert, ref = simular_rampa_con_perturbacion(
        G, C=None, amp_pert=0.1, t_final=30
    )
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), height_ratios=[3, 1])
    ax1.plot(t, ref, 'k--', linewidth=1, label='Referencia (rampa)')
    ax1.plot(t, y_total, 'b-', linewidth=2, label='Salida total $\\theta(t)$')
    ax1.plot(t, y_rampa, 'g--', linewidth=1, alpha=0.7, label='Componente rampa')
    ax1.plot(t, y_pert, 'r--', linewidth=1, alpha=0.7, label='Componente perturbación')
    ax1.set_ylabel('$\\theta$ [rad]')
    ax1.set_title('Respuesta ante rampa unitaria + perturbación de carga (0.1)')
    ax1.legend()
    ax1.set_xlim([0, 30])
    error = ref - y_total
    ax2.plot(t, error, 'm-', linewidth=1.5)
    ax2.set_xlabel('Tiempo [s]')
    ax2.set_ylabel('Error [rad]')
    ax2.set_title('Error de seguimiento $e(t) = r(t) - \\theta(t)$')
    ax2.set_xlim([0, 30])
    plt.tight_layout()
    guardar(fig, 'rampa_perturbacion.pdf')

    # =========================================================================
    # PARTE B - Diseño del controlador
    # =========================================================================

    # --- B.a) Respuesta sin compensar con bandas ---
    T_sin = ctrl.feedback(G, 1)
    t_step = np.linspace(0, 30, 3000)
    t_out, y_out = ctrl.step_response(T_sin, t_step)
    y_out = np.array(y_out).flatten()

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t_out, y_out, 'b-', linewidth=2, label='Salida $\\theta(t)$')
    ax.axhline(1.0, color='k', linestyle='--', linewidth=1, label='Referencia')
    ax.axhline(1.02, color='r', linestyle=':', alpha=0.5)
    ax.axhline(0.98, color='r', linestyle=':', alpha=0.5, label='Banda $\\pm 2\\%$')
    ax.axvline(SPECS['ts_max'], color='g', linestyle='--', alpha=0.5,
               label=f'$t_{{s,max}}$ = {SPECS["ts_max"]} s')
    ax.set_xlabel('Tiempo [s]')
    ax.set_ylabel('$\\theta$ [rad]')
    ax.set_title('Respuesta al escalón — Sin compensar ($C(s) = 1$)')
    ax.legend()
    ax.set_xlim([0, 30])
    plt.tight_layout()
    guardar(fig, 'escalon_sin_compensar.pdf')

    # --- B.b) Root locus sin compensar ---
    fig = graficar_root_locus_con_region(
        G, C=None, titulo='Root Locus — $G(s)$ sin compensar'
    )
    guardar(fig, 'rlocus_sin_compensar.pdf')

    # --- B.b) Diseño del compensador (reproducir la lógica del notebook) ---
    sigma_d = 0.7
    zeta_d = 0.6
    wd_d = sigma_d * np.sqrt(1 - zeta_d**2) / zeta_d
    s_d = complex(-sigma_d, wd_d)

    G_eval = ctrl.evalfr(G, s_d)
    angulo_G = np.degrees(np.angle(G_eval))
    fase_necesaria = -180 - angulo_G
    while fase_necesaria > 0: fase_necesaria -= 360
    while fase_necesaria < -360: fase_necesaria += 360

    zc = -2.0
    angulo_cero = np.degrees(np.angle(s_d - zc))
    angulo_polo_nec = angulo_cero - fase_necesaria
    pc = np.real(s_d) - np.imag(s_d) / np.tan(np.radians(angulo_polo_nec))

    C_lead = ctrl.tf([1, -zc], [1, -pc])
    CG_eval = ctrl.evalfr(C_lead * G, s_d)
    K = 1.0 / np.abs(CG_eval)
    C_final = K * C_lead

    # --- B.b) Root locus compensado ---
    fig = graficar_root_locus_con_region(
        G, C=C_final, titulo='Root Locus — Sistema compensado $C(s) \\cdot G(s)$'
    )
    T_comp = ctrl.feedback(C_final * G, 1)
    polos_lc = ctrl.poles(T_comp)
    ax = fig.axes[0]
    for p in polos_lc:
        ax.plot(np.real(p), np.imag(p), 'k*', markersize=15, markeredgewidth=1.5)
    guardar(fig, 'rlocus_compensado.pdf')

    # --- B.b) Escalón compensado vs sin compensar + perturbación ---
    t_c, y_total_c, y_ref_c, y_pert_c = simular_escalon_con_perturbacion(
        G, C=C_final, amp_ref=1.0, amp_pert=0.2, t_final=30
    )
    t_sc, y_total_sc, _, _ = simular_escalon_con_perturbacion(
        G, C=None, amp_ref=1.0, amp_pert=0.2, t_final=30
    )

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t_sc, y_total_sc, 'b--', linewidth=1, alpha=0.6, label='Sin compensar ($C=1$)')
    ax.plot(t_c, y_total_c, 'r-', linewidth=2, label='Compensado')
    ax.axhline(1.0, color='k', linestyle='--', linewidth=0.8, label='Referencia')
    ax.axhline(1.02, color='gray', linestyle=':', alpha=0.5)
    ax.axhline(0.98, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Tiempo [s]')
    ax.set_ylabel('$\\theta$ [rad]')
    ax.set_title('Respuesta al escalón + perturbación (0.2): Compensado vs Sin Compensar')
    ax.legend()
    ax.set_xlim([0, 30])
    plt.tight_layout()
    guardar(fig, 'escalon_comp_vs_sin.pdf')

    # --- B.c) Bode comparativo ---
    fig = comparar_bode(G, C_final,
                        titulo='Diagrama de Bode — Compensado vs Sin Compensar')
    guardar(fig, 'bode_comparativo.pdf')

    # --- B.d) Rampa compensada ---
    t_ramp = np.linspace(0, 30, 3000)
    ref_ramp = t_ramp
    T_comp_r = ctrl.feedback(C_final * G, 1)
    _, y_rampa_comp = ctrl.forced_response(T_comp_r, t_ramp, ref_ramp)
    y_rampa_comp = np.array(y_rampa_comp).flatten()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), height_ratios=[3, 1])
    ax1.plot(t_ramp, ref_ramp, 'k--', linewidth=1, label='Referencia (rampa)')
    ax1.plot(t_ramp, y_rampa_comp, 'r-', linewidth=2, label='Salida compensada')
    ax1.set_ylabel('$\\theta$ [rad]')
    ax1.set_title('Seguimiento de rampa — Compensador lead (sin integrador)')
    ax1.legend()
    ax1.set_xlim([0, 30])
    error_r = ref_ramp - y_rampa_comp
    ax2.plot(t_ramp, error_r, 'm-', linewidth=1.5)
    ax2.set_xlabel('Tiempo [s]')
    ax2.set_ylabel('Error [rad]')
    ax2.set_xlim([0, 30])
    plt.tight_layout()
    guardar(fig, 'rampa_lead.pdf')

    # --- B.e) Con integrador ---
    # Verificar si el integrador simple funciona o hay que rediseñar
    C_int_simple = C_final * ctrl.tf([1], [1, 0])
    T_int_test = ctrl.feedback(C_int_simple * G, 1)
    polos_test = ctrl.poles(T_int_test)
    estable = all(np.real(polos_test) < 0)

    gm_t, pm_t, _, _ = ctrl.margin(C_int_simple * G)

    if not estable or pm_t < 30:
        # Rediseño
        G_aug = G * ctrl.tf([1], [1, 0])
        G_aug_eval = ctrl.evalfr(G_aug, s_d)
        angulo_G_aug = np.degrees(np.angle(G_aug_eval))
        fase_nec_2 = -180 - angulo_G_aug
        while fase_nec_2 > 0: fase_nec_2 -= 360
        while fase_nec_2 < -360: fase_nec_2 += 360

        # Cero más cercano al origen para no requerir polo en SPD.
        # Con zc=-0.5 la condición angular se satisface con pc<0 (polo en LHP).
        zc2 = -0.5
        angulo_cero2 = np.degrees(np.angle(s_d - zc2))
        angulo_polo_nec2 = angulo_cero2 - fase_nec_2
        pc2 = np.real(s_d) - np.imag(s_d) / np.tan(np.radians(angulo_polo_nec2))

        C_lead2 = ctrl.tf([1, -zc2], [1, -pc2])
        K2 = 1.0 / np.abs(ctrl.evalfr(C_lead2 * G_aug, s_d))
        C_con_int = K2 * C_lead2 * ctrl.tf([1], [1, 0])
    else:
        C_con_int = C_int_simple

    fig = graficar_respuesta_rampa_compensada(
        G, C_final, C_con_int, t_final=40,
        titulo='Seguimiento de rampa — Con y sin integrador'
    )
    guardar(fig, 'rampa_con_integrador.pdf')

    print(f"\n✓ Todas las figuras generadas en {os.path.abspath(FIG_DIR)}")

    # Devolver controladores para uso externo
    return C_final, C_con_int


if __name__ == '__main__':
    main()
