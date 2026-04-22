"""
generar_figuras_ej2.py - Genera figuras del Ejercicio 2 para el informe.

Ejecutar desde src/:
    python3 generar_figuras_ej2.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import control as ctrl

from planta2 import crear_planta2, B_NOM, B_FALLA
from robustez import (graficar_nyquist, simular_falla, graficar_respuesta_falla,
                      graficar_rlocus_comparativo, disenar_lead_por_margen_fase,
                      comparar_bode_multiples, graficar_nyquist_comparativo)

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
    fig.savefig(os.path.join(FIG_DIR, nombre))
    plt.close(fig)
    print(f"  OK {nombre}")


def main():
    print("Generando figuras Ejercicio 2...\n")

    G_nom = crear_planta2(B_NOM)

    # A.a) Nyquist nominal
    fig, _ = graficar_nyquist(G_nom, titulo='Nyquist — $G(s)$ nominal ($B = 0{,}5$)')
    guardar(fig, 'ej2_nyquist_nominal.pdf')

    # A.b) Simulación falla sin compensar — t_final=100s para ver la divergencia
    t, theta, _ = simular_falla(C_tf=None, t_falla=15, t_final=100)
    fig = graficar_respuesta_falla(t, theta,
        titulo='Respuesta al escalón con falla en $t=15$ s (sin compensar)')
    guardar(fig, 'ej2_falla_sin_comp.pdf')

    # A.c) Root locus comparativo
    fig = graficar_rlocus_comparativo([B_NOM, B_FALLA],
        labels=[f'$B = {B_NOM}$ (nominal)', f'$B = {B_FALLA}$ (falla)'],
        colores=['b', 'r'])
    guardar(fig, 'ej2_rlocus_comparativo.pdf')

    # B.a) Diseño de los 3 compensadores
    PM_objs = [35, 50, 65]
    compensadores = []
    infos = []
    for PM in PM_objs:
        C, info = disenar_lead_por_margen_fase(G_nom, PM, B_val=B_NOM)
        compensadores.append(C)
        infos.append(info)

    # B.b) Bode comparativo
    labels = [f'PM = {i["PM_deseado"]}° (obt: {i["PM_obtenido"]:.1f}°)' for i in infos]
    fig = comparar_bode_multiples(G_nom, compensadores, labels,
        colores=['b', 'r', 'g'],
        titulo='Bode del sistema nominal compensado')
    guardar(fig, 'ej2_bode_3compensadores.pdf')

    # B.c) Respuesta con falla para los 3 compensadores
    fig, ax = plt.subplots(figsize=(10, 7))
    for C, info, color in zip(compensadores, infos, ['b', 'r', 'g']):
        t, theta, _ = simular_falla(C_tf=C, t_falla=15, t_final=50)
        ax.plot(t, theta, color=color, linewidth=1.5,
                label=f'PM = {info["PM_deseado"]}°')
    t0, theta0, _ = simular_falla(C_tf=None, t_falla=15, t_final=50)
    ax.plot(t0, theta0, 'k--', linewidth=1, alpha=0.5, label='Sin compensar')
    ax.axhline(1.0, color='k', linestyle=':', linewidth=0.8, alpha=0.3)
    ax.axvline(15, color='gray', linestyle=':', linewidth=1, alpha=0.5, label='Falla ($t=15$ s)')
    ax.set_xlabel('Tiempo [s]')
    ax.set_ylabel('$\\theta$ [rad]')
    ax.set_title('Respuesta al escalón con falla ($B$: 0,5 → −0,1 en $t=15$ s)')
    ax.legend(loc='best')
    ax.set_xlim([0, 50])
    ax.set_ylim([-2, 5])
    plt.tight_layout()
    guardar(fig, 'ej2_falla_3compensadores.pdf')

    # B.d) Nyquist comparativo
    labels_nyq = [f'PM = {i["PM_deseado"]}°' for i in infos]
    fig = graficar_nyquist_comparativo(G_nom, compensadores, labels_nyq,
        colores=['b', 'r', 'g'],
        titulo='Nyquist del sistema nominal compensado')
    guardar(fig, 'ej2_nyquist_compensado.pdf')

    # =========================================================================
    # Generar valores numéricos para el informe LaTeX
    # =========================================================================
    from analisis import calcular_metricas_temporales as _met

    def fmt(x, d=2):
        v = float(x)
        if not np.isfinite(v):
            return r'\infty'
        return f"{v:.{d}f}".replace('.', '{,}')

    prefijos = ['TreintaCinco', 'Cincuenta', 'SesentaCinco']
    val_path = os.path.join(os.path.dirname(__file__), '..', 'informe', 'valores_ej2.tex')
    with open(val_path, 'w') as vf:
        def w(cmd, val):
            vf.write(rf'\def\{cmd}{{{val}}}' + '\n')
        vf.write('% AUTO-GENERADO por generar_figuras_ej2.py — no editar a mano\n')
        for pref, C, info in zip(prefijos, compensadores, infos):
            _, pm_comp, _, wcp_comp = ctrl.margin(C * G_nom)
            met = _met(G_nom, C, t_final=30)
            w(f'PM{pref}Obt',   fmt(info['PM_obtenido'], 1))
            w(f'PM{pref}Alpha', fmt(info['alpha'], 4))
            w(f'PM{pref}Zcero', fmt(info['z'], 4))
            w(f'PM{pref}Polo',  fmt(info['p'], 4))
            w(f'PM{pref}Kc',    fmt(info['Kc'], 4))
            w(f'PM{pref}Wcp',   fmt(wcp_comp))
            w(f'PM{pref}Mp',    fmt(met['Mp'], 1))
            w(f'PM{pref}Ts',    fmt(met['ts']))
    print("  OK valores_ej2.tex")

    print(f"\nTodas las figuras Ej.2 generadas en {os.path.abspath(FIG_DIR)}")
    return compensadores, infos


if __name__ == '__main__':
    main()
