"""
planta.py - Modelo del motor DC para control de rotor de aspas (Ejercicio 1)

Define los parámetros físicos del motor, la función de transferencia G(s) = Θ(s)/V(s)
y la representación en variables de estado.

Curso: IOR442 Sistemas de Control - Dr. Mariano Scaramal
"""

import control as ctrl
import numpy as np

# =============================================================================
# Parámetros del motor DC
# =============================================================================
L = 0.15      # Inductancia de armadura [H]
R = 0.7       # Resistencia de armadura [Ω]
Ke = 0.1      # Constante de fuerza electromotriz [V·s/rad]
J = 0.02      # Momento de inercia total [kg·m²]
B = 0.25      # Coeficiente de fricción viscosa [N·m·s/rad]
Kt = 0.25     # Constante de torque [N·m/A]

# =============================================================================
# Coeficientes del denominador de G(s)
# El denominador proviene de expandir s·[(Ls+R)(Js+B) + Ke·Kt]
#   (Ls+R)(Js+B) = LJs² + (LB+RJ)s + RB
#   Sumando Ke·Kt: LJs² + (LB+RJ)s + (RB + Ke·Kt)
# =============================================================================
a2 = L * J                    # 0.003  - coef. de s²
a1 = L * B + R * J            # 0.0515 - coef. de s
a0 = R * B + Ke * Kt          # 0.2    - término independiente


def crear_planta():
    """
    Construye la función de transferencia del motor DC:
        G(s) = Θ(s)/V(s) = Kt / [s · (LJs² + (LB+RJ)s + (RB+KeKt))]

    Es un sistema de 3er orden con un integrador puro (Tipo 1).

    Returns:
        G: función de transferencia (control.TransferFunction)
    """
    num = [Kt]
    den = [a2, a1, a0, 0]  # el 0 final corresponde al factor 's' (integrador)
    G = ctrl.tf(num, den)
    return G


def crear_planta_ss():
    """
    Construye la representación en variables de estado del motor DC.

    Estados: x = [θ, ω, i]ᵀ
        - θ: posición angular [rad]
        - ω: velocidad angular [rad/s]
        - i: corriente de armadura [A]

    Ecuaciones:
        dθ/dt = ω
        dω/dt = -(B/J)·ω + (Kt/J)·i
        di/dt = -(Ke/L)·ω - (R/L)·i + (1/L)·V

    Salida: y = θ

    Returns:
        sys_ss: sistema en espacio de estados (control.StateSpace)
    """
    A = np.array([
        [0,    1,       0      ],
        [0,   -B/J,     Kt/J   ],
        [0,   -Ke/L,   -R/L    ]
    ])

    B_mat = np.array([
        [0],
        [0],
        [1/L]
    ])

    C = np.array([[1, 0, 0]])
    D = np.array([[0]])

    sys_ss = ctrl.ss(A, B_mat, C, D)
    return sys_ss


def crear_planta_perturbacion():
    """
    Construye la función de transferencia de perturbación:
        Θ(s)/N(s) con realimentación unitaria

    La perturbación de carga (torque) entra en la ecuación mecánica:
        J·dω/dt + B·ω = Kt·i + n(t)

    En lazo abierto, la transferencia de la perturbación al ángulo es:
        G_n(s) = 1 / [s · (Js + B)]

    En lazo cerrado con realimentación unitaria y controlador C(s):
        Θ(s)/N(s) = G_n(s) / [1 + C(s)·G(s)]

    Para C(s) = 1 (sin compensador):
        Θ(s)/N(s) = G_n(s) / [1 + G(s)]

    Returns:
        G_n_ol: transferencia perturbación→ángulo en lazo abierto
    """
    # Perturbación → ω: 1/(Js + B), luego ω → θ: 1/s
    num_n = [1]
    den_n = [J, B, 0]  # s·(Js + B) = Js² + Bs
    G_n_ol = ctrl.tf(num_n, den_n)
    return G_n_ol


def info_sistema():
    """Imprime un resumen de los parámetros y polos del sistema."""
    G = crear_planta()
    polos = ctrl.poles(G)
    ceros = ctrl.zeros(G)

    print("=" * 60)
    print("  Motor DC - Parámetros del Sistema")
    print("=" * 60)
    print(f"  L  = {L} H      (inductancia)")
    print(f"  R  = {R} Ω      (resistencia)")
    print(f"  Ke = {Ke} V·s/rad (cte. FEM)")
    print(f"  J  = {J} kg·m²  (inercia)")
    print(f"  B  = {B} N·m·s/rad (fricción)")
    print(f"  Kt = {Kt} N·m/A   (cte. torque)")
    print("-" * 60)
    print(f"  G(s) = {Kt} / [s · ({a2}s² + {a1}s + {a0})]")
    print(f"  Polos: {polos}")
    print(f"  Ceros: {ceros}")
    print(f"  Tipo del sistema: 1 (un integrador)")
    print("=" * 60)
    return G, polos, ceros
