"""
planta2.py - Modelo simplificado del rotor (Ejercicio 2)

Planta: G(s) = K / [s(Js + B)]  con K=1, J=1.
Se desprecia la dinámica eléctrica (constante de tiempo mucho menor que la mecánica).

El parámetro B puede degradarse de 0.5 (nominal) a -0.1 (falla).

Curso: IOR442 Sistemas de Control - Dr. Mariano Scaramal
"""

import control as ctrl
import numpy as np

# Parámetros simplificados
K = 1       # ganancia (simplificada)
J = 1       # momento de inercia [normalizado]
B_NOM = 0.5   # fricción nominal
B_FALLA = -0.1  # fricción con falla (amortiguamiento negativo → inestable)


def crear_planta2(B_val=B_NOM):
    """
    Crea G(s) = K / [s(Js + B)] para un valor dado de B.

    Con K=1, J=1:
        G(s) = 1 / [s(s + B)]

    - B > 0: dos polos reales en s=0 y s=-B (estable en LA, marginalmente)
    - B < 0: polo en s=|B| (semiplano derecho → inestable en LA)

    Args:
        B_val: valor del coeficiente de fricción

    Returns:
        G: función de transferencia
    """
    num = [K]
    den = [J, B_val, 0]  # Js² + Bs + 0 → s(Js + B)
    G = ctrl.tf(num, den)
    return G


def crear_planta2_ss(B_val=B_NOM):
    """
    Representación en espacio de estados para simulación LTV.

    Estados: x = [θ, ω]ᵀ
        dθ/dt = ω
        dω/dt = -(B/J)·ω + (K/J)·u

    Salida: y = θ

    Args:
        B_val: coeficiente de fricción

    Returns:
        A, B_mat, C, D: matrices del espacio de estados (numpy arrays)
    """
    A = np.array([
        [0,      1],
        [0, -B_val/J]
    ])
    B_mat = np.array([
        [0],
        [K/J]
    ])
    C = np.array([[1, 0]])
    D = np.array([[0]])
    return A, B_mat, C, D


def info_planta2(B_val=B_NOM):
    """Imprime resumen del sistema simplificado."""
    G = crear_planta2(B_val)
    polos = ctrl.poles(G)

    print(f"{'='*50}")
    print(f"  Planta simplificada: G(s) = 1/[s(s + {B_val})]")
    print(f"  Polos: {polos}")
    estable_la = all(np.real(polos) <= 0)
    print(f"  Estable en LA: {'Sí (marginalmente)' if estable_la else 'NO'}")
    print(f"{'='*50}")
    return G, polos
