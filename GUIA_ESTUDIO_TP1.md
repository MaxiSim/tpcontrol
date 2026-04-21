# Guía de Estudio — TP1 Sistemas de Control

**Para qué sirve este documento**: para que entiendan qué hicimos, por qué, y qué les van a preguntar en la defensa. Está pensado para ingenieros de último año que necesitan entender rápido los fundamentos sin leer 300 páginas de Ogata.

---

## Índice

1. [Conceptos fundamentales que necesitan sí o sí](#1-conceptos-fundamentales)
2. [Ejercicio 1: qué pide, qué hicimos, por qué](#2-ejercicio-1)
3. [Ejercicio 2: qué pide, qué hicimos, por qué](#3-ejercicio-2)
4. [Preguntas típicas de defensa](#4-preguntas-típicas-de-defensa)

---

## 1. Conceptos fundamentales

### 1.1 Función de transferencia

Es la relación entre la salida y la entrada de un sistema **en el dominio de Laplace** (con condiciones iniciales nulas). Si tenés una ecuación diferencial, la transformás con Laplace y despejás:

```
G(s) = Salida(s) / Entrada(s)
```

**Ejemplo concreto del TP**: el motor DC tiene 3 ecuaciones diferenciales (eléctrica, mecánica, cinemática). Las transformamos con Laplace, sustituimos una en otra, y llegamos a:

```
G(s) = Θ(s)/V(s) = Kt / [s · (LJs² + (LB+RJ)s + RB+KeKt)]
```

Esto dice: "si le meto tensión V al motor, la posición angular θ que obtengo es V multiplicado por G(s)".

### 1.2 Polos y ceros

- **Polos**: raíces del denominador de G(s). Determinan la **estabilidad** y la **forma** de la respuesta.
- **Ceros**: raíces del numerador. Afectan la **amplitud** de la respuesta.

**Regla de oro**: si algún polo tiene parte real positiva (está a la derecha del eje imaginario en el plano s), el sistema es **inestable** → la salida diverge al infinito.

```
Plano s:

     Im
      |
      |   ← ESTABLE    INESTABLE →
      |
------+------→ Re
      |
      |
```

### 1.3 Tipo del sistema

El "tipo" es la **cantidad de integradores puros** (polos en s=0) en la función de transferencia de lazo abierto. Determina qué entradas puede seguir con error nulo:

| Tipo | Polos en s=0 | Escalón | Rampa | Parábola |
|------|:---:|:---:|:---:|:---:|
| 0 | 0 | Error finito | Error → ∞ | Error → ∞ |
| **1** | **1** | **Error = 0** | **Error finito** | Error → ∞ |
| 2 | 2 | Error = 0 | Error = 0 | Error finito |

**En el TP**: nuestra planta del Ej.1 es Tipo 1 (tiene un `s` suelto en el denominador). Por eso sigue el escalón perfecto pero no la rampa. Para seguir la rampa, le agregamos un integrador (1/s) → pasa a Tipo 2.

### 1.4 Lazo cerrado con realimentación unitaria

Este es el diagrama que aparece en todo el TP:

```
          error                    salida
R(s) ──►(+)──► C(s) ──► G(s) ──► Y(s)
          -▲                        │
           │                        │
           └────────────────────────┘
```

- **R(s)**: referencia (lo que queremos)
- **Y(s)**: salida (lo que obtenemos)
- **C(s)**: controlador (lo que diseñamos)
- **G(s)**: planta (el motor, lo que nos dan)
- **Error**: e(s) = R(s) - Y(s)

La **transferencia de lazo cerrado** es:

```
T(s) = C(s)·G(s) / [1 + C(s)·G(s)]
```

Si C(s) = 1 (sin controlador), queda T(s) = G(s) / [1 + G(s)].

### 1.5 Márgenes de estabilidad (GM y PM)

Son medidas de **cuánto aguanta** el sistema antes de volverse inestable.

**Margen de Ganancia (GM)**: cuánto podemos **multiplicar** la ganancia del sistema antes de que se vuelva inestable. Se mide en la frecuencia donde la fase cruza -180°.
- GM > 0 dB → estable
- GM > 6 dB → aceptable
- Referencia: "¿Cuántas veces puedo amplificar la señal antes de que todo explote?"

**Margen de Fase (PM)**: cuánta **fase adicional** puede perder el sistema antes de ser inestable. Se mide en la frecuencia de cruce de ganancia (donde |G| = 1 = 0 dB).
- PM > 0° → estable
- PM > 45° → buen diseño
- Referencia: "¿Cuánto retraso puedo tener antes de que todo explote?"

**Diagrama de Bode**: son dos gráficos (magnitud y fase vs frecuencia) donde se leen estos márgenes. El GM se lee en el gráfico de magnitud y el PM en el de fase.

### 1.6 Root Locus (Lugar de las raíces)

Es un gráfico que muestra **cómo se mueven los polos de lazo cerrado** cuando variamos la ganancia K de 0 a ∞.

- Empieza en los **polos** de lazo abierto (cuando K=0)
- Termina en los **ceros** de lazo abierto (cuando K→∞) o se va al infinito
- Si alguna rama cruza al semiplano derecho → el sistema se vuelve inestable para esa ganancia

**Para qué nos sirve**: para elegir dónde queremos que estén los polos de lazo cerrado. Si queremos respuesta rápida y amortiguada, necesitamos polos en una zona específica del plano s.

### 1.7 Diagrama de Nyquist

Es el gráfico de G(jω) en el plano complejo cuando ω varía de 0 a ∞ (y su reflejo para -∞ a 0).

**Criterio de Nyquist**: Z = N + P
- P = polos de G(s) en el semiplano derecho
- N = encirculamientos horarios del punto (-1, 0)
- Z = polos de lazo cerrado en el semiplano derecho

**Para estabilidad**: necesitamos Z = 0.

### 1.8 Compensador Lead (adelanto de fase)

```
C(s) = Kc · (s + z) / (s + p),    donde p > z
```

**Qué hace**: aporta fase positiva en un rango de frecuencias.

**Analogía**: es como adelantar el timing en un motor de combustión. Le estás diciendo al controlador "reaccioná un poco antes".

**Cuándo se usa**: cuando necesitás mejorar el **transitorio** (reducir sobrepico, tiempo de establecimiento) o aumentar el **margen de fase**.

**Cuándo NO se usa**: cuando necesitás mejorar solo la precisión en régimen permanente (para eso se usa un lag o un integrador).

### 1.9 Perturbaciones y superposición

Si el sistema tiene dos entradas (referencia + perturbación), por **superposición** (porque es lineal):

```
Salida = T_ref(s) · R(s) + T_pert(s) · N(s)
```

Cada componente se analiza por separado y después se suman. Esto es exactamente lo que hacemos en los puntos A.c y A.d del Ejercicio 1.

---

## 2. Ejercicio 1: Sistema de Control de un Rotor de Aspas

### Qué pide la consigna (en criollo)

Tenés un motor DC que mueve las aspas de un VTOL. Te dan las ecuaciones y los parámetros. Tenés que:

**Parte A — Entender el sistema:**
1. Sacar la función de transferencia y las ecuaciones de estado
2. Ver los márgenes de estabilidad (¿el sistema se banca realimentación directa?)
3. Simular con escalón + perturbación de torque
4. Simular con rampa + perturbación
5. Proponer qué controlador conviene

**Parte B — Diseñar un controlador:**
1. Verificar que sin controlador no se cumplen las specs (ts, Mp, error)
2. Diseñar un lead por root locus
3. Comparar Bode con y sin compensador
4. Ver si sigue una rampa
5. Si no → agregar integrador

### Qué hicimos y por qué

#### Parte A.a — Función de transferencia

Tomamos las 3 ecuaciones diferenciales del motor:
- **Eléctrica**: L·di/dt + R·i + Ke·ω = V (ley de Kirchhoff)
- **Mecánica**: J·dω/dt + B·ω = Kt·i (Newton rotacional)
- **Cinemática**: ω = dθ/dt

Las transformamos con Laplace (derivada → multiplicar por s), despejamos y llegamos a:

```
G(s) = 0.25 / [s · (0.003s² + 0.0515s + 0.2)]
```

**Resultado clave**: sistema de 3er orden, Tipo 1, con polos en s=0, s≈-5.94, s≈-11.23 (todos reales, no complejos).

**Variables de estado**: elegimos x = [θ, ω, i] porque son las variables físicas del sistema. Las matrices A, B, C, D salen directo de las ecuaciones.

#### Parte A.b — Márgenes

Graficamos el Bode de G(s) y leemos:
- GM ≈ 22.8 dB (muy bueno, >6 dB)
- PM ≈ 72.2° (excelente, >45°)

**Conclusión**: el sistema sin compensar ya es estable en lazo cerrado. Los márgenes son holgados. Pero esto no significa que cumple las specs de transitorio (ts, Mp).

#### Parte A.c — Escalón + perturbación

La perturbación de carga entra como torque en la ecuación mecánica. En lazo cerrado:

```
θ(s) = [G/(1+G)] · R(s) + [Gn/(1+G)] · N(s)
```

donde Gn(s) = 1/[s(Js+B)] es la transferencia de perturbación.

**Resultado**: el sistema sigue el escalón con error nulo (Tipo 1 → Kp = ∞) pero la perturbación agrega un pequeño offset.

#### Parte A.d — Rampa + perturbación

Con rampa r(t) = t, el error es:

```
e_ss = 1/Kv    (Kv = constante de velocidad, finita para Tipo 1)
```

**Resultado**: el sistema NO puede seguir la rampa con error nulo. Siempre queda retrasado.

#### Parte A.e — Propuesta de controlador

| Caso | Controlador | Por qué |
|------|-------------|---------|
| Escalón + pert. | **Lead** | Mejora transitorio (ts, Mp) y agrega ganancia |
| Rampa + pert. | **Lead + integrador** | El integrador da Tipo 2 → error de rampa = 0 |

**No usamos lag** porque empeoraría el transitorio (reduce ancho de banda).

#### Parte B.a — Sistema sin compensar vs specs

Specs:
- ts ≤ 7.5s (2%)
- Mp ≤ 15%
- Error por perturbación (0.2) < 5%

Estas specs se traducen a restricciones en el plano s:
- De Mp ≤ 15% → ζ ≥ 0.517 (líneas de amortiguamiento)
- De ts ≤ 7.5s → σ ≥ 0.533 (línea vertical en el plano s)

Los polos de lazo cerrado deben estar en la intersección de ambas regiones.

#### Parte B.b — Diseño del lead por root locus

**Método de contribución angular**:

1. **Elegimos un polo deseado** sd = -0.7 + j·0.933 (dentro de la región admisible, con margen)
2. **Calculamos** ∠G(sd): cuánta fase aporta la planta sola en ese punto
3. **La diferencia** con -180° es lo que debe aportar el compensador
4. **Colocamos el cero** del lead en -2 (para atraer el root locus)
5. **Calculamos el polo** para que la contribución angular sea la correcta
6. **Ajustamos K** con la condición de magnitud: |C(sd)·G(sd)| = 1

**Por qué funciona**: el root locus muestra las trayectorias de los polos. Si agregamos un cero y un polo (el lead), cambiamos la forma del root locus para que pase por donde queremos.

#### Parte B.c — Comparación Bode

Graficamos el Bode con y sin compensador. El lead:
- Aumenta la fase cerca de la frecuencia de cruce → PM sube
- Modifica el ancho de banda

#### Parte B.d — Rampa

Con el lead solo, el sistema sigue siendo Tipo 1 → error finito ante rampa. El lead mejoró el transitorio pero no la precisión ante rampa.

#### Parte B.e — Agregar integrador

```
C_nuevo(s) = C_lead(s) · (1/s)
```

Esto agrega un polo en s=0 → Tipo 2 → error de rampa = 0.

**Riesgo**: el 1/s agrega -90° de fase en todas las frecuencias → reduce PM. Por eso probablemente hay que rediseñar el lead para compensar.

---

## 3. Ejercicio 2: Análisis y Control frente a Perturbaciones

### Qué pide la consigna (en criollo)

El mismo rotor pero simplificado: G(s) = 1/[s(s+B)]. El coeficiente de fricción B se puede degradar de 0.5 a -0.1 (falla física real por calentamiento). Tenés que:

**Parte A — Entender la falla:**
1. Nyquist del sistema nominal → ¿es estable?
2. Simular qué pasa cuando B cambia de 0.5 a -0.1 en t=15s
3. Root locus para ambos casos → ¿por qué se vuelve inestable?
4. Proponer un compensador

**Parte B — Diseñar compensadores robustos:**
1. Diseñar leads para PM = 35°, 50°, 65° (con ganancia DC = 1)
2. Verificar en Bode
3. Simular la falla con cada compensador
4. Nyquist compensado + analizar desventajas de mucho PM

### Qué hicimos y por qué

#### Parte A.a — Nyquist nominal

Con B = 0.5, G(s) = 1/[s(s+0.5)]. Polos en s=0 y s=-0.5 (ambos en el semiplano izquierdo o en el borde).

Criterio de Nyquist: P=0 (no hay polos en SPD), N=0 (no rodea -1,0), Z=0 → **estable en LC**.

**Concepto clave**: Nyquist nos dice si el lazo cerrado es estable mirando solo el lazo abierto. No necesitamos cerrar el lazo para saberlo.

#### Parte A.b — Simulación de la falla

Acá hay algo importante: **cuando B cambia en el tiempo, el sistema deja de ser LTI** (lineal tiempo-invariante). No podemos usar funciones de transferencia porque estas asumen parámetros constantes.

**Solución**: usamos espacio de estados con integración numérica (solve_ivp de scipy). Las matrices A cambian cuando t cruza los 15 segundos.

**Resultado**: con B = -0.1, la planta tiene un polo en s = +0.1. La realimentación unitaria no alcanza para estabilizar → la salida diverge.

**Por qué B negativo es catastrófico**: B es fricción. Con B > 0, la fricción frena el sistema (estabilizante). Con B < 0, es como si el sistema se autoalimentara: en vez de frenarse, se acelera solo.

#### Parte A.c — Root locus comparativo

- **B = 0.5**: polos de LA en s=0 y s=-0.5. El root locus se queda en el semiplano izquierdo para K razonable → estable.
- **B = -0.1**: polos de LA en s=0 y s=+0.1. Una rama del root locus sale del SPD. Para K=1, los polos de LC están cerca del SPD o dentro → inestable.

#### Parte A.d — Propuesta: compensador lead

**Lead**, porque:
- Aporta fase positiva → aleja la curva de Nyquist del punto (-1,0)
- Mayor margen de fase = mayor tolerancia a que B se degrade
- Lag empeoraría el transitorio
- Lead-lag es innecesariamente complejo

#### Parte B.a — Diseño para PM = 35°, 50°, 65°

**Método de diseño por Bode** (distinto al Ej.1 que fue por root locus):

1. Calculamos PM actual del sistema sin compensar (≈28°)
2. Fase extra necesaria: φ = PM_deseado - PM_actual + 5° (margen)
3. Parámetro α = (1 - sin φ) / (1 + sin φ) → cuánto "comprime" el lead
4. Frecuencia de máxima fase: donde |G| = √α
5. Cero z = ωm·√α, polo p = ωm/√α
6. Ganancia Kc = p/z → esto normaliza C(0) = 1

**¿Por qué normalizar la ganancia DC?**

```
C(0) = Kc · z/p = (p/z) · (z/p) = 1
```

Si C(0) = 1, el compensador no afecta la ganancia en frecuencia cero → no cambia el error de posición del sistema. Esto es deseable porque solo queremos mejorar la fase, no alterar la precisión que ya tenemos.

**¿Por qué es posible?** Porque un lead siempre tiene |z/p| < 1 (el cero está más cerca del origen que el polo). Al multiplicar por Kc = p/z, compensamos exactamente.

#### Parte B.b — Verificación Bode

Graficamos el Bode de C(s)·G(s) para los 3 compensadores. Verificamos que:
- El PM obtenido ≈ PM deseado
- La ganancia en baja frecuencia es igual para los 3 (normalización DC)
- El ancho de banda (ωcp) disminuye al aumentar PM

#### Parte B.c — Respuesta con falla

Simulamos escalón unitario con B cambiando en t=15s. Resultados:

| PM | Comportamiento post-falla |
|----|--------------------------|
| 35° | Sobrevive pero con oscilaciones |
| 50° | Se recupera con algo de sobrepico |
| 65° | Recuperación suave, casi primer orden |
| Sin comp. | **Diverge** |

**Observación clave**: a mayor PM, el sistema "se parece más a uno de primer orden". Esto es porque mayor margen de fase = mayor amortiguamiento equivalente = polos más alejados del eje imaginario.

#### Parte B.d — Nyquist compensado + desventajas

El compensador rota la curva de Nyquist alejándola de (-1,0). A más PM, más lejos queda.

**Desventajas de mucho PM**:
1. **Más lento**: menor ancho de banda → el sistema tarda más en responder
2. **ts mayor**: el tiempo de establecimiento crece
3. **Más ruido**: el polo del lead se aleja → amplifica ruido en alta frecuencia
4. **Peor rechazo de perturbaciones**: menos ganancia en frecuencias medias

**Bottom line**: hay un trade-off. No se puede tener máxima robustez Y máxima velocidad al mismo tiempo. El ingeniero elige según el peor caso de variación paramétrica esperado.

---

## 4. Preguntas típicas de defensa

### Sobre la teoría

**P: ¿Por qué el sistema Tipo 1 no puede seguir una rampa con error nulo?**
R: Porque la constante de velocidad Kv es finita. El error ante rampa es e_ss = 1/Kv. Para que sea cero, necesitaríamos Kv = ∞, lo cual requiere Tipo 2 (dos integradores).

**P: ¿Qué pasa si le pongo ganancia infinita al controlador?**
R: El sistema se vuelve inestable. El root locus muestra que para ganancia muy alta, los polos cruzan al semiplano derecho.

**P: ¿Cuál es la diferencia entre lead y lag?**
R: El lead aporta fase positiva (mejora transitorio/estabilidad). El lag aporta ganancia en baja frecuencia (mejora precisión) pero reduce ancho de banda. Lead = rápido y estable. Lag = preciso pero lento.

**P: ¿Por qué no usar un PID directamente?**
R: Un PID es lead + lag + integrador. Funciona, pero el enunciado pide diseñar paso a paso por root locus y frecuencia. Además, para entender qué hace cada parte, conviene diseñarlos por separado.

### Sobre el Ejercicio 1

**P: ¿Por qué los polos son reales y no complejos?**
R: Porque el discriminante del polinomio cuadrático es positivo (Δ ≈ 2.52×10⁻⁴ > 0). Si fuera negativo, tendríamos polos complejos conjugados (oscilaciones naturales).

**P: ¿Cómo elegiste la ubicación del polo deseado?**
R: A partir de las specs. Mp ≤ 15% da ζ ≥ 0.517 (cono de amortiguamiento). ts ≤ 7.5s da σ ≥ 0.533 (línea vertical). El polo deseado debe estar en la intersección, con algo de margen.

**P: ¿Por qué el cero del lead se pone en -2?**
R: Para atraer el root locus hacia la región admisible. El cero "tira" las ramas del root locus hacia él. Lo ponemos entre el origen y el polo dominante de la planta (-5.94).

**P: ¿Qué pasa si el integrador desestabiliza el sistema?**
R: Hay que rediseñar el lead con más fase para compensar los -90° que agrega el integrador. En nuestro código, esto se hace automáticamente.

### Sobre el Ejercicio 2

**P: ¿Por qué B negativo hace al sistema inestable?**
R: Con B > 0, la fricción disipa energía (estabilizante). Con B < 0, el término -B·ω en la ecuación mecánica se convierte en un término que inyecta energía → retroalimentación positiva → inestabilidad. En la función de transferencia, el polo pasa de s = -B (izquierda) a s = |B| (derecha).

**P: ¿Por qué no podés usar funciones de transferencia para simular la falla?**
R: Porque la función de transferencia asume que los parámetros son constantes (LTI). Cuando B cambia en t=15s, el sistema es LTV (variante en el tiempo). Hay que usar espacio de estados con integración numérica paso a paso.

**P: ¿Qué significa "normalizar la ganancia DC"?**
R: Que C(0) = 1. A frecuencia cero (DC = señales constantes), el compensador no amplifica ni atenúa. Esto preserva el error de posición del sistema original. Es posible solo con un lead (no con un integrador, que tiene ganancia infinita en DC).

**P: ¿Por qué no poner PM = 90° y listo?**
R: Porque a mayor PM, más lento es el sistema. El ancho de banda cae, el tiempo de establecimiento sube, y el rechazo de perturbaciones empeora. Hay que balancear robustez vs desempeño. PM = 45°-60° suele ser el sweet spot.

**P: ¿Cómo interpreta el diagrama de Nyquist compensado?**
R: El lead rota la curva de Nyquist en sentido antihorario, alejándola del punto (-1,0). Cuanto más lejos esté la curva de (-1,0), más robusta es la estabilidad. Si la curva llega a rodear (-1,0), el sistema es inestable.

---

## Glosario rápido

| Término | Qué es | Unidad |
|---------|--------|--------|
| G(s) | Función de transferencia de la planta | - |
| C(s) | Función de transferencia del controlador | - |
| PM | Margen de fase (Phase Margin) | grados |
| GM | Margen de ganancia (Gain Margin) | dB |
| ts | Tiempo de establecimiento (settling time) | segundos |
| Mp | Sobrepico (overshoot) | % |
| ζ (zeta) | Coeficiente de amortiguamiento | - |
| ωn | Frecuencia natural no amortiguada | rad/s |
| σ | Parte real del polo (= ζ·ωn) | rad/s |
| Kp | Constante de posición (Tipo 0+) | - |
| Kv | Constante de velocidad (Tipo 1+) | 1/s |
| SPD | Semiplano derecho (Re > 0) | - |
| LA | Lazo abierto | - |
| LC | Lazo cerrado | - |
| LTI | Lineal tiempo-invariante | - |
| LTV | Lineal tiempo-variante | - |

---

## Archivos del proyecto y qué hace cada uno

```
src/
├── planta.py        → Define G(s) del motor DC completo (Ej.1)
├── analisis.py      → Funciones de análisis: Bode, márgenes, simulaciones (Ej.1)
├── controlador.py   → Diseño del lead, evaluación de specs, root locus (Ej.1)
├── planta2.py       → Define G(s) simplificada (Ej.2)
├── robustez.py      → Nyquist, simulación LTV, diseño de leads por Bode (Ej.2)
└── notebooks/
    ├── ejercicio1.ipynb  → Todo el Ej.1 ejecutable con explicaciones
    └── ejercicio2.ipynb  → Todo el Ej.2 ejecutable con explicaciones
```

**Tip para la defensa**: ejecuten los notebooks antes. Si les preguntan algo numérico, pueden mostrar la celda correspondiente.
