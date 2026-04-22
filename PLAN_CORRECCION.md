# Plan de corrección — TP1 Sistemas de Control

**Audiencia:** agente Sonnet que ejecutará este plan de forma iterativa y autónoma.
**Objetivo:** que el código y el informe cumplan *todas* las consignas del TP1, en particular las que hoy no se cumplen, y que los valores numéricos del informe reflejen exactamente lo que produce el código.

Antes de empezar: leé este archivo completo. La regla de oro es **iterar**: modificar código → correr script → leer valores → comparar con specs → volver a modificar si algo no cumple. No pasar a la etapa del informe hasta que el código produzca valores que cumplan todas las condiciones numéricas.

---

## 0. Contexto rápido

- Repo: `/Users/maxi/Documents/Udesa/Sistemas de control/tpcontrol/`
- Venv: `./.venv/` (usar `../.venv/bin/python3` desde `src/`)
- Scripts principales: `src/generar_figuras.py` (Ej.1), `src/generar_figuras_ej2.py` (Ej.2)
- Informe LaTeX: `informe/informe.tex`. Valores auto-generados: `informe/valores_ej1.tex`, `informe/valores_ej2.tex`
- Consigna: `TP1_Sistemas_de_Control.pdf` (ya leído, contenido resumido abajo).

### Resumen de la consigna relevante

**Ej.1 Parte B — Especificaciones (simultáneas)**:
1. $t_s \le 7{,}5$ s al 2 %
2. $M_p \le 15\%$
3. $e_{pert} < 5\%$ en régimen permanente ante una perturbación de carga constante del 20 % de la referencia ($N_0 = 0{,}2$)

**B.b**: "Diseñe un controlador utilizando el método del root locus que permita cumplir **simultáneamente** con los requerimientos del punto anterior." — *los tres*.

**B.d / B.e**: evaluar rampa; si no hay seguimiento con error nulo, modificar el controlador.

**Ej.2 Parte B — Tres compensadores lead** para MF = 35°, 50°, 65°, normalizados a $C(0)=1$. Esto hoy ya cumple.

---

## 1. Diagnóstico del problema principal (Ej.1 B.b)

El diseño actual es un **lead puro**. Matemáticamente, para este plante $G(s) = K_t / [s\,(Ls+R)(Js+B) + K_e K_t\,s]$ el error permanente ante perturbación constante $N_0$ vale:

$$
e_{pert} = \frac{0{,}64}{C(0)}\quad\text{(en \%, con }N_0=0{,}2\text{)}
$$

Para cumplir $e_{pert} < 5\%$ se necesita $C(0) > 12{,}8$. El lead actual da $C(0) \approx 7{,}6$, → $e_{pert} \approx 8{,}4\%$. **No cumple**.

Un lead puro no puede subir $C(0)$ arbitrariamente manteniendo $t_s$ y $M_p$: cuando el polo se aleja lo suficiente para aportar fase grande, la ganancia DC se degrada y $K$ no compensa sin perder $M_p$.

**El informe actual intenta justificar por qué "no se puede"**: esa justificación es incorrecta respecto a la consigna. La consigna pide *resolverlo*, no rendirse.

### Solución propuesta: **Lead-Lag**

Estructura:

$$
C(s) = K_c \cdot \underbrace{\frac{s+z_L}{s+p_L}}_{\text{lead: }p_L>z_L}\cdot\underbrace{\frac{s+z_{lag}}{s+p_{lag}}}_{\text{lag: }z_{lag}>p_{lag}}
$$

- **Lead**: aporta fase positiva en la banda media → mantiene $M_p$ y $t_s$.
- **Lag**: sube la ganancia en DC sin alterar la fase cerca de $\omega_{cp}$ (poniendo $z_{lag},\,p_{lag} \ll \omega_{cp}$). La ganancia DC del lag es $z_{lag}/p_{lag} > 1$.
- $C(0) = K_c \cdot (z_L/p_L) \cdot (z_{lag}/p_{lag})$. Elegir $z_{lag}/p_{lag}$ lo suficientemente grande para que $C(0) > 12{,}8$ (típicamente $C(0) \approx 20$ da buen margen, $e_{pert}\approx 3{,}2\%$).

Alternativa aceptable: **PI-lead** (integrador + lead). Hace $e_{pert}=0$ trivialmente, pero vuelve el sistema Tipo 2 y entonces B.e queda trivial. Se recomienda **preferir lead-lag** porque preserva el flujo pedagógico del TP (Tipo 1 en B.b, se vuelve Tipo 2 en B.e).

---

## 2. Plan de ejecución — Fase 1: Código (iterativo)

**Meta**: producir una función `disenar_lead_lag_gridsearch(G)` en `src/controlador.py` que devuelva un compensador que cumpla los tres specs, y que sea llamada desde `src/generar_figuras.py` y desde el notebook `src/notebooks/ejercicio1.ipynb`.

### 2.1 Implementar `disenar_lead_lag_gridsearch` en `src/controlador.py`

Reglas:

- **NO borrar** `disenar_lead_gridsearch` (lo usa el notebook y conviene mantenerlo como comparación).
- Agregar una función nueva con esta firma:
  ```python
  def disenar_lead_lag_gridsearch(G, C0_min=13.5, verbose=True):
      """
      Diseña un compensador lead-lag para cumplir ts<=7.5, Mp<=15, e_pert<5 (=>C(0)>12.8).
      Estrategia:
        1) Grid search sobre el lead (como disenar_lead_gridsearch).
        2) Para cada lead factible, grid search sobre el lag:
           - p_lag en [0.001, 0.05]
           - ratio_lag = z_lag/p_lag en [C0_min/C0_lead * 1.05, ...] (lo necesario para llegar a C0_min)
           - z_lag = ratio_lag * p_lag, con z_lag << omega_cp (e.g. z_lag <= omega_cp/10)
        3) Para cada candidato C = Kc * (s+z_L)/(s+p_L) * (s+z_lag)/(s+p_lag):
           - Recalcular Kc por condición de magnitud en s_d (el mismo de la etapa lead).
           - Simular lazo cerrado: ts, Mp, e_pert (numéricamente con step_response sobre G_n).
           - Filtrar: ts<=7.5, Mp<=15, e_pert<5, estable.
        4) Devolver el factible con mínimo Mp (o mínimo ts*Mp).
      """
  ```

- **IMPORTANTE**: calcular `e_pert` **numéricamente** (no usando `64/C(0)`). La fórmula analítica es una aproximación del valor final; con lag el transitorio de la perturbación puede ser lento y `step_response` debe correrse con `t_final` suficientemente largo (p.ej. 60-100 s) para atrapar el régimen permanente. Usar exactamente el mismo bloque que `evaluar_specs` para consistencia.

- La etapa lead del paso (1) puede reutilizar el código de `disenar_lead_gridsearch` pero **relajando la selección**: no descartar candidatos con $e_{pert}>5\%$ (porque ahora el lag lo va a resolver). En lugar de eso, guardar todos los leads factibles en $t_s, M_p$ y estables.

- Valores sugeridos de grilla:
  - Lead: `z_L ∈ arange(0.2, 8.0, 0.4)`, `sigma_d ∈ arange(0.7, 9.0, 0.4)`, `zeta_d ∈ arange(0.55, 0.92, 0.06)` (misma que hoy).
  - Lag: `p_lag ∈ [0.002, 0.005, 0.01, 0.02, 0.05]`, `ratio_lag ∈ [15, 20, 25, 30]`.
  - Son ~5000 combinaciones, corre en < 2 minutos.

### 2.2 Actualizar `src/generar_figuras.py`

Reemplazar el bloque de B.b que usa `disenar_lead_gridsearch` por el nuevo:

```python
from controlador import disenar_lead_lag_gridsearch
best, _ = disenar_lead_lag_gridsearch(G, C0_min=13.5, verbose=True)
z_L  = best['z_L']; p_L  = best['p_L']
z_lag = best['z_lag']; p_lag = best['p_lag']
Kc    = best['Kc']
C_final = best['C']
```

Y ajustar los `w('LeadZcero', ...)` en el bloque de valores LaTeX para que escriban los 5 parámetros del lead-lag: `\LeadZcero`, `\LeadPolo`, `\LagZcero`, `\LagPolo`, `\LeadGanancia`.

Mantener el bloque B.e (lead + integrador) tal cual: ya funciona.

### 2.3 Actualizar el notebook `src/notebooks/ejercicio1.ipynb`

Reemplazar la celda B.b (`947f7038` según el transcript) para que llame a `disenar_lead_lag_gridsearch` y muestre los 5 parámetros. Usar `NotebookEdit` tool.

### 2.4 Ciclo de iteración

Esta es **la parte crítica**. Loop hasta que todo cumpla:

```
REPEAT:
  1. Correr:  cd src && ../.venv/bin/python3 generar_figuras.py
  2. Leer el output impreso (NO mirar sólo "✓"):
     - ¿ts <= 7.5? ¿Mp <= 15? ¿e_pert < 5?
     - Si el script crashea: leer el traceback, arreglar el bug, volver a correr.
  3. SI CUMPLE: break.
  4. SI NO CUMPLE ts o Mp: aflojar grilla del lag (reducir ratio_lag máximo) o del lead.
  5. SI NO CUMPLE e_pert: aumentar ratio_lag mínimo (ej. empezar en 20 en vez de 15), o aumentar C0_min.
  6. SI no encuentra ningún factible: leer el mensaje de debug, ampliar grilla (más valores de sigma_d, de zeta_d, de p_lag).
UNTIL (ts <= 7.5 AND Mp <= 15 AND e_pert < 5).
```

No intentes "calcular" los valores a ojo — confiá en la grilla. Si la grilla no encuentra nada, amplíala. No toques el informe hasta que el código pase.

**Sanity check final** (antes de pasar a la fase 2): abrir `informe/valores_ej1.tex` y verificar que:
- `\LeadTs` sea `<= 7,50`
- `\LeadMp` sea `<= 15,0`
- `\LeadEpert` sea `< 5,0`
- `\LeadCumpleTs` y `\LeadCumpleMp` sean ambos `\checkmark`

Si alguno no, volver a iterar.

---

## 3. Plan de ejecución — Fase 2: Ejercicio 2

Ej.2 ya cumple numéricamente (verificar mentalmente con la última corrida):
- PM=35° obtenido: 37,3° ≥ 35 ✓
- PM=50° obtenido: 51,5° ≥ 50 ✓
- PM=65° obtenido: 65,6° ≥ 65 ✓
- $C(0)=1{,}00$ en los tres (por construcción, $K_c = p/z$)

**No hay nada que arreglar en el código de Ej.2.** Verificar que `valores_ej2.tex` existe y tiene valores razonables.

Si por algún motivo la última celda no coincide (ej. si se regenera y cambian las grillas), re-correr `../.venv/bin/python3 generar_figuras_ej2.py`.

---

## 4. Plan de ejecución — Fase 3: Actualizar el informe

Sólo arrancar esta fase cuando la Fase 1 cumpla todos los specs.

### 4.1 Archivos a editar

Sólo `informe/informe.tex`. Los `valores_ej*.tex` se autogeneran.

### 4.2 Secciones que necesitan re-escritura conceptual (no sólo números)

**Son las partes donde el texto actual explica/justifica algo que era consecuencia del diseño viejo (lead puro que no cumplía). Ahora que el diseño SÍ cumple, esa narrativa está obsoleta.**

Ubicar cada una leyendo el tex y reemplazar:

#### B.b — Descripción del método (líneas ~300-330)

Reemplazar el párrafo "Método: búsqueda en grilla con contribución angular" para que diga:

- Estructura del compensador: lead-lag ($C = K_c \frac{s+z_L}{s+p_L}\frac{s+z_{lag}}{s+p_{lag}}$).
- El lead fija transitorio ($t_s$, $M_p$) por contribución angular en el root locus.
- El lag sube la ganancia DC (ratio $z_{lag}/p_{lag}$ > 1) sin afectar la fase cerca de $\omega_{cp}$ (pues $z_{lag}, p_{lag} \ll \omega_{cp}$). Esto permite llegar a $C(0) > 12{,}8$ y cumplir $e_{pert} < 5\%$.
- Cómo se eligen: grid search anidado, primero lead factible en $t_s$ y $M_p$, luego lag con ratio suficiente.

Agregar la ecuación del compensador obtenido con los valores numéricos reales (vía macros):

```latex
C(s) = \LeadGanancia \cdot \frac{s + \LeadZcero}{s + \LeadPolo} \cdot \frac{s + \LagZcero}{s + \LagPolo}
```

#### B.b — Tabla de cumplimiento (ya parametrizada, sólo verificar que la columna "Cumple" del $e_{pert}$ use macro)

Cambiar la tercera fila para que use `\LeadCumpleEpert` (crear este macro también en `generar_figuras.py`):

```latex
$e_{pert}$ & \LeadEpert\,\% & $< 5\,\%$ & \LeadCumpleEpert \\
```

Agregar al script:
```python
w('LeadCumpleEpert', ck(res['cumple_pert']))
```

#### B.b — Párrafo "¿Por qué el lead no puede cumplir la especificación de perturbación?" (líneas ~371-386)

**Borrar completo este párrafo.** Ya no aplica: el diseño actual SÍ cumple. Reemplazarlo por un párrafo corto que explique *por qué el lead-lag sí cumple*:

> El lag agrega ganancia en DC ($z_{lag}/p_{lag} \approx X$) sin afectar la fase cerca de $\omega_{cp}$, lo que eleva $C(0)$ de $\approx 7{,}6$ (lead puro) a $\approx \LeadCeroDC$, cumpliendo $C(0) > 12{,}8$ y por ende $e_{pert} < 5\%$. La condición $z_{lag}, p_{lag} \ll \omega_{cp}$ asegura que el lag contribuye fase negativa despreciable en la banda de cruce, preservando $t_s$ y $M_p$.

(Crear macro `\LeadCeroDC` = $C(0)$ del diseño final.)

#### B.c — Bode comparativo (líneas ~388-412)

Ya está parametrizado. Verificar que el párrafo final que explica "el GM también mejora porque el polo del lead..." sea coherente con el lead-lag. Probablemente reemplazar por una descripción cualitativa del Bode lead-lag: pico de ganancia en baja frecuencia (del lag), fase positiva cerca de $\omega_{cp}$ (del lead).

#### B.d — Seguimiento de rampa (líneas ~414-429)

Sigue siendo válido: el sistema lead-lag mantiene Tipo 1, la rampa tiene error finito. Actualizar la fórmula de $K_v$ para incluir el lag:

$$K_v = K_c \cdot \frac{z_L}{p_L} \cdot \frac{z_{lag}}{p_{lag}} \cdot \frac{K_t}{a_0}$$

Usar `\LeadKv` y `\LeadEss` como ya están.

#### B.e — Integrador (líneas ~431-472)

Sin cambios estructurales. El flujo "añadir integrador para rampa" sigue siendo válido.

#### Conclusiones (líneas ~653-660)

Reemplazar el párrafo de Ej.1:

> **Ejercicio 1**: se partió de un modelo completo del motor DC (3er orden, Tipo 1) y se diseñó un compensador **lead-lag** por root locus que satisface **simultáneamente** las especificaciones de transitorio ($t_s = \LeadTs$ s, $M_p = \LeadMp$ \%) y de rechazo de perturbación ($e_{pert} = \LeadEpert$ \% $< 5\%$). El lag fue clave para elevar la ganancia DC sin degradar los márgenes transitorios. Para seguimiento de rampa con error nulo se añadió un integrador (Tipo 2); el sobrepico resultante ($M_p = \IntMp$ \%) supera el límite del 15 %, evidenciando el compromiso entre seguimiento de alta tipo y amortiguamiento.

### 4.3 Reglas al editar

- **NO** hardcodear números en el tex. Siempre usar macros `\Lead*`, `\PM*`, etc.
- Si hacés falta un macro nuevo, agregalo en `generar_figuras*.py` primero, corré el script, verificá que aparezca en `valores_ej*.tex`, y después referencialo en el tex.
- Usá `Edit` del tex con contextos grandes para evitar ambigüedad.
- No reescribas secciones que ya están bien (A.a, A.b, A.c, A.d, A.e del Ej.1; casi todo el Ej.2). Focalizate en las secciones listadas arriba.

### 4.4 Compilar y verificar

```
cd informe
pdflatex -interaction=nonstopmode informe.tex
pdflatex -interaction=nonstopmode informe.tex   # segunda pasada para refs
```

Buscar en el stdout "Undefined control sequence" (macro usado sin definir): si aparece, agregar el macro al script Python y re-correr.

Abrir el PDF (o usar `pdftotext` + grep) y verificar que:
- La tabla de B.b muestra los tres ✓ en la columna Cumple.
- Las tablas del Ej.2 tienen valores reales (no quedan `(código)` sueltos).
- Las frases del informe no contradicen los valores (p.ej. que no diga "$e_{pert} > 5\%$" cuando el valor real es 3 %).

---

## 5. Plan de ejecución — Fase 4: verificación cruzada final

Al terminar las fases 1-3, hacer un barrido:

1. `grep -n "(ver salida código)\|(ver código)\|(código)\|ver salida"` en `informe/informe.tex` → debe devolver 0 resultados.
2. `grep -n "TODO\|FIXME\|XXX"` en el tex → debe devolver 0 o sólo los TODO de portada (nombres).
3. Abrir `informe/valores_ej1.tex` y `informe/valores_ej2.tex`: todos los `\def` deben tener valores finitos y razonables.
4. Re-compilar el PDF dos veces. No debe haber errores ni "Undefined control sequence".
5. Mirar el PDF (al menos spot-checks): tablas B.b, B.c, B.e, tabla PM del Ej.2, conclusiones.

---

## 6. Reglas de seguridad (para no romper cosas)

- **No** hacer `git commit`, `git push`, ni borrar archivos.
- **No** tocar `src/planta.py`, `src/planta2.py`, `src/analisis.py` (modelos físicos, no cambian).
- **No** regenerar figuras que ya están bien (bode_planta, escalon_perturbacion, rampa_perturbacion, rlocus_sin_compensar, todas las de Ej.2). Sólo las que dependen de `C_final` deben regenerarse al cambiar el compensador — `generar_figuras.py` ya las sobreescribe todas, está OK.
- Si surge un error en `ctrl.margin` o `ctrl.feedback`, es probable que sea porque la simplificación del control dio un polo en el origen u otro caso degenerado. Envolver en `try/except` en el grid search y seguir.
- No toques los `\def` auto-generados manualmente — siempre regenerar vía el script.
- **No** aumentar la granularidad de la grilla sin haber visto primero la salida: si la grilla actual ya produce >100 factibles, buscar está bien afinado, basta con reordenar el criterio de selección.

---

## 7. Orden de ejecución recomendado

1. Fase 1.1: implementar `disenar_lead_lag_gridsearch` en `controlador.py`.
2. Fase 1.2: modificar `generar_figuras.py`.
3. Fase 1.4: correr + iterar hasta que los tres specs cumplan.
4. Fase 1.3: actualizar el notebook (una sola vez, al final de la fase 1).
5. Fase 2: verificar que Ej.2 sigue OK (una corrida).
6. Fase 3: reescribir las secciones del informe listadas.
7. Fase 4: verificación cruzada, compilar PDF.

Tiempo estimado: 45-90 minutos si todo va bien.

---

## 8. Criterio de terminación

Listo cuando, ejecutando desde cero:

```
cd src
../.venv/bin/python3 generar_figuras.py
../.venv/bin/python3 generar_figuras_ej2.py
cd ../informe
pdflatex -interaction=nonstopmode informe.tex
pdflatex -interaction=nonstopmode informe.tex
```

...produce `informe.pdf` sin errores, con:
- Tabla B.b del Ej.1: ts ≤ 7,5 ✓, Mp ≤ 15 ✓, e_pert < 5 ✓ (los tres ticks).
- Tabla B.c: márgenes del compensado poblados con números reales.
- Tablas Ej.2 B.a y B.b: seis filas con los parámetros reales de los tres compensadores.
- Conclusiones: narrativa coherente con los números nuevos (el Ej.1 "cumple los tres simultáneamente").
- Cero ocurrencias de `(ver código)`, `(ver salida código)`, `(código)` en el tex.
