# Acoplamiento Magnético en 7-AGNR

Este repositorio contiene los scripts y cálculos para estudiar el acoplamiento magnético de intercambio ($J$) en nanodiscos y cintas de grafeno finitas de tipo 7-AGNR (Armchair Graphene Nanoribbons).

## Contenido

- **`Main/sisl_hubbard_7agnr.py`**: Es el script principal del proyecto. Se encarga de construir la geometría de las cintas de diferentes longitudes, definir el Hamiltoniano de Tight-Binding y calcular $J$ usando dos aproximaciones distintas para poder compararlas:
  1. **Dímero de Hubbard Efectivo**: Un modelo simplificado que extrae $t_{\mathrm{eff}}$ y $U_{\mathrm{eff}}$ para dar $J = 4t_{\mathrm{eff}}^2 / U_{\mathrm{eff}}$.
  2. **Hubbard en Campo Medio (MFH)**: Un cálculo self-consistent iterativo que evalúa la diferencia de energía real entre los estados Ferromagnético (FM) y Antiferromagnético (AFM).
- **`Main/23-06-2026.py`**: Una versión más de "andar por casa" que usa solo `numpy` y `scipy` para comprobar la regla del dímero y ver la probabilidad de los estados topológicos.
- **`Main/plot_spin_polarization_en.py`**: Script para visualizar cómo se distribuyen los espines en los bordes en zigzag de la cinta.
- **`Main/test_geom.py`**: Un pequeño test para asegurar que al cortar la cinta no nos hemos dejado enlaces colgantes indeseados.

## Resultados

Al ejecutar el script principal, genera dos gráficas directamente en tu carpeta:
1. **`J_comparison_vs_L.png`**: Compara cómo decae el acoplamiento magnético $J$ exponencialmente al hacer la cinta más larga ($L$). Le hemos añadido líneas de referencia súper útiles como el gap superconductor $2\Delta$ del Nb(110) a 3.0 meV y la energía térmica de un STM operando a 1.2 K. ¡Así se ve a simple vista si $J$ cae dentro del gap o no!
2. **`J_diff_vs_L.png`**: Una gráfica para ver exactamente cuánta diferencia de energía (el error o la desviación) hay entre el modelo sencillo del Dímero y la simulación completa de MFH.

## Requisitos

Para poder ejecutar todo sin problemas, tu entorno virtual de Python necesita tener instalados:
- `numpy` y `scipy`
- `matplotlib` (para las figuras)
- `sisl` (es la librería clave que usamos para definir las geometrías complejas y los hamiltonianos de forma sencilla)

¡Espero que el código quede claro! Tienes un montón de comentarios dentro explicando el paso a paso.
