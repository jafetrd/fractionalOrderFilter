# fractionalOrderFilter

Esta librería realiza el cálculo de la respuesta impulso de operadores fraccionarios utilizando aproximaciones de primer orden como Euler, Al-Alaoui y Tustin. A partir de estos coeficientes, que pueden ser usados para un filtro FIR, se generan los coeficientes para un filtro tipo IIR mediante tres métodos: Padé, Shank y Prony.

## Características

- **Aproximaciones de Primer Orden**:
  - Euler
  - Al-Alaoui
  - Tustin

- **Métodos de Conversión FIR a IIR**:
  - Padé: La peor aproximación en términos de precisión, pero la menos costosa computacionalmente y sin cota de error.
  - Shank: Calcula los coeficientes del numerador de forma óptima en el sentido cuadrático.
  - Prony: Calcula los coeficientes del numerador y denominador en el sentido cuadratico, siendo el más costoso computacionalmente.

## Requisitos

- C++11 o superior

## Instalación

Para utilizar esta librería, simplemente clona el repositorio y compila el archivo `main.cpp`.

```bash
git clone <https://github.com/jafetrd/fractionalOrderFilter>
cd <fractionalOrderFilter>
g++ -std=c++11 main.cpp -o main

