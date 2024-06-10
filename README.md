# fractionalOrderFilter

Esta librería realiza el cálculo de la respuesta impulso es decir los coeficientes de la serie de MacLaurin de un integrador sintonizable o de Smith elevado a una potencia real entre -1 y 1, haciendo variar un parametro de este, se obtienen las aproximaciones tipo Euler atrasado, Al-Alaoui y Tustin/bilineal que corresponde a la familia de aproximaciones discretas de primer orden para la integral o la derivada. A partir de estos coeficientes, que pueden ser usados para un filtro FIR, se generan los coeficientes para un filtro tipo IIR mediante tres métodos: Padé, Shank y Prony.

## Características

- **Aproximaciones de Primer Orden**:
  - Euler
  - Al-Alaoui
  - Tustin

- **Métodos de Conversión FIR a IIR**:
  - Padé: La peor aproximación en términos de precisión, pero la menos costosa computacionalmente y sin cota de error.
  - Shank: Calcula los coeficientes del numerador de forma óptima en el sentido cuadrático.
  - Prony: Calcula los coeficientes del numerador y denominador en el sentido cuadratico, siendo el más costoso computacionalmente.

# Documentación
Para hacer uso hay que agregarla al proyecto con:
 ``` #include "fractionalOrderFilter.h" ```

 Se declara el operador con 3 parametros
 ```FractionalOrderFilter(double orden, double periodo, double tipo)```
 
- orden: es un valor entre -1 y 1, donde valores negativos devuelven coeficientes para aproximar integrales y valores positivos son para derivadas
- periodo: es la base de tiempo en segundos de la señal a la que se le aplicará el filtrado, debe ser un valor mayor a cero
- tipo: es un valor entre 0.5 y 1 para elegir el tipo de operador

  * tipo = 1 -> Euler atrasado, usar solo para derivadas 
  * tipo = 0.5 -> Tustin/Bilineal, usar solo para integrales 
  * tipo = 0.875 -> Al-Alaoui, se puede usar en ambos al ser una interpolación entre Euler y Tustin

## Métodos
```C 
void _respuestaImpulso(int cantidad, double* respuesta);
```
Asigna al array **respuesta** el numero de coeficientes de MacLaurin indicado por **cantidad**, cuando tipo=1 los coeficientes generados son los de la conocida definición de Grunwald - Letnikov.

```C
void _Pade(double* num, int ordenNum, double* den, int ordenDen);
```
Asigna a los arrays **num** y **den** la cantidad de coeficientes mas uno, indicado por sus ordenes respectivos **ordenNum** y **ordenDen**, estos ordenes corresponden a los polinomios del numerador y denominador del filtro IIR $P_N^M(z^{-1})$ en el que serán usados. Internamente ejecuta el metodo **_respuestaImpulso** y mediante el método de Padé que extrae la mayor información posible de una serie de Taylor, genera los coeficientes. Es el menos costoso computacionalmente pero sin cota de error despues de ordenNum+ordenDen+1 valores generados por su respuesta impulso. 

```C
void _Shank(double* num, int ordenNum, double* den, int ordenDen, int cantidad);
```
El parametro **cantidad** recibe la el numero de valores de la repuestaImpulso que serán tomado en cuenta para generar una minimización cuadratica que optimizará los coeficientes del denominador del filtro IIR, los coeficientes del numerador se generan son el método de Padé.

```C
void _Prony(double* num, int ordenNum, double* den, int ordenDen, int cantidad);
```
Tanto los coeficientes del numerador como los del denominador son generados mediante una minimización cuadratica, representada como una pseudo inversa, es el mas costoso computacionalmente pero con mejores resultados. 

# Ejemplo de uso

```C
#include "fractionalOrderFilter.h"
int main() {
    // Definir un integrador de orden medio tipo Tustin 
    //con un tiempo entre muestras de 0.01 segundos 
    FractionalOrderFilter filter(-0.5, 0.01, 0.5);
    // Aproximarlo para un filtro IIR de orden 9 tanto en el 
    //numerador como el denominador tomando en cuenta 
    //1000 coeficientes de la respuesta impulso
    int ordenNum = 9, ordenDen = 9, cantidad = 1000;
    double num[ordenNum + 1], den[ordenDen + 1];
    filter._Prony(num, ordenNum, den, ordenDen, cantidad);
    return 0;
}
```

# Posibles problemas

Los tres metodos Padé, Shank y Prony son calculados usando el método de Gauss Jordan con pivoteo parcial para resolver ecuaciones simultáneas y pueden surgir sistemas que no tengan solución, los métodos de Shank y Prony hacen uso de pseudoinversas las cuales para resolverlas se representan como un sistema lineal $A^T A x = A^T b$ y se resuelven para $x$ y solo tiene solución si $A^T A$ se pude convertir en una matriz unitaria mediante operaciones por renglones. Otro problema que surge es cuando los ordenes del numerador y denominador se asignan a valores mayores a 15, lo que se traduce a un sistema de ecuaciones simultaneas de 15 x 15, la acumulación de errores numéricos provoca que no se llegue a ningun resultado. 

# Recomendaciones 
Para integrales hacer uso del operador Tustin y para las derivadas el operador de Euler atrasado, son los que generan mejores resultados 

# Contribuciones
Las contribuciones son bienvenidas. Por favor, abre un issue o envía un pull request para contribuir al desarrollo de esta librería.

# Bibliografia
El codigo es una implementación de la representación matricial mostrada por: 
Ramiro S. Barbosa, J.A. Tenreiro Machado, Manuel F. Silva,
Time domain design of fractional differintegrators using least-squares,
Signal Processing,
Volume 86, Issue 10,
2006,
Pages 2567-2581,
ISSN 0165-1684,
<https://doi.org/10.1016/j.sigpro.2006.02.005.>
<https://www.sciencedirect.com/science/article/pii/S0165168406000430>

## Requisitos

- C++11 o superior

## Instalación

Para utilizar esta librería, simplemente clona el repositorio y compila el archivo `main.cpp`.

```bash
git clone https://github.com/jafetrd/fractionalOrderFilter
cd fractionalOrderFilter
g++ -std=c++11 main.cpp -o main
