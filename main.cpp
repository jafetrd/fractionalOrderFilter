#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include "math.h"
#include "fractionalOrderFilter.h"

// Función para guardar un vector de datos en un archivo
void guardarVector(const std::string& filename, const double* data, int size) {
    std::ofstream file(filename);
    file << "{";
    for (int i = 0; i < size - 1; ++i) {
        file << std::setprecision(15) << data[i] << ",";
    }
    file << std::setprecision(15) << data[size - 1] << "};";
}

// Función para guardar los coeficientes de numerador y denominador en un archivo
void guardarCoeficientes(const std::string& filename, const double* num, int ordenNum, const double* den, int ordenDen) {
    std::ofstream file(filename);
    file << "{{";
    for (int i = 0; i < ordenNum; ++i) {
        file << std::setprecision(15) << num[i] << ",";
    }
    file << std::setprecision(15) << num[ordenNum] << "},{";
    for (int i = 0; i < ordenDen; ++i) {
        file << std::setprecision(15) << den[i] << ",";
    }
    file << std::setprecision(15) << den[ordenDen] << "}}";
}

int main() {
    // Inicialización del filtro con parámetros comunes
    FractionalOrderFilter filter(-0.5, 0.01, 0.5);

   // Ejemplo de coeficientes de la respuesta impulso
    int cantidad1 = 100;
    double h[cantidad1];
    filter.respuestaImpulso(cantidad1, h);
    guardarVector("datos1.txt", h, cantidad1);

    // Ejemplo de coeficientes de Padé
    int ordenNum1 = 9, ordenDen1 = 9;
    double num1[ordenNum1 + 1], den1[ordenDen1 + 1];
    filter.Pade(num1, ordenNum1, den1, ordenDen1);
    guardarCoeficientes("datos2.txt", num1, ordenNum1, den1, ordenDen1);

    // Ejemplo de coeficientes de Shank
    int ordenNum3 = 9, ordenDen3 = 9, cantidad3 = 1000;
    double num3[ordenNum3 + 1], den3[ordenDen3 + 1];
    filter.Shank(num3, ordenNum3, den3, ordenDen3, cantidad3);
    guardarCoeficientes("datos3.txt", num3, ordenNum3, den3, ordenDen3);
    
    // Ejemplo de coeficientes de Prony
    int ordenNum2 = 9, ordenDen2 = 9, cantidad2 = 1000;
    double num2[ordenNum2 + 1], den2[ordenDen2 + 1];
    filter.Prony(num2, ordenNum2, den2, ordenDen2, cantidad2);
    guardarCoeficientes("datos4.txt", num2, ordenNum2, den2, ordenDen2); 

    return 0;
}
