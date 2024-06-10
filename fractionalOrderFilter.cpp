#include "fractionalOrderFilter.h"
#include <cmath>
#include <algorithm>

FractionalOrderFilter::FractionalOrderFilter(double orden, double periodo, double tipo) {
     if (orden < _MIN_ORDER || orden > _MAX_ORDER || tipo < _MIN_TYPE || tipo > _MAX_TYPE) {
        return;
    }
    _gC._orden = orden;
    _gC._periodo = periodo;
    _gC._tipo = tipo;
    _Escala._valor = pow(_gC._periodo, -_gC._orden);
    _Escala._aplicar = true;
}

void FractionalOrderFilter::_swapRows(double* A, double* b, int N, int row1, int row2) {
    if (row1 == row2) return;
    for (int i = 0; i < N; i++) {
        double temp = A[row1 * N + i];
        A[row1 * N + i] = A[row2 * N + i];
        A[row2 * N + i] = temp;
    }
    double temp = b[row1];
    b[row1] = b[row2];
    b[row2] = temp;
}

void FractionalOrderFilter::_gaussJordan(double* A, double* b, int N) {
    for (int i = 0; i < N; i++) {
         // Pivoteo parcial
        int maxRow = i;
        for (int k = i + 1; k < N; k++) {
            if (fabs(A[k * N + i]) > fabs(A[maxRow * N + i])) {
                maxRow = k;
            }
        }
        _swapRows(A, b, N, i, maxRow);
         
        double pivot = A[i * N + i];
        for (int j = i; j < N; j++) {
            A[i * N + j] /= pivot;
        }
        b[i] /= pivot;
        for (int j = 0; j < N; j++) {
            if (j != i) {
                double factor = A[j * N + i];
                for (int k = i; k < N; k++) {
                    A[j * N + k] -= factor * A[i * N + k];
                }
                b[j] -= factor * b[i];
            }
        }
    }
}

void FractionalOrderFilter::_calcularMatriz(const _MatrixParams& params) {
    for (int i = 0; i < params._filas; i++) {
        for (int j = 0; j <= i; j++) {
            double sumMatrix = 0.0;
            double sumResult = 0.0;
            for (int k = params._inicio; k < params._fin; k++) {
                double val_i = (k - i >= 0) ? params._g_o_h[k - i] : 0.0;
                double val_j = (k - j >= 0) ? params._g_o_h[k - j] : 0.0;
                sumMatrix += val_i * val_j;
                if (j == 0) {
                    sumResult += val_i * (params._positivo ? params._h[k] : params._h[k + 1]);
                }
            }
            params._matriz[params._filas * i + j] = sumMatrix;
            if (i != j) {
                params._matriz[params._filas * j + i] = sumMatrix;
            }
            if (j == 0) {
                params._result[i] = (params._positivo ? sumResult : -sumResult);
            }
        }
    }
}

void FractionalOrderFilter::_ponerCeros(double* array, int size){
    std::fill(array, array + size, 0.0);
}

void FractionalOrderFilter::_respuestaImpulso(int cantidad, double* respuesta) {
    if (cantidad <= 0) {
        return;
    }
    double k1 = pow(_gC._tipo, -_gC._orden);
    double k2 = (1.0 - _gC._tipo) / _gC._tipo;
    double inv_k2 = 1/k2;
    double k2_power = 1.0;
    double c1 = 1 + _gC._orden;
    double c2 = 1 - _gC._orden;
    double EE1 = 1.0;
    double DD[cantidad];
    double EE[cantidad];
    DD[0] = 1.0;
    EE[0] = 1.0;
    
    for (int n = 0; n < cantidad; n++) {
        if (n > 0){
            EE[n] = (1 - c1 / n) * EE[n - 1];
            EE1 *= (1 - c2 / n);
            DD[n] = EE1 *  (1 - 2 * (n & 1));
            k2_power *= k2;
        }
        double sum = 0.0;
        if (_gC._tipo == 1){
            sum = EE[n];
        }else{
            double k2_local = k2_power;
            for (int k = 0; k <= n; k++) {
                sum += k2_local * EE[k] * DD[n - k];
                k2_local *=inv_k2;
            }
        }
        respuesta[n] =k1*sum*((_Escala._aplicar == true) ? _Escala._valor : 1);
    }
}

void FractionalOrderFilter::_Pade(double* num, int ordenNum, double* den, int ordenDen) {
    _Escala._aplicar = false;
    double h[ordenNum + ordenDen + 1];
    _respuestaImpulso(ordenNum + ordenDen + 1, h);
    _Pade(h, num, ordenNum, den, ordenDen);
}

void FractionalOrderFilter::_Pade(double* h, double* num, int ordenNum, double* den, int ordenDen) {
    double dH[ordenDen * ordenDen];
    for (int i = 0; i < ordenDen; i++) {
        for (int j = 0; j < ordenDen; j++) {
            dH[ordenDen * i + j] = h[ordenNum + i - j];
        }
    }
    for (int i = 0; i < ordenDen; i++) {
        den[i] = -h[ordenNum + 1 + i];
    }
    _gaussJordan(dH,den,ordenDen);
    std::copy(den, den + ordenDen, &den[1]);
    den[0] = 1;
    _numPade(h,num,ordenNum,den,ordenDen);
}

void FractionalOrderFilter::_numPade(double* h, double* num,double ordenNum,double* den,double ordenDen){
    for (int n = 0; n <= ordenNum; n++) {
        int limite = (n < ordenDen) ? n : ordenDen;
        for (int k = 0; k <= limite; k++) {
            num[n] += den[k] * h[n - k];
        }
       num[n] = num[n] * _Escala._valor;
    }
}

void FractionalOrderFilter::_Shank(double* num, int ordenNum, double* den, int ordenDen, int cantidad) {
    _ponerCeros(num,ordenNum);
    _ponerCeros(den,ordenDen);
    if (cantidad<=ordenNum+ordenDen){
        return;
    }else if(cantidad == ordenNum+ordenDen +1){
        _Pade(num,ordenNum,den,ordenDen);
    }else{
        _Escala._aplicar = false;
        double h[cantidad];
        _respuestaImpulso(cantidad,h);
        _Shank(h, num, ordenNum, den, ordenDen,cantidad);
        _numPade(h, num, ordenNum, den, ordenDen);
    }
}

void FractionalOrderFilter::_Shank(double* h, double* num, int ordenNum, double* den, int ordenDen, int cantidad) {
    int sHTH = ordenDen * ordenDen;
    double HTH[sHTH];
    _ponerCeros(HTH,sHTH);
    _MatrixParams params = {HTH, h, h, ordenDen, ordenNum, cantidad - 1, den, false};
    _calcularMatriz(params);
    _gaussJordan(HTH, den, ordenDen);
    std::copy(den, den + ordenDen, &den[1]);
    den[0] = 1.0;
}


void FractionalOrderFilter::_Prony(double* num, int ordenNum, double* den, int ordenDen, int cantidad) {
    _ponerCeros(num,ordenNum);
    _ponerCeros(den,ordenDen);
    if (cantidad<=ordenNum+ordenDen){
        return;
    }else if(cantidad == ordenNum+ordenDen +1){
        _Pade(num,ordenNum,den,ordenDen);
    }else{
        _Escala._aplicar = false;
        double h[cantidad];
        _respuestaImpulso(cantidad, h);
        _Shank(h, num, ordenNum, den, ordenDen,cantidad);
        _Prony(h, num, ordenNum, den, ordenDen,cantidad);
        _numPade(h, num, ordenNum, den, -1);
    }
}


void FractionalOrderFilter::_Prony(double* h, double* num, int ordenNum, double* den, int ordenDen, int cantidad) {
    double g[cantidad];
    g[0] = 1.0;  
    for (int n = 1; n < cantidad; n++) {  
        g[n] = 0.0; 
        for (int k = 0; k <= ((n < ordenDen) ? n : ordenDen); k++) {
            g[n] -= den[k] * g[n - k];
        }
    }

    int sGTG = (ordenNum + 1) * (ordenNum + 1);
    double GTG[sGTG];
    _ponerCeros(GTG,sGTG);
    _MatrixParams params = {GTG, h, g, ordenNum + 1, 0, cantidad, num, true};
    _calcularMatriz(params);
    _gaussJordan(GTG, num, ordenNum + 1);
}
