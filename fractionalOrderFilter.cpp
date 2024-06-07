#include "fractionalOrderFilter.h"
#include <cmath>
#include <algorithm>
#include <cstring> // Para memmove

FractionalOrderFilter::FractionalOrderFilter(double orden, double periodo, double tipo) {
     if (orden < MIN_ORDER || orden > MAX_ORDER || tipo < MIN_TYPE || tipo > MAX_TYPE) {
        return;
    }
    gC.orden = orden;
    gC.periodo = periodo;
    gC.tipo = tipo;
    Escala.valor = pow(gC.periodo, -gC.orden);
    Escala.aplicar = true;
}

void FractionalOrderFilter::gaussJordan(double* A, double* b, int N) {
    for (int i = 0; i < N; i++) {
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

void FractionalOrderFilter::calcularMatriz(const MatrixParams& params) {
    for (int i = 0; i < params.filas; i++) {
        for (int j = 0; j <= i; j++) {
            double sumMatrix = 0.0;
            double sumResult = 0.0;
            for (int k = params.inicio; k < params.fin; k++) {
                double val_i = (k - i >= 0) ? params.g_o_h[k - i] : 0.0;
                double val_j = (k - j >= 0) ? params.g_o_h[k - j] : 0.0;
                sumMatrix += val_i * val_j;
                if (j == 0) {
                    sumResult += val_i * (params.positivo ? params.h[k] : params.h[k + 1]);
                }
            }
            params.matriz[params.filas * i + j] = sumMatrix;
            if (i != j) {
                params.matriz[params.filas * j + i] = sumMatrix;
            }
            if (j == 0) {
                params.result[i] = (params.positivo ? sumResult : -sumResult);
            }
        }
    }
}

void FractionalOrderFilter::ponerCeros(double* arreglo,int size){
    std::fill(arreglo, arreglo + size, 0.0);
}

void FractionalOrderFilter::respuestaImpulso(int cantidad, double* resultado) {
    if (cantidad <= 0) {
        return;
    }
    double k1 = pow(gC.tipo, -gC.orden);
    double k2 = (1.0 - gC.tipo) / gC.tipo;
    double inv_k2 = 1/k2;
    double k2_power = 1.0;
    double c1 = 1 + gC.orden;
    double c2 = 1 - gC.orden;
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
        if (gC.tipo == 1){
            sum = EE[n];
        }else{
            double k2_local = k2_power;
            for (int k = 0; k <= n; k++) {
                sum += k2_local * EE[k] * DD[n - k];
                k2_local *=inv_k2;
            }
        }
        resultado[n] =k1*sum*((Escala.aplicar == true) ? Escala.valor : 1);
    }
}

void FractionalOrderFilter::Pade(double* num, int ordenNum, double* den, int ordenDen) {
    Escala.aplicar = false;
    double h[ordenNum + ordenDen + 1];
    respuestaImpulso(ordenNum + ordenDen + 1, h);
    Pade(h, num, ordenNum, den, ordenDen);
}

void FractionalOrderFilter::Pade(double* h, double* num, int ordenNum, double* den, int ordenDen) {
    double dH[ordenDen * ordenDen];
    for (int i = 0; i < ordenDen; i++) {
        for (int j = 0; j < ordenDen; j++) {
            dH[ordenDen * i + j] = h[ordenNum + i - j];
        }
    }
    for (int i = 0; i < ordenDen; i++) {
        den[i] = -h[ordenNum + 1 + i];
    }
    gaussJordan(dH,den,ordenDen);
    //memmove(&den[1], den, ordenDen * sizeof(double));
    std::copy(den, den + ordenDen, &den[1]);
    den[0] = 1;
    numPade(h,num,ordenNum,den,ordenDen);
}

void FractionalOrderFilter::numPade(double* h, double* num,double ordenNum,double* den,double ordenDen){
    for (int n = 0; n <= ordenNum; n++) {
        int limite = (n < ordenDen) ? n : ordenDen;
        for (int k = 0; k <= limite; k++) {
            num[n] += den[k] * h[n - k];
        }
       num[n] = num[n] * Escala.valor;
    }
}

void FractionalOrderFilter::Shank(double* num, int ordenNum, double* den, int ordenDen, int cantidad) {
    ponerCeros(num,ordenNum);
    ponerCeros(den,ordenDen);
    if (cantidad<=ordenNum+ordenDen){
        return;
    }else if(cantidad == ordenNum+ordenDen +1){
        Pade(num,ordenNum,den,ordenDen);
    }else{
        Escala.aplicar = false;
        double h[cantidad]={0.0};
        respuestaImpulso(cantidad,h);
        Shank(h, num, ordenNum, den, ordenDen,cantidad);
        numPade(h, num, ordenNum, den, ordenDen);
    }
}

void FractionalOrderFilter::Shank(double* h, double* num, int ordenNum, double* den, int ordenDen, int cantidad) {
    int sHTH = ordenDen * ordenDen;
    double HTH[sHTH];
    ponerCeros(HTH,sHTH);
    MatrixParams params = {HTH, h, h, ordenDen, ordenNum, cantidad - 1, den, false};
    calcularMatriz(params);
    gaussJordan(HTH, den, ordenDen);
    //memmove(&den[1], den, ordenDen * sizeof(double));
    std::copy(den, den + ordenDen, &den[1]);
    den[0] = 1.0;
}


void FractionalOrderFilter::Prony(double* num, int ordenNum, double* den, int ordenDen, int cantidad) {
    ponerCeros(num,ordenNum);
    ponerCeros(den,ordenDen);
    if (cantidad<=ordenNum+ordenDen){
        return;
    }else if(cantidad == ordenNum+ordenDen +1){
        Pade(num,ordenNum,den,ordenDen);
    }else{
        Escala.aplicar = false;
        double h[cantidad];
        respuestaImpulso(cantidad, h);
        Shank(h, num, ordenNum, den, ordenDen,cantidad);
        Prony(h, num, ordenNum, den, ordenDen,cantidad);
        numPade(h, num, ordenNum, den, -1);
    }
}


void FractionalOrderFilter::Prony(double* h, double* num, int ordenNum, double* den, int ordenDen, int cantidad) {
    double g[cantidad];
    g[0] = 0.0;
    for (int n = 0; n < cantidad; n++) {
        int limite = (n < ordenDen) ? n : ordenDen;
        for (int k = 0; k <= limite; k++) {
            g[n] -= den[k] * g[n - k];
        }
        if (n==0)
        g[n] += 1;
    }
    int sGTG = (ordenNum + 1) * (ordenNum + 1);
    double GTG[sGTG];
    ponerCeros(GTG,sGTG);
    MatrixParams params = {GTG, h, g, ordenNum + 1, 0, cantidad, num, true};
    calcularMatriz(params);
    gaussJordan(GTG, num, ordenNum + 1);
}