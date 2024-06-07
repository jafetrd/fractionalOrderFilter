#ifndef FRACTIONAL_ORDER_FILTER_H
#define FRACTIONAL_ORDER_FILTER_H

#include <vector>
#include <string>

const double _MIN_ORDER = -1.0;
const double _MAX_ORDER = 1.0;
const double _MIN_TYPE = 0.0;
const double _MAX_TYPE = 1.0;

struct _GlobalConfig {
    double _orden;
    double _periodo;
    double _tipo;
};

struct _Scale {
    double _valor;
    bool _aplicar;
};

struct _MatrixParams {
    double* _matriz;
    double* _h;
    double* _g_o_h;
    int _filas;
    int _inicio;
    int _fin;
    double* _result;
    bool _positivo;
};

class FractionalOrderFilter {
public:
    FractionalOrderFilter(double orden, double periodo, double tipo);
    void _respuestaImpulso(int cantidad, double* respuesta);
    void _Pade(double* num, int ordenNum, double* den, int ordenDen);
    void _Shank(double* num, int ordenNum, double* den, int ordenDen, int cantidad);
    void _Prony(double* num, int ordenNum, double* den, int ordenDen, int cantidad);
    
private:
    _GlobalConfig _gC;
    _Scale _Escala;
    
    void _numPade(double* h, double* num,double ordenNum,double* den,double ordenDen);
    void _Pade(double* h, double* num, int ordenNum, double* den, int ordenDen);
    void _Shank(double* h, double* num, int ordenNum, double* den, int ordenDen, int cantidad);
    void _Prony(double* h, double* num, int ordenNum, double* den, int ordenDen, int cantidad);
    void _gaussJordan(double* A, double* b, int N);
    void _calcularMatriz(const _MatrixParams& params);
    void _ponerCeros(double* array, int size);
};

#endif
