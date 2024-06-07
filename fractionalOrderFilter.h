#ifndef FRACTIONALORDERFILTER_H
#define FRACTIONALORDERFILTER_H

#define MIN_ORDER -1.0
#define MAX_ORDER 1.0
#define MIN_TYPE 0
#define MAX_TYPE 1.5

class FractionalOrderFilter {
public:

    FractionalOrderFilter(double orden, double periodo, double tipo);
    void respuestaImpulso(int cantidad, double* resultado);
    void Pade(double* num, int ordenNum, double* den, int ordenDen);
    void Prony(double* num, int ordenNum, double* den, int ordenDen,int L);
    void Shank(double* num, int ordenNum, double* den, int ordenDen,int L);
    
private:

    bool desdeShank;
    
    struct Config {
        double orden;
        double periodo;
        double tipo;
    };

    Config gC;
    
    struct ganancia {
        double valor;
        bool aplicar;
    };
    
    ganancia Escala;
    
    struct MatrixParams {
        double* matriz;
        double* h;
        double* g_o_h;
        int filas;
        int inicio;
        int fin;
        double* result;
        bool positivo;
    };
    
    void calcularMatriz(const MatrixParams& params);
    void gaussJordan(double* A, double* b, int N);
    void Pade(double* h, double* num, int ordenNum, double* den, int ordenDen);
    void numPade(double*h, double* num,double ordenNum,double* den,double ordenDen);
    void Prony(double* h, double* num, int ordenNum, double* den, int ordenDen,int L);
    void Shank(double* h, double* num, int ordenNum, double* den, int ordenDen,int L);
    void ponerCeros(double* arrelgo,int size);
};

#endif // FRACTIONALORDERFILTER_H
