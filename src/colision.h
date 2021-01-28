#ifndef COLISION_H
#define COLISION_H

#define VECTOR_SIZE 3
#define ARRAY_SIZE 13

struct Punto
{
    double x[3];
    double y[3];
    long D;
    
    double n;
    
    double a;
    double b;
    
    double c;
    double d;
    
    double e;
    double f;
};


struct Tren
{
    
    /* ESTOS ATRIBUTOS DEBEN INICIARSE DESDE PYTHON */
    double M;
    double Ixx;
    double Iyy;
    double Izz;
    double Ixz;
    
    double FC[VECTOR_SIZE];
    double MC[VECTOR_SIZE];
    double FB[VECTOR_SIZE];
    double MB[VECTOR_SIZE];

    double x[ARRAY_SIZE];
    double t;

    double T0;

    struct Punto ** P;
    long nP;
    /* FIN DE ATRIBUTOS */

    // Estado del integrador
    long colision;
    
    double * x0;
    double * x1;
    double * x2;
    
    double * xa;
    double * xi;

    double * F0;
    double * F1;
    double * F2;
    double * Fi;
    
    double  * Fa;

    double ti;
    double t0;
    double t1;
    double t2;

    double dt1;
    double dt2;

};

struct Tren * malloc_tren(long);
void free_tren(struct Tren *);
void step_tren(struct Tren *, double );
void F(double *, double, double *, struct Tren *);
#endif
