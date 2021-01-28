#include "algebra.h"

/*****************************************************************************/
/*  CREACION DE VECTORES Y MATRICES                                          */
/*****************************************************************************/
double *newvector() {
    return malloc(sizeof(double)*3);
}

double *newmatrix() {
    return malloc(sizeof(double)*9);
}

/*****************************************************************************/
/*	OPERACIONES DE MATRICES						                             */
/*****************************************************************************/
void setM(double *M, int i, int j, double v) {
    M[3*i + j] = v;
}

double getM(double *M, int i, int j) {
    return M[3*i + j];
}

void identity(double *R)
{
    R[0] = 1.0;
    R[1] = 0.0;
    R[2] = 0.0;
    
    R[3] = 0.0;
    R[4] = 1.0;
    R[5] = 0.0;
    
    R[6] = 0.0;
    R[7] = 0.0;
    R[8] = 1.0;
}

void dotSM(double s, double *M, double *R)
{
    int i;
    for (i=0;i<9;++i) {
        R[i] = s*M[i];
    }
}

void dotMM(double *M, double *N, double *R)
{
    R[0] = M[0]*N[0] + M[1]*N[3] + M[2]*N[6];
    R[1] = M[0]*N[1] + M[1]*N[4] + M[2]*N[7];
    R[2] = M[0]*N[2] + M[1]*N[5] + M[2]*N[8];

    R[3] = M[3]*N[0] + M[4]*N[3] + M[5]*N[6];
    R[4] = M[3]*N[1] + M[4]*N[4] + M[5]*N[7];
    R[5] = M[3]*N[2] + M[4]*N[5] + M[5]*N[8];

    R[6] = M[6]*N[0] + M[7]*N[3] + M[8]*N[6];
    R[7] = M[6]*N[1] + M[7]*N[4] + M[8]*N[7];
    R[8] = M[6]*N[2] + M[7]*N[5] + M[8]*N[8];
}

void addMM(double *M, double *N, double *R)
{
    int i;
    for (i=0; i<9; ++i) {
        R[i] = M[i] + N[i];
    }
}

void subMM(double *M, double *N, double *R)
{
    int i;
    for (i=0; i<9; ++i) {
        R[i] = M[i] - N[i];
    }
}

void invertM(double *M, double *R)
{
    double M0 = M[0];
    double M1 = M[1];
    double M2 = M[2];
    double M3 = M[3];
    double M4 = M[4];
    double M5 = M[5];
    double M6 = M[6];
    double M7 = M[7];
    double M8 = M[8];

    double det = M0*(M4*M8 - M5*M7) +
                 M1*(M5*M6 - M3*M8) + 
                 M2*(M3*M7 - M4*M6);
    
    R[0] = (M4*M8 - M5*M7)/det;
    R[1] = (M2*M7 - M1*M8)/det;
    R[2] = (M1*M5 - M2*M4)/det;
    
    R[3] = (M5*M6 - M3*M8)/det;
    R[4] = (M0*M8 - M2*M6)/det;
    R[5] = (M2*M3 - M0*M5)/det;
    
    R[6] = (M3*M7 - M4*M6)/det;
    R[7] = (M1*M6 - M0*M7)/det;
    R[8] = (M0*M4 - M1*M3)/det;
}

void transposeM(double *M, double *R)
{
    R[0] = M[0];
    R[1] = M[3];
    R[2] = M[6];
    
    R[3] = M[1];
    R[4] = M[4];
    R[5] = M[7];
    
    R[6] = M[2];
    R[7] = M[5];
    R[8] = M[8];
}


/*****************************************************************************/
/*	OPERACIONES DE VECTORES						                             */
/*****************************************************************************/
void copyVV(double *v1, double *v2)
{
    v2[0] = v1[0];
    v2[1] = v1[1];
    v2[2] = v1[2];
}

void divVS(double *v, double s, double *r)
{
    r[0] = v[0]/s;
    r[1] = v[1]/s;
    r[2] = v[2]/s;
}

void dotSV(double s, double *v, double *r)
{
    r[0] = s*v[0];
    r[1] = s*v[1];
    r[2] = s*v[2];
}

double dotVV(double *v1, double *v2) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

void antiV(double *v, double *R)
{
    R[0] =  0.0;
    R[1] = -v[2];
    R[2] =  v[1];
    R[3] =  v[2];
    R[4] =  0.0;
    R[5] = -v[0];
    R[6] = -v[1];
    R[7] =  v[0];
    R[8] =  0.0;
}

void cross(double *v1, double *v2, double *r)
{
    r[0] =  v1[1]*v2[2] - v1[2]*v2[1];
    r[1] =  v1[2]*v2[0] - v1[0]*v2[2];
    r[2] =  v1[0]*v2[1] - v1[1]*v2[0];
}

void subVV(double *v1, double *v2, double *r)
{
    r[0] = v1[0] - v2[0];
    r[1] = v1[1] - v2[1];
    r[2] = v1[2] - v2[2];
}

// Suma v1-v2, almacena en r
void addVV(double *v1, double *v2, double *r)
{
    r[0] = v1[0] + v2[0];
    r[1] = v1[1] + v2[1];
    r[2] = v1[2] + v2[2];
}

// Modulo del vector V
double modV(double *v) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/*****************************************************************************/
/*	OPERACIONES MIXTAS						                                 */
/*****************************************************************************/
void dotMV(double *M, double *v, double *r)
{
    r[0] = M[0]*v[0] + M[1]*v[1] + M[2]*v[2];
    r[1] = M[3]*v[0] + M[4]*v[1] + M[5]*v[2];
    r[2] = M[6]*v[0] + M[7]*v[1] + M[8]*v[2];
}
/*****************************************************************************/
/*	OPERACIONES DE ARRAYS						                             */
/*****************************************************************************/
double * newarray(long n) {
    return (double *) malloc(sizeof(double)*n);
}

void copyAA(double *a, double *b, int n)
{
    int i;
    for (i=0; i<n; ++i) {
        b[i] = a[i];
    }
}

void dotSA(double k, double *a, double *r, int n)
{
    int i;
    for (i=0; i<n; ++i) {
        r[i] = k*a[i];
    }
}

void addAA(double *a, double *b, double *r, int n)
{
    int i;
    for (i=0; i<n; ++i) {
        r[i] = a[i] + b[i];
    }
}

void subAA(double *a, double *b, double *r, int n)
{
    int i;
    for (i=0; i<n; ++i) {
        r[i] = a[i] - b[i];
    }
}

double maxA(double *a, int n)
{
    double m = fabs(a[0]);
    int i;
    for (i=1; i<n; ++i) {
        double ma = fabs(a[i]);
        if (ma>m) {
            m = ma;
        }
    }
    return m;
}
