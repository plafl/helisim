#ifndef ALGEBRA_H
#define ALGEBRA_H
/*****************************************************************************/
/*  CREACION DE VECTORES Y MATRICES                                          */
/*****************************************************************************/

// Crea un nuevo vector, como un puntero de 3 doubles
double * newvector();

// Crea una nuva matriz, como un puntero de 9 doubles
double * newmatrix();

/*****************************************************************************/
/*	OPERACIONES DE MATRICES						                             */
/*****************************************************************************/

// Escritura de elementos
void setM(double *M, int i, int j, double v);

// Lectura de elementos
double getM(double *M, int i, int j);

// Convierte a identidad
void identity(double *R);

// Multiplica por escalar
void dotSM(double s, double *M, double *R);

// Multiplica dos matrices
void dotMM(double *M, double *N, double *R);

// Suma dos matrices
void addMM(double *M, double *N, double *R);

// Resta dos matrices
void subMM(double *M, double *N, double *R);

// Invierte la matriz
void invertM(double *M, double *R);

// Calcula la transpuesta
void transposeM(double *M, double *R);

/*****************************************************************************/
/*	OPERACIONES DE VECTORES						                             */
/*****************************************************************************/
// Copia v1 a v2
void copyVV(double *v1, double *v2);

// Divide vector por escalar
void divVS(double *v, double s, double *r);

// Multiplica escalar por vector
void dotSV(double s, double *v, double *r);

// Producto escalar de los dos vectores
double dotVV(double *v1, double *v2);

// Matriz antisimetrica de V
void antiV(double *v, double *R);

// Producto vectorial de dos vectores
void cross(double *v1, double *v2, double *r);

// Resta dos vectores
void subVV(double *v1, double *v2, double *r);

// Suma dos vectores
void addVV(double *v1, double *v2, double *r);

/*****************************************************************************/
/*	OPERACIONES MIXTAS						                                 */
/*****************************************************************************/
// Multiplica matriz por vector
void dotMV(double *M, double *v, double *r);


/*****************************************************************************/
/*	ARRAYS                                                                   */
/*****************************************************************************/
// Crea un nuevo array de tamno n
double * newarray(long n);

// Copia a en b, ambos miden n
void copyAA(double *a, double *b, int n);

// Multiplica por el escalar
void dotSA(double k, double *a, double *r, int n);

// Suma arrays
void addAA(double *a, double *b, double *r, int n);

// Resta arrays
void subAA(double *a, double *b, double *r, int n);

// Devuelve el maximo
double maxA(double *a, int n);

#endif
