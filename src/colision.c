/*
 * Compilar con gcc -fPIC -shared algebra.c colision.c -o c_colision
 */

//#define DEBUG_MSG

#include <math.h>
#include "algebra.h"
#include "colision.h"

struct Tren * malloc_tren(long nP)
{
    struct Tren * T = malloc(sizeof(struct Tren));
    
    T->x0 = newarray(ARRAY_SIZE);
    T->x1 = newarray(ARRAY_SIZE);
    T->x2 = newarray(ARRAY_SIZE);

    T->F0 = newarray(ARRAY_SIZE);
    T->F1 = newarray(ARRAY_SIZE);
    T->F2 = newarray(ARRAY_SIZE);
    T->Fi = newarray(ARRAY_SIZE);

    T->P = malloc(nP*sizeof(struct Punto *));
    T->nP = nP;
#ifdef DEBUG_MSG
    printf("malloc_tren: memoria inicializada para %p\n", T);
#endif
    return T;
}

void free_tren(struct Tren * T)
{
    free(T->x0); free(T->x1); 
    free(T->x2); 
    free(T->F0); free(T->F1);
    free(T->F2); free(T->Fi);
    free(T->P);
#ifdef DEBUG_MSG
    printf("malloc_tren: memoria liberada para %p", T);
#endif
}

void step_tren(struct Tren * T, double t)
{
#ifdef DEBUG_MSG
    printf("step_tren llamada desde %p, t = %f\n", T, t);
#endif
    // Paso minimo de integracion
    double dtmin = 1e-8;
    // Paso maximo de integracion
    double dtmax = 0.01;

    // Arrays temporales
    double * a1 = newarray(ARRAY_SIZE);
    double * a2 = newarray(ARRAY_SIZE);
    double * a3 = newarray(ARRAY_SIZE);
    
    if (T->colision == 0 && t>T->t) 
    {
#ifdef DEBUG_MSG
        printf("step_tren: iniciando integrador\n");
#endif
        double tf = T->t + dtmin;
        // Condiciones iniciales
        T->xi = T->x;
        T->ti = T->t;
        F(T->xi, T->ti, T->Fi, T);

        copyAA(T->xi, T->x1, ARRAY_SIZE);
        copyAA(T->Fi, T->F1, ARRAY_SIZE);
        T->t1 = T->ti;
        T->dt1 = dtmin;
        
        // Metodo de Euler: x0 = x1 + dt*F1
        dotSA( T->dt1, T->F1, a1,    ARRAY_SIZE );
        addAA( a1,     T->x1, T->x0, ARRAY_SIZE );

        // Guardamos otros resultados
        T->t = T->t0 = tf;
        F(T->x0, T->t0, T->F0, T);
        copyAA(T->x0, T->x, ARRAY_SIZE);
//        T->colision = 1;
#ifdef DEBUG_MSG
        printf("step_tren: integrador iniciado\n");
#endif
    } 
    if (t>(T->t0 + dtmin))
    {
#ifdef DEBUG_MSG
        printf("step_tren: dando un paso\n");
#endif
        double * xp = newarray(ARRAY_SIZE);
        double * Fp = newarray(ARRAY_SIZE);

        // Se ha producido una colision en algun punto del intervalo?
        long colision = 0;

        // guardamos la memoria de x2 y F2
        T->xa = T->x2;
        T->Fa = T->F2;
        
        // desplazamos las magnitudes
        T->x2 = T->x1;
        T->F2 = T->F1;
        T->t2 = T->t1;

        T->x1 = T->x0;
        T->F1 = T->F0;
        T->t1 = T->t0;

        T->dt2 = T->dt1;
        T->dt1 = t - T->t0;
        
        int STOP = 0;
        // STOP si Tn<T0 y t0 > t - dtmin
        while (!STOP) 
        {
            T->t0  = T->t1 + T->dt1;
            // Predictor: xp, Fp
#ifdef DEBUG_MSG
            printf("Predictor\n");
#endif
            double bp1 = (2*T->dt2 + T->dt1)/(2*T->dt2);
            double bp2 = -T->dt1/(2*T->dt2);
            
            dotSA(T->dt1*bp1, T->F1, a1,  ARRAY_SIZE);
            dotSA(T->dt1*bp2, T->F2, a2,  ARRAY_SIZE);
            addAA(a1,     a2,    a3,      ARRAY_SIZE);
            addAA(T->x1,  a3,    xp,      ARRAY_SIZE); 
           
            F(xp, T->t0, Fp, T);
            if (T->colision == 1) {
                colision = 1;
            }
                
            // Corrector: x0, F0
#ifdef DEBUG_MSG
            printf("Corrector\n");
#endif
            double bc0 = 0.5;
            double bc1 = 0.5;
            
            T->x0 = T->xa;
            T->F0 = T->Fa;

            dotSA(T->dt1*bc0, Fp,    a1,    ARRAY_SIZE);
            dotSA(T->dt1*bc1, T->F1, a2,    ARRAY_SIZE);
            addAA(a1,         a2,    a3,    ARRAY_SIZE);
            addAA(T->x1,      a3,    T->x0, ARRAY_SIZE); 
           
            F(T->x0, T->t0, T->F0, T);
            if (T->colision == 1) {
                colision = 1;
            }
                
            // Error: Tn = |Cc/(Cp + Cc)*(x0p - x00)|
            double Cp = (2*T->dt1 + 3*T->dt2)/(12*T->dt1);
            double Cc = -1./12.;
            
            subAA(xp, T->x0, a1, ARRAY_SIZE);
            double Tn = fabs((Cc/(Cp + Cc)*maxA(a1, ARRAY_SIZE))); 
            
            if (Tn>1.1*T->T0) 
            {
                // Repetimos el paso
                T->t0 = T->t1;
                // Nuevo paso de integracion
                if (Tn>0) {
                    T->dt1 *= pow((T->T0/Tn),(1./3.));
                    if (T->dt1>dtmax) {
                        T->dt1 = dtmax;
                    }
                }
            } else if (T->t0<t-dtmin) 
            {
                // No hace falta repetir paso, pero
                // hay que seguir dandolos --> STOP = 0

                // Avanzamos las variables
                T->xa = T->x2;
                T->Fa = T->F2;

                T->x2 = T->x1;
                T->F2 = T->F1;
                T->t2 = T->t1;

                T->x1 = T->x0;
                T->F1 = T->F0;
                T->t1 = T->t0;

                T->dt2 = T->dt1;
                // Nuevo paso de integracion
                if (Tn>0) {
                    T->dt1 *= pow((T->T0/Tn),(1./3.));
                    if (T->dt1>dtmax) {
                        T->dt1 = dtmax;
                    }
                }
            } else {
                STOP = 1;
            }
#ifdef DEBUG_MSG
             printf("step_tren: paso dado con Tn, T0 = %f, %f\n", Tn, T->T0);
#endif
        } 

        // Estamos en colision si en alguna parte del intervalo ha
        // habido colision
        T->colision = colision;
        copyAA(T->x0, T->x, ARRAY_SIZE);
        T->t = T->t0;
        free(xp); free(Fp);

#ifdef DEBUG_MSG
        printf("step_tren: paso dado, integrador finalizado\n");
#endif
    }
    free(a1); free(a2); free(a3);
}
        

/*
 * De Tren necesita:
 *      FC, MC, FB, MB
 *      puntos
 * 
 * En Tren coloca:
 *      colision
 */
void F(double *x, double t, double *y, struct Tren * T)
{
#ifdef DEBUG_MSG
    printf("Funcion F llamada para t=%f\n", t);
#endif
    // Matrices de cambio de ejes cuerpo - ejes colision
    double * LBC = newmatrix();
    double * LCB = newmatrix();
        
    // Fuerzas de reaccion
    double *  FR = newvector();
    double *  MR = newvector();
        
    // Fuerzas totales
    double *  FT = newvector();
    double *  MT = newvector();
        
    // Posicion del punto en ejes colision
    double *  dC = newvector();

    // Desplzamiento y velocidad de desplazamiento
    double *  eC = newvector();
    double *  vC = newvector();
        
    // Posicion, velocidad y vel. angular CM
    double * rC = newvector();
    double * vB = newvector();
    double * wB = newvector();
        
    // Vectores temporales
    double * v1 = newvector();
    double * v2 = newvector();
    double * v3 = newvector();
    double * v4 = newvector();
    double * v5 = newvector();
        
        
    // Componentes de la velocidad en ejes cuerpo
    double u = vB[0] = x[0];
    double v = vB[1] = x[1];
    double w = vB[2] = x[2];
        
    // Componentes de la velocidad angular en ejes cuerpo
    double p = wB[0] = x[3];
    double q = wB[1] = x[4];
    double r = wB[2] = x[5];

    // Componentes del cuaternio
    double q0 = x[6];
    double q1 = x[7];
    double q2 = x[8];
    double q3 = x[9];

    // Componentes del vector de posicion
    double xi = rC[0] = x[10];
    double yi = rC[1] = x[11];
    double zi = rC[2] = x[12];

    setM(LCB, 0, 0, 1 - 2*( q2*q2 + q3*q3));
    setM(LCB, 0, 1,     2*(-q0*q3 + q1*q2));
    setM(LCB, 0, 2,     2*( q0*q2 + q1*q3));
    setM(LCB, 1, 0,     2*( q0*q3 + q1*q2));
    setM(LCB, 1, 1, 1 - 2*( q1*q1 + q3*q3));
    setM(LCB, 1, 2,     2*(-q0*q1 + q2*q3));
    setM(LCB, 2, 0,     2*(-q0*q2 + q1*q3));
    setM(LCB, 2, 1,     2*( q0*q1 + q2*q3));
    setM(LCB, 2, 2, 1 - 2*( q1*q1 + q2*q2));

    transposeM(LCB, LBC);
#ifdef DEBUG_MSG
    printf("F: matriz LBC y LCB calculadas\n");
#endif
    FR[0] = FR[1] = FR[2] = 0.0;
    MR[0] = MR[1] = MR[2] = 0.0;

    // Algun punto contacta?
    long contacto = 0;
    // Margen para el contacto
    double margen = 0.1;

    // Para cada punto del tren
    long iP;
    for (iP = 0; iP<(T->nP); ++iP)
    {
#ifdef DEBUG_MSG
        printf("F: punto %ld de %ld\n", iP + 1, T->nP);
#endif
        struct Punto *P = T->P[iP];

        // Calculamos dC = r + LCB*P.x
        dotMV(LCB, P->x, v1);
        addVV(v1,  rC,  dC);

        if (dC[2]>0)
        {
#ifdef DEBUG_MSG
            printf("F: punto sin contacto\n");
#endif
            // Si no hay colision quitamos el punto de contacto
            P->D = 0;
            // De todas formas miramos el margen
//            if (dC[2]<margen) {
//                contacto = 1;
//            }
            // y pasamos al siguiente punto
            continue;
        } else  
        {
            contacto = 1;
            // Hay colision
            if (P->D == 0)
            {
#ifdef DEBUG_MSG
            printf("F: creando contacto\n");
#endif
                // Pero no habia punto de contacto, lo creamos
                P->y[0] = dC[0];
                P->y[1] = dC[1];
                P->y[2] = 0;
                // Recordamos que ya esta creado
                P->D = 1;
            }
            // Calculamos eC = dC - P.y 
            subVV(dC,  P->y, eC);
            
            // Calculamos vC = LCB*(vB + wB^P.x)
            cross(wB,  P->x, v1);
            addVV(vB,  v1,  v2);
            dotMV(LCB, v2,  vC);
             
            double es = eC[0];
            double et = eC[1];
            double en = eC[2];

            double vs = vC[0];
            double vt = vC[1];
            double vn = vC[2];

            // Fuerza normal
            double N1 = -P->a*en;
            double N2;
            if (vn<0) {
                N2 = -P->b*vn;
            } else {
                N2 = 0;
            }
            double N = N1 + N2;

            // Fuerzas tangenciales
            double S = -(P->c*es + P->d*vs);
            double T = -(P->e*et + P->f*vt);
#ifdef DEBUG_MSG
            printf("c*es = %f, d*vs = %f\n", P->c*es, P->d*vs);
            printf("e*et = %f, f*vt = %f\n", P->e*et, P->f*vt);
#endif

            // Siempre se oponen a la velocidad
            int mismo_signo(double a, double b){
                return ((a<0 && b<0) || (a>=0 && b>=0));
            }
//            if (mismo_signo(S, vs)){
//                S = 0;
//            }
//            if (mismo_signo(T, vt)){
//                T = 0;
//            }
            // Comprobamos que la fuerza de rozamiento esta dentro
            // de los limietes dados por la normal y el coeficiente
            // de rozamiento
            double R = sqrt(S*S + T*T);
            double Rmax = P->n*N;
            // Se pueden producir errores numericos
            if (Rmax<0) {
                Rmax = 0;
            }

#ifdef DEBUG_MSG
            printf("F: R = %f, Rmax = %f, n = %f, N=%f\n", R, Rmax, P->n, N);
#endif
            if (R>Rmax)
            {
                // Hay deslizamiento, desplazamos el punto de contacto
                double k;
                if (es!=0 || et!=0) {
                    k = Rmax/sqrt(P->c*P->c*es*es + P->e*P->e*et*et);
                } else {
                    k = 1;
                }
                v1[0] = (1-k)*es;
                v1[1] = (1-k)*et;
                v1[2] = 0;
                addVV (P->y, v1,  v2 );
                copyVV(v2,   P->y    );
                // Trucamos el rozamiento
                S = Rmax*S/R;
                T = Rmax*T/R;
            }
            // Sumamos las fuerzas y momentos de reaccion
            v1[0] = S;
            v1[1] = T;
            v1[2] = N;
#ifdef DEBUG_MSG
            printf("F: S,T,N = %f,%f,%f\n", S, T, N);
#endif
            dotMV (LBC,  v1, v2);
            cross (P->x, v2, v3);
            addVV (FR,   v2, v4);
            addVV (MR,   v3, v5);
            copyVV(v4,   FR    );
            copyVV(v5,   MR    );
        }    
    } // Fin para el punto
    // Ya hemos calculado todas las reacciones, calculamos
    // las fuerzas y momentos totales FT, MT
    dotMV (LBC, T->FC, v1);
    addVV (v1,  FR,    v2);
    addVV (v2,  T->FB, v3);
    copyVV(v3,  FT);
            
    dotMV (LBC, T->MC, v1);
    addVV (v1,  MR,    v2);
    addVV (v2,  T->MB, v3);
    copyVV(v3,  MT);

#ifdef DEBUG_MSG
    printf("F: FR = %f, %f, %f\n", FR[0], FR[1], FR[2]);
    printf("F: MR = %f, %f, %f\n", MR[0], MR[1], MR[2]);
    printf("F: FT = %f, %f, %f\n", FT[0], FT[1], FT[2]);
    printf("F: MT = %f, %f, %f\n", MT[0], MT[1], MT[2]);
#endif
    // Unos alias
    double M   = T->M;
    double Ixx = T->Ixx;
    double Iyy = T->Iyy;
    double Izz = T->Izz;
    double Ixz = T->Ixz;
            
    // Aceleraciones de inercia
    double axI = -( w*q - v*r );
    double ayI = -( u*r - w*p );
    double azI = -( v*p - u*q ); 
    // Momentos de inercia
    double LI =  ( Iyy - Izz )*q*r + Ixz*p*q;
    double MI =  ( Izz - Ixx )*r*p + Ixz*( r*r - p*p );
    double NI =  ( Ixx - Iyy )*p*q + Ixz*q*r;
            
    double l = LI + MT[0];
    double m = MI + MT[1];
    double n = NI + MT[2];
    double k = Ixx*Izz - Ixz*Ixz;

    // Ecuaciones de fuerzas y momentos
    y[0] = axI + FT[0]/M;
    y[1] = ayI + FT[1]/M;
    y[2] = azI + FT[2]/M;
    y[3] = (Izz*l + Ixz*n)/k;
    y[4] = m/Iyy;
    y[5] = (Ixz*l + Ixx*n)/k;

    dotMV(LCB, wB, v1);
    double pi = v1[0];
    double qi = v1[1];
    double ri = v1[2];
            
    // Relacion cinematica para los cuaternios
    y[6] = -0.5*(  q1*pi + q2*qi + q3*ri );
    y[7] = -0.5*( -q0*pi - q3*qi + q2*ri );
    y[8] = -0.5*(  q3*pi - q0*qi - q1*ri ); 
    y[9] = -0.5*( -q2*pi + q1*qi - q0*ri );

    dotMV(LCB, vB, v1);
    double ui = v1[0];
    double vi = v1[1];
    double wi = v1[2];

    // Relacion cinematica para la posicion del CM
    y[10] = ui;
    y[11] = vi;
    y[12] = wi;
        
    T->colision = contacto;
#ifdef DEBUG_MSG
    printf("F: Fin\n");
#endif
    // Liberamos memoria usada
    free(v1);  free(v2);  free(v3); free(v4); free(v5);
    free(FR);  free(MR);  free(FT); free(MT);
    free(dC);  free(eC);  free(vC);
    free(rC);  free(vB);  free(wB);
    free(LBC); free(LCB);
} // Fin de F

