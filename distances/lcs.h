#ifndef LCS_H
#define LCS_H

#include "distance.h"

/** Longest Common Subsequence, mide el costo de la subsecuencia comun mas larga
 * entre dos series temporales
 * Soporta series discretas, continuas y series multidimensioanles **/

class LCS: public Distance
{
private:
    //selecciona la tecnica de distance lineal
    // 0=Valores discretos, 1=Valores continuos y para series multidimensionales
    int type_lcs;    
    //banda de aceptacion, si type_lcs==1
    int band;    
    //fraccion de la dimension que se toma como tamanio de la banda
    double fraccion_b;
    //actualizar band de acuerdo a la dimension del par de objects entrantes
    bool actualizar_b;
    //umbral de igualdad, rango <0,1>
    double epsilon;
    

public:
    LCS(int _type_lcs=1, double _fraccion_b=0.25, double _epsilon=0.05);
    virtual ~LCS();
    bool cumple_condicion(Object* p, Object* q, int i, int j);
    
    double d(Object* p,Object* q);
    double lower_bound(Object* ,Object* )
    {
        throw ("** Function not implement **");
    }
    double upper_bound(Object* ,Object* )
    {
        throw ("** Function not implement **");
    }

};

#endif // LINEAL_H
