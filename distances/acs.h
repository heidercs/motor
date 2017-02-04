#ifndef ACS_H
#define ACS_H

#include "distance.h"

/** All common Subsequences, retorna el numero de subsecuencias similares
 * entre dos series temporales
 * Soporta series discretas, continuas y series multidimensioanles **/

class ACS: public Distance
{
private:
    //selecciona la tecnica de distance lineal
    // 0=Valores discretos, 1=Valores continuos y para series multidimensionales
    int type_acs;       
    //umbral de igualdad, rango <0,1>
    double epsilon;
    

public:
    ACS(int _type_acs=1, double _epsilon=0.05);
    virtual ~ACS();
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
