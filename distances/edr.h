#ifndef EDR_H
#define EDR_H

#include "distance.h"

/** Edit Distance Real Sequence, mide el costo de transformacions de una serie a otra
  * Soporta series discretas, continuas y series multidimensioanles **/

class EDR: public Distance
{
private:
    //selecciona la tecnica de distance lineal
    // 0=Valores discretos, 1=Valores continuos y para series multidimensionales
    int type_edr;       
    //umbral de igualdad, rango <0,1>
    double epsilon;
    

public:
    EDR(int _type_edr=1, double _epsilon=0.05);
    virtual ~EDR();
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
