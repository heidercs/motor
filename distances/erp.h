#ifndef ERP_H
#define ERP_H

#include "distance.h"

/** Edit distance with Real Penalty, une lo mejor de EDR y DTW
  * Soporta series continuas unidimensionales y multidimensioanles **/

class ERP: public Distance
{
private:    
    // brecha de compensacion en valores faltantes para el matching, rango [0, inf>
    double gap;
    

public:
    ERP(double _gap=0);
    virtual ~ERP();
    double dist(Object* p, Object* q, int i, int j);
    double dist(Object* p, int i, int gap);
    

    double d(Object* p,Object* q);
    double lower_bound(Object* ,Object* );
    double upper_bound(Object* ,Object* )
    {
        throw ("** Function not implement **");
    }

};

#endif // LINEAL_H
