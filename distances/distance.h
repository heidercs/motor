#ifndef DISTANCIA_H
#define DISTANCIA_H
#include "../objects/object.h"

class Distance
{
public:
    Distance() {}
    virtual ~Distance() {}
    virtual double d(Object* p, Object* q)=0;
    virtual double lower_bound(Object* p, Object* q)=0;
    virtual double upper_bound(Object* p, Object* q)=0;

    virtual double d(Object* p, Object* q, double _)
    {
        return this->d(p, q);
    }
};

#endif // DISTANCIA_H
