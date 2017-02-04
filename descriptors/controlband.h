#ifndef CONTROLBAND_H
#define CONTROLBAND_H

class ControlBand
{

public:
    static bool* get_control_band(double* p,int n, double fact_k=0.1, double fact_x=0.2);
    static double* suavizar(double* p, int n, int k);
};

#endif // CONTROLBAND_H
