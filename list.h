#ifndef list_h
#define list_h
#include <iostream>
using namespace std;

class list
{
public:
    int n,mini,maxi;
    double *a,min,max;
    list();
    ~list();
    void append(double);
};

#include "list.cpp"
#endif
