#ifndef VEC_H
#define VEC_H
#include <iostream>
#include <cmath>
using namespace std;

class vec
{
public:
    
    double x,y,z;
    vec();
    ~vec();
    vec(const double, const double, const double);
    vec(const double, const double);
    vec(const double);
    vec& operator()(const double,const double,const double);
    vec& operator()(const double, const double);
    vec& operator()(const double);
    vec& operator()();
    friend istream& operator>>(istream&,vec&);
    friend ostream& operator<<(ostream&,const vec&);
    vec& operator=(const vec&);
    vec operator+(const vec&)const;
    vec operator-(const vec&)const;
    vec operator*(const double&)const;
    friend vec operator*(const double&,const vec&);
    vec operator/(const double&)const;
    vec& operator+=(const vec&);
    vec& operator-=(const vec&);
    vec& operator*=(const double&);
    vec& operator/=(const double&);
    bool operator==(const vec&)const;
    bool operator!=(const vec&)const;
    double dot(const vec&)const;
    double operator*(const vec&)const; //(arbitrary) defined as dot product
    vec cross(const vec&)const;
    vec operator%(const vec&)const; //(arbitrary) defined as cross product
    vec& operator*=(const vec&); //necessarily infers cross product definition
    double mag2()const;
    double mag()const;
    vec hat()const;
    double angle(const vec&)const; //(radians) [0,+Pi]
    vec rotate(const double&)const;
};

#include "vec.cpp"
#endif
