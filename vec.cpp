vec::vec()
{
    x=0.;
    y=0.;
    z=0.;
}

vec::~vec() {}

vec::vec(const double a,const double b,const double c)
{
    x=a;
    y=b;
    z=c;
}

vec::vec(const double a, const double b)
{
    x=a;
    y=b;
    z=0.;
}

vec::vec(const double a)
{
    x=a;
    y=0.;
    z=0.;
}

vec& vec::operator()(const double a,const double b,const double c)
{
    x=a;
    y=b;
    z=c;
    return (*this);
}

vec& vec::operator()(const double a, const double b)
{
    x=a;
    y=b;
    z=0;
    return (*this);
}

vec& vec::operator()(const double a)
{
    x=a;
    y=0.;
    z=0;
    return (*this);
}

vec& vec::operator()()
{
    x=0.;
    y=0.;
    z=0.;
    return (*this);
}

istream& operator>>(istream &min,vec &a)
{
    min>>a.x>>a.y>>a.z;
    return min;
}

ostream& operator<<(ostream &mout,const vec &a)
{
    mout<<"<"<<a.x<<", "<<a.y<<", "<<a.z<<">";
    return mout;
}

vec& vec::operator=(const vec &a)
{
    x=a.x;
    y=a.y;
    z=a.z;
    return (*this);
}

vec vec::operator+(const vec &a)const
{
    return vec(this->x+a.x,this->y+a.y,this->z+a.z);
}

vec vec::operator-(const vec &a)const
{
    return vec(this->x-a.x,this->y-a.y,this->z-a.z);
}

vec vec::operator*(const double &a)const
{
    return vec(this->x*a,this->y*a,this->z*a);
}

vec operator*(const double &a,const vec &v)
{
    return v*a;
}

vec vec::operator/(const double &a)const
{
    return vec(this->x/a,this->y/a,this->z/a);
}

vec& vec::operator+=(const vec &a) {return (*this=*this+a);}
vec& vec::operator-=(const vec &a) {return (*this=*this-a);}
vec& vec::operator*=(const double &a) {return (*this=*this*a);}
vec& vec::operator/=(const double &a) {return (*this=*this/a);}

bool vec::operator==(const vec &a)const
{
    return (this->x!=a.x ? false : this->y!=a.y ? false : this->z!=a.z ? false : true);
}

bool vec::operator!=(const vec &a)const
{
    return !(*this==a);
}

double vec::dot(const vec &a)const
{
    return ((this->x*a.x)+(this->y*a.y)+(this->z*a.z));
}

double vec::operator*(const vec &a)const
{
    return (this->dot(a));
}

vec vec::cross(const vec &a)const
{
    return vec((this->y*a.z)-(this->z*a.y),(this->z*a.x)-(this->x*a.z),(this->x*a.y)-(this->y*a.x));
}

vec vec::operator%(const vec &a)const
{
    return (this->cross(a));
}

vec& vec::operator*=(const vec &a) {return (*this=this->cross(a));}

double vec::mag2()const
{
    return (this->dot(*this));
}

double vec::mag()const
{
    return (sqrt(this->mag2()));
}

vec vec::hat()const
{
    return (*this/(this->mag()));
}

double vec::angle(const vec &a)const
{
    return acos(((this->dot(a))/((this->mag())*a.mag())));
}

vec vec::rotate(const double &a)const
{
    double sa,ca;
    sa=sin(a);
    ca=cos(a);
    return vec(x*ca+y*sa,-x*sa+y*ca);
}
