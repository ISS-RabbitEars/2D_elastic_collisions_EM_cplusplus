#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include "ran3.h"
#include "vec.h"
#include "list.h"
using namespace std;
#define PI 3.141592653589793238462643383

void init_pos2D(vec*,vec&,double*,int,int);
void init_vel2D(vec*,double,double,int,int);
void init_mrsq(double*,double*,double*,double,double,int,int);
void int_dt(vec*,vec*,double,int);
void cb2D_pbc(vec*,vec&,int);
void cb2D_pbc(vec&,vec&);
void cc2D_pbc(vec*,vec*,vec&,double*,double*,double,int,bool);
void cb2D_hbc(vec*,vec*,vec&,double*,double,int);
void cb2D_hbc(vec&,vec&,vec&,double&,double);
void cc2D_hbc(vec*,vec*,vec&,double*,double*,double,int,bool);
vec af(int,vec*,double*,double,int);
void int_dt(vec*,vec*,double*,double,double,int);
vec qf(int,vec*,double*,double,int);
void int_dt(vec*,vec*,vec&,vec&,double*,double*,double,int);
void int_dt(vec*,vec*,vec&,vec&,double*,double*,double,double,int);
void int_dt(vec*,vec*,vec&,vec&,double*,double*,double,double,double,int);
vec bs(int,vec*,vec*,double*,double,double,int);

int main()
{
    int seed=5763943,N=1000,fps,ns,fr;
	vec r[N],s,v[N],E,B;
	double rs[N],va,vb,m[N],ma,mb,dt,ar,ay,q[N],k,vc,Emag,Eth;
    bool c;
    ofstream file;
    
    ar=4./3.;
    ay=30;
    s(ar*ay,ay);
    k=1;
    vc=1;
    va=0;
    vb=0;
    ma=1;
    mb=1;
    dt=0.01;
    fps=60;
    ns=3*60;
    fr=fps*ns;
    c=true;
    
    Emag=1;
    Eth=0;
    E=Emag*vec(cos(Eth),sin(Eth));
    B=vec(0,0,1);
    
    init_mrsq(m,rs,q,ma,mb,N,seed);
    init_pos2D(r,s,rs,N,seed);
    init_vel2D(v,va,vb,N,seed);
    
    file.open("pos/N.dat");
    file<<N<<endl;
    file.close();
    
    file.open("pos/nframes.dat");
    file<<fr<<endl;
    file.close();
    
    file.open("pos/rs.dat");
    for(int i=0;i<N;i++)
    {
        file<<rs[i]<<endl;
        if(i<N-1) file<<",";
        file<<endl;
    }
    file.close();
    
    file.open("pos/q.dat");
    for(int i=0;i<N;i++)
    {
        file<<q[i]<<endl;
        if(i<N-1) file<<",";
        file<<endl;
    }
    file.close();
    
    file.open("pos/param.dat");
    file<<"N = "<<N<<endl;
    file<<"s = { "<<s.x<<" , "<<s.y<<" }"<<endl;
    file<<"Aspect Ratio = "<<ar<<endl;
    file<<"k = "<<k<<endl;
    file<<"C = "<<vc<<endl;
    file<<"E = "<<Emag<<endl;
    file<<"Etheta = "<<Eth<<endl;
    file<<"B = "<<B.z<<endl;
    file<<"v0 = { "<<va<<" , "<<vb<<" }"<<endl;
    file<<"m = { "<<ma<<" , "<<mb<<" }"<<endl;
    file.close();
    
    string pts="pos/pos_";
    for(int i=0;i<fr;i++)
    {
        //int_dt(r,v,E,B,q,m,dt,N);
        int_dt(r,v,E,B,q,m,k,dt,N);
        //int_dt(r,v,E,B,q,m,k,vc,dt,N);
        do
        {
            c=false;
            cb2D_hbc(r,v,s,rs,dt,N);
            cc2D_hbc(r,v,s,rs,m,dt,N,c);
        } while(c);
        
        string pu,snum;
        stringstream ssnum;
        ssnum<<i+1;
        snum=ssnum.str();
        ssnum.str("");
        pu=pts+snum;
        file.open(pu.c_str());
        for(int j=0;j<N;j++)
        {
            file<<r[j];
            if(j<N-1) file<<",";
            file<<endl;
        }
        file.close();
    }
    
	return 0;
}

void init_pos2D(vec* r, vec& s, double* rs, int N, int seed)
{
    vec rt,dr;
    bool col;
    double ds;
    
    r[0]((s.x-rs[0])*(1-2*ran3(&seed)),(s.y-rs[0])*(1-2*ran3(&seed)));
    for(int i=1;i<N;i++)
    {
        col=false;
        do
        {
            rt((s.x-rs[i])*(1-2*ran3(&seed)),(s.y-rs[i])*(1-2*ran3(&seed)));
            int j=0;
            do
            {
                col=false;
                dr=rt-r[j];
                ds=rs[i]+rs[j];
                if(dr.mag()<ds) {col=true;}
                j++;
            }while(j<i&&!col);
        }while(col);
        r[i]=rt;
    }
}

void init_vel2D(vec* v,double a,double b, int N, int seed)
{
    double alpha,theta,vi;
    
    alpha=b-a;
    for(int i=0;i<N;i++)
    {
        theta=2*PI*ran3(&seed);
        vi=alpha*ran3(&seed)+a;
        v[i](vi*cos(theta),vi*sin(theta));
    }
}

void init_mrsq(double* m, double* rs, double* q, double a, double b, int N, int seed)
{
    double alpha=b-a;
    for(int i=0;i<N;i++)
    {
        m[i]=alpha*ran3(&seed)+a;
        rs[i]=m[i]/(2*b);
        if(ran3(&seed)<0.5)
            q[i]=-1;
        else
            q[i]=1;
    }
}

void int_dt(vec* r, vec* v, double dt, int N)
{
    for(int i=0;i<N;i++)
    {
        r[i]+=v[i]*dt;
    }
}

void cb2D_pbc(vec* r, vec& s, int N)
{
    for(int i=0;i<N;i++)
    {
        cb2D_pbc(r[i],s);
    }
}

void cb2D_pbc(vec& r, vec& s)
{
    if(r.x<-s.x) {r.x+=2*s.x;}
    if(r.x>=s.x) {r.x-=2*s.x;}
    if(r.y<-s.y) {r.y+=2*s.y;}
    if(r.y>=s.y) {r.y-=2*s.y;}
}

void cc2D_pbc(vec* r, vec* v, vec& s, double* rs, double* m, double t, int n, bool c)
{
    for(int i=0;i<n;i++)
    {
        list drc,drj,drlmi;
        bool col=false;
        for(int j=0;j<n;j++)
        {
            if(j!=i)
            {
                vec drt,drp1,drp2,temp;
                list drl;
                
                drt=r[j]-r[i];
                drl.append(drt.mag());
                drp1=drt+2*s;
                drp2=drt-2*s;
                temp(drp1.x,drt.y);
                drl.append(temp.mag());
                temp(drp2.x,drt.y);
                drl.append(temp.mag());
                temp(drt.x,drp1.y);
                drl.append(temp.mag());
                temp(drt.x,drp2.y);
                drl.append(temp.mag());
                drl.append(drp1.mag());
                temp(drp1.x,drp2.y);
                drl.append(temp.mag());
                temp(drp2.x,drp1.y);
                drl.append(temp.mag());
                drl.append(drp2.mag());
                if(drl.min<=rs[i]+rs[j])
                {
                    drc.append(drl.min);
                    drj.append(j);
                    drlmi.append(drl.mini+1);
                    col=true;
                    c=true;
                }
            }
        }
        if(col)
        {
            int j,id;
            vec dr,dv,drh,vri,vrj;
            double qa,qb,qc,dt1,dt2,tr,sm,dm,vifx,vjfx;
            
            j=drj.a[drc.mini];
            id=drlmi.a[drc.mini];
            r[i]-=v[i]*t;
            r[j]-=v[j]*t;
            switch(id)
            {
                case 1: dr=r[j]-r[i]; break;
                case 2: dr=r[j]-r[i]+vec(2*s.x,0); break;
                case 3: dr=r[j]-r[i]+vec(-2*s.x,0); break;
                case 4: dr=r[j]-r[i]+vec(0,2*s.y); break;
                case 5: dr=r[j]-r[i]+vec(0,-2*s.y); break;
                case 6: dr=r[j]-r[i]+2*s; break;
                case 7: dr=r[j]-r[i]+vec(2*s.x,-2*s.y); break;
                case 8: dr=r[j]-r[i]+vec(-2*s.x,2*s.y); break;
                case 9: dr=r[j]-r[i]-2*s; break;
            }
            dv=v[j]-v[i];
            qa=dv.mag2();
            qb=2*(dr*dv);
            qc=dr.mag2()-pow(rs[i]+rs[j],2);
            dt1=(-qb-sqrt(qb*qb-4*qa*qc))/(2*qa);
            dt2=t-dt1;
            r[i]+=v[i]*dt1;
            r[j]+=v[j]*dt1;
            switch(id)
            {
                case 1: dr=r[j]-r[i]; break;
                case 2: dr=r[j]-r[i]+vec(2*s.x,0); break;
                case 3: dr=r[j]-r[i]+vec(-2*s.x,0); break;
                case 4: dr=r[j]-r[i]+vec(0,2*s.y); break;
                case 5: dr=r[j]-r[i]+vec(0,-2*s.y); break;
                case 6: dr=r[j]-r[i]+2*s; break;
                case 7: dr=r[j]-r[i]+vec(2*s.x,-2*s.y); break;
                case 8: dr=r[j]-r[i]+vec(-2*s.x,2*s.y); break;
                case 9: dr=r[j]-r[i]-2*s; break;
            }
            drh=dr.hat();
            tr=acos(drh.x);
            if(drh.y<0) {tr=-tr;}
            vri=v[i].rotate(tr);
            vrj=v[j].rotate(tr);
            sm=m[i]+m[j];
            dm=m[j]-m[i];
            vifx=(-dm/sm)*vri.x+(2*m[j]/sm)*vrj.x;
            vjfx=(2*m[i]/sm)*vri.x+(dm/sm)*vrj.x;
            vri.x=vifx;
            vrj.x=vjfx;
            v[i]=vri.rotate(-tr);
            v[j]=vrj.rotate(-tr);
            r[i]+=v[i]*dt2;
            r[j]+=v[j]*dt2;
            cb2D_pbc(r[i],s);
            cb2D_pbc(r[j],s);
        }
    }
}

void cb2D_hbc(vec* r, vec* v, vec& s, double* rs, double t, int n)
{
    double dt1,dt2;
    for(int i=0;i<n;i++)
    {
        cb2D_hbc(r[i],v[i],s,rs[i],t);
    }
}

void cb2D_hbc(vec& r, vec& v, vec& s, double& rs, double t)
{
    double dt1,dt2;
    if(r.x-rs<-s.x)
    {
        r-=v*t;
        dt1=(-s.x+rs-r.x)/v.x;
        dt2=t-dt1;
        r+=v*dt1;
        v.x=-v.x;
        r+=v*dt2;
    }
    if(r.x+rs>=s.x)
    {
        r-=v*t;
        dt1=(s.x-rs-r.x)/v.x;
        dt2=t-dt1;
        r+=v*dt1;
        v.x=-v.x;
        r+=v*dt2;
    }
    if(r.y-rs<-s.y)
    {
        r-=v*t;
        dt1=(-s.y+rs-r.y)/v.y;
        dt2=t-dt1;
        r+=v*dt1;
        v.y=-v.y;
        r+=v*dt2;
    }
    if(r.y+rs>=s.y)
    {
        r-=v*t;
        dt1=(s.y-rs-r.y)/v.y;
        dt2=t-dt1;
        r+=v*dt1;
        v.y=-v.y;
        r+=v*dt2;
    }
}

void cc2D_hbc(vec* r, vec* v, vec& s, double* rs, double* m, double t, int n, bool c)
{
    for(int i=0;i<n;i++)
    {
        list drc,drj;
        bool col=false;
        for(int j=0;j<n;j++)
        {
            if(j!=i)
            {
                vec drt;
                double drtm,rss;
                
                drt=r[j]-r[i];
                drtm=drt.mag();
                rss=rs[i]+rs[j];
                if(drtm<=rss)
                {
                    drc.append(drtm/rss);
                    drj.append(j);
                    col=true;
                    c=true;
                }
            }
        }
        if(col)
        {
            int j=drj.a[drc.mini];
            vec dr,dv,drh,vri,vrj;
            double qa,qb,qc,dt1,dt2,tr,sm,dm,vifx,vjfx;
            
            r[i]-=v[i]*t;
            r[j]-=v[j]*t;
            dr=r[j]-r[i];
            dv=v[j]-v[i];
            qa=dv.mag2();
            qb=2*(dr*dv);
            qc=dr.mag2()-pow(rs[i]+rs[j],2);
            dt1=(-qb-sqrt(qb*qb-4*qa*qc))/(2*qa);
            dt2=t-dt1;
            r[i]+=v[i]*dt1;
            r[j]+=v[j]*dt1;
            drh=dr.hat();
            tr=acos(drh.x);
            if(drh.y<0) {tr=-tr;}
            vri=v[i].rotate(tr);
            vrj=v[j].rotate(tr);
            sm=m[i]+m[j];
            dm=m[j]-m[i];
            vifx=(-dm/sm)*vri.x+(2*m[j]/sm)*vrj.x;
            vjfx=(2*m[i]/sm)*vri.x+(dm/sm)*vrj.x;
            vri.x=vifx;
            vrj.x=vjfx;
            v[i]=vri.rotate(-tr);
            v[j]=vrj.rotate(-tr);
            r[i]+=v[i]*dt2;
            r[j]+=v[j]*dt2;
            cb2D_hbc(r[i],v[i],s,rs[i],t);
            cb2D_hbc(r[j],v[j],s,rs[j],t);
        }
    }
}

vec af(int i, vec* r, double* m, double G, int n)
{
    vec dr,sum;
    for(int j=0;j<n;j++)
    {
        if(j!=i)
        {
            dr=r[j]-r[i];
            sum+=(G*m[j]/dr.mag2())*dr.hat();
        }
    }
    return sum;
}

void int_dt(vec* r, vec* v, double* m, double G, double dt, int n)
{
    for(int i=0;i<n;i++)
    {
        v[i]+=af(i,r,m,G,n)*dt;
        r[i]+=v[i]*dt;
    }
}

vec qf(int i, vec* r, double* q, double k, int n)
{
    vec dr,sum;
    for(int j=0;j<n;j++)
    {
        if(j!=i)
        {
            dr=r[j]-r[i];
            sum+=-(k*q[i]*q[j]/dr.mag2())*dr.hat();
        }
    }
    return sum;
}

void int_dt(vec* r, vec* v, vec& E, vec& B, double* q, double* m, double t, int n)
{
    for(int i=0;i<n;i++)
    {
        v[i]+=(q[i]/m[i])*(E+(v[i]%B))*t;
        r[i]+=v[i]*t;
    }
}

void int_dt(vec* r, vec* v, vec& E, vec& B, double* q, double* m, double k, double dt, int n)
{
    for(int i=0;i<n;i++)
    {
        v[i]+=((qf(i,r,q,k,n)+q[i]*(E+(v[i]%B)))/m[i])*dt;
        r[i]+=v[i]*dt;
    }
}

void int_dt(vec* r, vec* v, vec& E, vec& B, double* q, double* m, double k, double c, double dt, int n)
{
    vec bq;
    for(int i=0;i<n;i++)
    {
        bq=bs(i,r,v,q,k,c,n); // augmented for 2D (use z component only)
        v[i]+=((qf(i,r,q,k,n)+q[i]*(E+(v[i]%(B+vec(0,0,bq.z)))))/m[i])*dt;
        r[i]+=v[i]*dt;
    }
}

vec bs(int i, vec* r, vec* v, double* q, double k, double c, int n)
{
    vec dr,sum;
    for(int j=0;j<n;j++)
    {
        if(j!=i)
        {
            dr=r[i]-r[j]; //to get proper direction for dB
            sum+=(k/(c*c*dr.mag2()))*q[j]*(v[j]%dr.hat());
        }
    }
    return sum;
}

#undef PI
