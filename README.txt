1. the following must be installed:
	. c++ compiler (gcc)
	. povray
	. ffmpeg

2. ./run_me.sh


my code is configured for 2D collisions with hard boundary conditions.

it contains 3 different (overloaded) integrators (104-106):
        
	int_dt(r,v,E,B,q,m,dt,N);
        int_dt(r,v,E,B,q,m,k,dt,N);
        int_dt(r,v,E,B,q,m,k,vc,dt,N);


*comment out the 2 you don't want to use


1. only includes an external electric and magnetic field (E,B)
2. + a pairwise coulomb potential (k) {vec qf(int,vec*,double*,double,int);}
3. + a discrete application of biot-savart (vc) {vec bs(int,vec*,vec*,double*,double,double,int);}
