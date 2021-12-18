#include "colors.inc"

background{color Black}

#declare sy=30;
#declare sx=(4/3)*sy;
#declare cz=70;

camera{
    location <0,0,cz>
    look_at <0,0,0>
    right <-4/3,0,0>
    //sky<0,0,1>
}

light_source {<0,0,cz> color White}

#fopen ndt "N.dat" read
#read(ndt,temp)
#declare N=temp;


#declare r=array[N];
#declare rs=array[N];
#declare q=array[N];

#fopen Rsd "rs.dat" read
#declare i=0;
#while (defined(Rsd))
    #read(Rsd,rs1)
    #declare rs[i]=rs1;
    #declare i=i+1;
#end

#fopen qd "q.dat" read
#declare i=0;
#while (defined(qd))
    #read(qd,q1)
    #declare q[i]=q1;
    #declare i=i+1;
#end

#fopen Pts "temp.dat" read
#declare i=0;
#while (defined(Pts))
    #read(Pts,vec1)
    #declare r[i]=vec1;
    #declare i=i+1;
#end

#declare j=0;
#while(j < i)
#if (q[j]=-1)
sphere {
  r[j] , rs[j]
  pigment { color Red }
}
#else
sphere {
  r[j] , rs[j]
  pigment { color Blue }
}
#end
#declare j=j+1;
#end

