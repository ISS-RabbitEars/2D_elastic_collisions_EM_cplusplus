list::list()
{
    n=mini=maxi=0;
    max=min=0;
    a=new double[n];
    //cout<<"list created"<<endl;
}

list::~list()
{
    delete[] a;
    //cout<<"list deleted"<<endl;
}

void list::append(double x)
{
    if(n==0)
    {
        delete[] a;
        n=1;
        a=new double[1];
        a[0]=x;
        min=max=x;
        mini=maxi=0;
    }
    else
    {
        double *na;
        na=new double[n];
        for(int i=0;i<n;i++) {na[i]=a[i];}
        delete[] a;
        n++;
        a=new double[n];
        for(int i=0;i<n-1;i++) {a[i]=na[i];}
        delete[] na;
        a[n-1]=x;
        if(x>max) {max=x; maxi=n-1;}
        if(x<min) {min=x; mini=n-1;}
    }
}
