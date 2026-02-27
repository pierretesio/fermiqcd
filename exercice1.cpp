// PROGRAM1: Use and test of random generator
//
#include "MyHeaders.h"          //voir Préliminaires
ofstream outfile("random.dat"); //open stream for output
//
float myrandom(float a,float b)      //random number in [a,b]
{
    float range=pow(2.,31);
    return a+(b-a)*random()/range;
}

//
//
int main()
{
float x,sum=0.,mean;
int i,hitnb=100000;                //make hitnb calls
cout<<"create "<< hitnb<<" rand. numb. in [0,1]\ 
        and compute mean"<<endl;

for (i=1;i<=hitnb;i++)
{
    x=myrandom(0.,1.);
    outfile<<i<<" "<<x<<endl;      //write file for histogram
    sum=sum+x;
    mean=sum/i;
    if(i%(hitnb/100)==0)
    {
        cout<<"i= "<<i<<" mean= "<<mean<<endl;
    }
}

outfile.close();
// calcul des intégrales...
//end of PROGRAM1
}