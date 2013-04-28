

#include <iostream>    //header for basic io
#include <cmath>       //header for math functions
#include <fstream>     //header for file io
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


using namespace std;

//This is the main driver of a MD code
class MODEL
{
   public:
   MODEL();
   ~MODEL();
   //model description
   int npart;               //number of particles
   double (*x)[3];      //position of particles (in Angstroms)
   double (*v)[3];      //velocity (in A/ps)
   double (*f)[3];       //force (in kcal/mol/A)
   double (*fm)[3];     //forces from previous time step
   double (*m);         //mass (in g/mol)
   double Do;      //Lenndard Jones parameter in kcal/mol
   double Ro;      //Lenndard Jones parameter in Angstrom
   double p0;      //p0= 12*Do*pow(Ro,12);  in unit kcal/mol*A12
   double p1;      //p1= 12*Do*pow(Ro,6);   in unit kcal/mol*A6
   double mass;    //particle mass (in g/mol)
   double cell[3]; //unit cell length (A)
   bool period;    //periodity : 1
   bool gen_gro;  //generate .gro file(for vmd)

   //MD parameters
   double dt;             //time step (in ps)
   double t;               //current time (in ps)
   double Tset;         //temperature used in velocity initialization
   int nstep;              //total steps to run
   int istep;               //current step
   int nsamp;            //sampling frequency
   double Rcut;          //cut off distance for vdw interaction
   double rc2;           //Rcut*Rcut


   //model properties
   double T;                 //temperature (in K)
   double P;                 //pressure    (in GPa)
   double V;                 //volume      (in A3)
   int df;                       //degree of freedom
   double Ek;               //kinetic energy in kcal/mol
   double Ep;               //potential energy in kcal/mol
   double Etot;             //Ek+Ep
   double Ptail; //tail correction for pressure (GPa)
   double Eptail; //tail correction for P.E. (kcal/mol)
   double strs[6]; //stress tensor in kcal/mol-A3

   //class functions
   int init();                    //initialization of variables
   int force();                //force calculation
   int integrate();          //verlet integration
   int sample();             //calculation of properties
   int cal_T();                //calculation of system temperature;
   int cal_P();                //calculation of system pressure;
   double myrand01(); //generate a random number between 0 and 1 [0,1]
   int dump_trj();    //output trajectory data to a file
   void output();
   char trjfname[1024];     //filename of the MD trajectory
   char logfname[1024];     //filename for the MD log

   ofstream outf;     //file stream of the trajectory
   ofstream logf;     //file streeam of log
   ifstream input;
};


int main()
{
                
        MODEL model;                                 //declare a system
        
        model.init();                                       //initialization
        

        while(model.istep<model.nstep) {  //MD loop
            model.integrate();                         //integrate equation of motion
        }
        
        model.output();
    
    
    return 0;
}

int MODEL::init()
{
    //Simulation Parameters
    
    cin >>nstep; //steps to run
    cin >>npart; //number of particles
    cin >>Tset; //target temperature
    cin >>dt; // time step in ps    

    istep=0;                                  //current step
    nsamp=100;                            //save data every nsamp steps
    Do= 0.302;    // in kcal/mol, for Ar
    Ro= 4.237;    // in Angstrom, for Ar
    p0= 12*Do*pow(Ro,12); //in unit kcal/mol*A12
    p1= 12*Do*pow(Ro,6); //in unit kcal/mol*A6
    mass=39.948;        //MW of Argon in g/mol
    cell[0] = cell[1] = cell[2] = pow(npart*60,1.0/3.0);   //length of the unit cell
    period = 1;        //flag for periodicity
    gen_gro = 0;    //flag for generating .gro file
    Rcut = 2.5*Ro;     //cut off distance for vdw interaction
    //(It must be larger than 2.5*Ro but cannot larger than cell[i])
    rc2 = Rcut*Rcut;   //Rcut^2
    //calculation volume
    V=cell[0]*cell[1]*cell[2]; //in A^3
    //Tail correction for Ep
    double Ro_Rc3=pow(Ro/Rcut,3.0);
    double Ro_Rc6=Ro_Rc3*Ro_Rc3;
    double Ro_Rc9=Ro_Rc6*Ro_Rc3;
    double Ro3=Ro*Ro*Ro;
    Eptail=(4.0/3.0)*3.1415926*npart*npart/V*Do*Ro3*(Ro_Rc9/6.0-Ro_Rc3);
    //Tail correction for pressure
    double unitc=4184*1E21/6.0221367E23; //conver from kcal/mol/A3 to GPa
    Ptail=(8.0/3.0)*3.1415926*npart*npart/V/V*Do*Ro3*(Ro_Rc9/3.0-Ro_Rc3)*unitc;
    
    if(gen_gro) strcpy(trjfname,"256Ar.gro");  //filename of the MD trajectory

    //allocation of memerory
    x=new double [npart][3];        //position in Angstrom
    v=new double [npart][3];         //velocity in Angstrom/ps
    f=new double [npart][3];          //force in kcal/mol/A
    fm=new double [npart][3];        //force in next step
    m=new double [npart];             //mass in g/mo


    //assign mass of particles
    int i,j,k;
    for(i=0;i<npart;i++) m[i]=mass;//in g/mo

    //One simple way to place particles in space
    double sep=3.5;  //separation
    int nt=0;
    if(period){
               double delta[3];
               int pps;
               pps = int(pow(npart,1.0/3.0)) + 1; // particles per side
               for(k=0;k<3;k++)delta[k] = cell[k]/pps; //spacing of particles
               for(i=0;i<pps;i++) {
                 for(j=0;j<pps;j++) {
                   for(k=0;k<pps;k++) {
                     x[nt][0]=i*delta[0];
                     x[nt][1]=j*delta[1];
                     x[nt][2]=k*delta[2];
                     nt++; //number of particles placed
                     if(nt==npart) i=j=k=pps; //nt has reached npart, exit loops
                     }
                     }
                     }
               }
    else{
      for(i=0;i<npart;i++) {
        for(j=0;j<=i;j++) {
            for(k=0;k<=j;k++) {
                x[nt][0]=i*sep;
                x[nt][1]=j*sep;
                x[nt][2]=k*sep;
                nt++;
                if(nt==npart) i=j=k=npart;
                }
            }
        }
     }

    //Assign velocities
    double cmv[3],sumv2,tmass,Ti,fs;
    srand(123);
    cmv[0]=cmv[1]=cmv[2]=sumv2=tmass=0;
    df= 3*npart; //degree of freedom 3N

    for(i=0;i<npart;i++) {
        tmass += m[i]; //total mass
        for(k=0;k<3;k++) {
            v[i][k]= myrand01()-0.5; //random number between -0.5 and 0.5
            cmv[k]+= m[i]*v[i][k];   //center of mass velocity
            sumv2 += m[i]*v[i][k]*v[i][k];
        }
    }

    for(k=0;k<3;k++) { cmv[k]/=tmass; sumv2 -= tmass*cmv[k]*cmv[k]; }
    Ti = sumv2/(df*8.314*0.1);
    fs = sqrt(Tset/Ti);  //scale factor

    for(i=0;i<npart;i++) {
        for(k=0;k<3;k++) {
            v[i][k]  = (v[i][k]-cmv[k]) *fs;  //initial velocity
        }
    }   

    force();   //initial force calculation

    return 0;
}

int MODEL::force()
{  //The function determines the net force on each particle

    int i,j,k;
    for(i=0;i<npart;i++) {
        for(k=0;k<3;k++) f[i][k]=0;   //set forces to zero
    }
    for(i=0;i<6;i++) strs[i]=0; //set stress to zero

    // consider atom j at origin, the force on atom i at some position r
    // fx = - dU/dx = -(dU/dr)(dr/dx)= - (x/r)(dU/dr)
    // U = Do ( (Ro/r)^12 - 2 (Ro/r)^6) )
    // dU/dr = -12 Do/r ( (Ro/r)^12 - (Ro/r)^6) )
    // fx = 12 x Do/r^2 ( (Ro/r)^12 - (Ro/r)^6) )
    //    =  x ( 12DoRo^12/r^6 - 12DoRo^6 )/r^6 /r^2

    double r2,r2i,r6i,ff,xr[3];
    Ep=0;
    double nbox;
    for(i=0;i<npart;i++) {
        for(j=i+1;j<npart;j++) {
            for(k=0;k<3;k++) xr[k]= x[i][k] - x[j][k];  //distance vector
              if(period) { //periodic system, find dist within one cell
                for(k=0;k<3;k++) {
                  xr[k] += 0.5*cell[k]; //shift j to box center
                  if(xr[k]<0) nbox = int( (xr[k]/cell[k]) -1); //determine image
                    else nbox=int(xr[k]/cell[k]);
                  xr[k] -= ( nbox+0.5)*cell[k]; //map to original box and shift
                }
            }
            r2= xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];   //distance squared
            if(r2<rc2 || period==0) { //within cutoff distance
              r2i= 1/r2;
              r6i= r2i*r2i*r2i;
              ff = r6i*(p0*r6i-p1)*r2i;
              for(k=0;k<3;k++) {
                    f[i][k]+= ff*xr[k];
                    f[j][k]-= ff*xr[k];  //Newton's 3rd law
              }
              //the stress tensor
              strs[0]+= ff*xr[0]*xr[0]; //xx in unit of kcal/mol
              strs[1]+= ff*xr[1]*xr[1]; //yy
              strs[2]+= ff*xr[2]*xr[2]; //zz
              strs[3]+= ff*xr[0]*xr[1]; //xy
              strs[4]+= ff*xr[0]*xr[2]; //xz
              strs[5]+= ff*xr[1]*xr[2]; //yz
              Ep += (p0*r6i - p1*2)*r6i/12.0;
            }
          }
    }
    return 0;
}

int MODEL::integrate()
{    //velocity Verlet integration for particle position and velocity

    int i,k;
    //double xx,*pv;
    for(i=0;i<npart;i++) {
        for(k=0;k<3;k++) {
           // force in kcal/mol/A, mass in g/mol, x in A, dt in ps, v in A/ps
           x[i][k] +=  v[i][k]*dt + f[i][k]*dt*dt*418.4/(2*m[i]);
           fm[i][k]= f[i][k];
        }

    }
    force();
    for(i=0;i<npart;i++) {
        for(k=0;k<3;k++) {
           v[i][k] += (f[i][k]+fm[i][k])*dt*418.4/(2*m[i]);
        }
    }
    istep++;    //current step
    t=istep*dt; //current time in ps


    return 0;
}

void MODEL::output()
{

//    if((istep%nsamp)) return 0;

    //calculation of system temperature and kinetic energy
    int i;
    T=Ek=0;
    for(i=0;i<npart;i++) {
       Ek += (m[i]*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]));
    }

    T = Ek/(df*8.314*0.1);  //0.1 is due to unit conversion
    Ek /= (2*418.4);        //in kcal/mol
    Ep += Eptail;
    Etot=Ek+Ep;
    //calcualate pressure
    double unitc=4184*1E21/6.0221367E23; //conver from kcal/mol/A3 to GPa
    P=((df*1.380658E-2*T)+(strs[0]+strs[1]+strs[2])*unitc)/(3.0*V); //in GPa
    P+= Ptail;
    char null[1024];
    sprintf(null,"%d     %.3f     %.2f     %.2f     %.2f\n",npart,t,T,P*1000,Etot);    
    cout << null;
}

int MODEL::sample()
{
    if((istep%nsamp)) return 0;

    //calculation of system temperature and kinetic energy
    int i;
    T=Ek=0;
    for(i=0;i<npart;i++) {
       Ek += (m[i]*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]));
    }

    T = Ek/(df*8.314*0.1);  //0.1 is due to unit conversion
    Ek /= (2*418.4);        //in kcal/mol
    Ep += Eptail;
    Etot=Ek+Ep;
    //calcualate pressure
    double unitc=4184*1E21/6.0221367E23; //conver from kcal/mol/A3 to GPa
    P=((df*1.380658E-2*T)+(strs[0]+strs[1]+strs[2])*unitc)/(3.0*V); //in GPa
    P+= Ptail;

    //Display current information
    char null[1024];
    sprintf(null,"Step %d, Remain %d, Time: %.3f ps, ",istep,nstep-istep,t);
//    cout<<null<<endl;
    logf<<null;
    sprintf(null,"T %.2f K, P %.2f MPa, V %.0f A3, ",T,P*1000,V);
//    cout<<null<<endl;;
    logf<<null;
    sprintf(null,"E(kcal/mol) Ek %.2f Ep %.2f Et %.2f",Ek,Ep,Etot);
 //   cout<<null<<endl;
    logf<<null<<endl;
//    sprintf(null,"Tail contribution %.0f%% in Ep %.0f%% in P ",Eptail/Ep*100,Ptail/P*100);
//    cout<<null<<endl;

    if(gen_gro) dump_trj();
    return 0;
}

int MODEL::dump_trj()
{
    char null[1024];
    int i;

    sprintf(null,"My MD trj: Current Step %d, Time: %f ps",istep,t);
    outf<<null<<endl;
    sprintf(null,"%5d",npart);
    outf<<null<<endl;
    for(i=0;i<npart;i++) {
        sprintf(null,"%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f",1,"O","O",i+1,x[i][0]/10,x[i][1]/10,x[i][2]/10,v[i][0]/10,v[i][1]/10,v[i][2]/10);
        outf<<null<<endl;
    }
    sprintf(null,"%10.5f%10.5f%10.5f",cell[0]/10,cell[1]/10,cell[2]/10);
    outf<<null<<endl;
    return 0;
}

double MODEL::myrand01()
{
    return rand()/double(RAND_MAX); //returns a number between 0 (inclusive) and 1 (inclusive)
}

MODEL::MODEL()
{

};

 MODEL::~MODEL()
{
    delete [] x;
    delete [] v;
    delete [] f;
    delete [] m;
    delete [] fm;

    outf.close();
    logf.close();
};

