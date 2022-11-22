#include <math.h>
#include <stdlib.h>

double runif(); 
double rnorm(double mu, double sigma);

double trap2(double x, double a, double b, double c, double d); 
double trap1(double x, double a, double b); 

double MyMin(double x, double y);
double MyMax(double x, double y);
double MyMin3(double x, double y, double z);

double F_dMs_dt(double Gs, double Kl, double KM, double Ms);
double F_Gs(double gs, double Ms, double Cs, double Ns);
double F_dMr_dt(double Gr, double Kl, double KM, double Mr);
double F_Gr(double gs, double Mr, double Cr, double Nr);
double F_P(double A0, double Ms, double KA, double Cs, double Jc);
double F_Un(double N0, double Mr, double KA, double Nr, double Jn);
double F_TAUc(double Cs, double Cr, double Ms, double Mr, double RsC, double RrC);
double F_TAUn(double Ns, double Nr, double Ms, double Mr, double RsN, double RrN);
double F_RsC(double RHOc, double Ms, double q);
double F_RrC(double RHOc, double Mr, double q);
double F_RsN(double RHOn, double Ms, double q);
double F_RrN(double RHOn, double Mr, double q);
double F_dCs_dt(double P, double Fc, double Gs, double TAUc);
double F_dCr_dt(double Fc, double Gr, double TAUc);
double F_dNs_dt(double Fn, double Gs, double TAUn);
double F_dNr_dt(double Un, double Fn, double Gr, double TAUn);


void ThornTimeARU( int *steps, double *initials, 
            double const *TAIR, double const *TSOIL, double const *M, double const *N, 
            double const *FIRE, double const *A,
            double *Kl, double *gs, double *gr, double *KM, double *A0, double *N0,
            double *KA, double *Jc, double *Jn, double *q, double *RHOc, double *RHOn,
            double *Fc, double *Fn, 
            //  to be fitted environmental dependencies
            double *ma1,double *ma2,
            double *tn1,double *tn2, 
            double *mn1,double *mn2, double *mn3, double *mn4,
            double *tg1,double *tg2, double *tg3,double *tg4,
            double *mg1,double *mg2,
            double *tr1,double *tr2,
            double *f1, double *f2,            
            // uncertainty parameters
            double *uMs, double *uMr, double *uCs, double *uCr,
            double *uNs, double *uNr,
            // return
            double *ABUTIME
);

double runif() {
  // return a uniformly distributed random value
  return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}

double rnorm(double mu, double sigma) {
  // return a normally distributed random value
  double v1=runif();
  double v2=runif();
  double nstd=cos(2*3.14*v2)*sqrt(-2.*log(v1));
  return (mu + sigma * (double) nstd );
}

double MyMin(double x, double y)
{
  if (x > y) return y; 
  else return x;
}

double MyMin3(double x, double y, double z)
{
  int i;
  double xyz[3] ={x,y,z};
  double min;
  min = xyz[0];
  for(i = 1; i < 3; i++)
  {
    if(xyz[i] < min)   min = xyz[i];
  }
  return min;
}

double MyMax(double x, double y)
{
  if (x > y) return x; 
  else return y;
}

double trap2(double x, double a, double b, double c, double d) 
{
  double ret;
  ret = MyMax( MyMin3( (x-a)/(b-a), 1, (d-x)/(d-c) ), 0 );
  return ret;
}

double trap1(double x, double a, double b) 
{
  double ret;
  ret = MyMax( MyMin( (x-a)/(b-a), 1 ), 0 );
  return ret;
}

//-----Thornley Model functions	
double F_dMs_dt(double Gs, double Kl, double KM, double Ms)
{
  double ret;
  ret = Gs - ( Kl*Ms ) / ( 1.0 + KM / Ms );
  return ret;
}

double F_Gs(double gs, double Ms, double Cs, double Ns)
{
  double ret;
  ret = gs * (Cs * Ns) / Ms;
  return ret;
}

double F_dMr_dt(double Gr, double Kl, double KM, double Mr)
{
  double ret;
  ret = Gr - ( Kl*Mr ) / ( 1.0 + KM / Mr );
  return ret;
}

double F_Gr(double gs, double Mr, double Cr, double Nr)
{
  double ret;
  ret = gs * (Cr * Nr) / Mr;
  return ret;
}

double F_P(double A0, double Ms, double KA, double Cs, double Jc)
{
  double ret;
  ret = ( A0*Ms ) / ( (1.0 + Ms / KA) * (1.0 + Cs / (Jc*Ms))  );
  return ret;
}

double F_Un(double N0, double Mr, double KA, double Nr, double Jn)
{
  double ret;
  ret = ( N0*Mr ) / ( (1 + Mr / KA) * (1 + Nr / (Jn*Mr))  );
  return ret;
}

double F_TAUc(double Cs, double Cr, double Ms, double Mr, double RsC, double RrC)
{
  double ret;
  ret = ( Cs/Ms - Cr/Mr ) / ( RsC + RrC);
  return ret;
}

double F_TAUn(double Ns, double Nr, double Ms, double Mr, double RsN, double RrN)
{
  double ret;
  ret = ( Nr/Mr - Ns/Ms ) / ( RsN + RrN);
  return ret;
}

double F_RsC(double RHOc, double Ms, double q)
{
  double ret;
  ret = RHOc / pow(Ms,q);
  return ret;
}

double F_RrC(double RHOc, double Mr, double q)
{
  double ret;
  ret = RHOc / pow(Mr,q);
  return ret;
}

double F_RsN(double RHOn, double Ms, double q)
{
  double ret;
  ret = RHOn / pow(Ms,q);
  return ret;
}

double F_RrN(double RHOn, double Mr, double q)
{
  double ret;
  ret = RHOn / pow(Mr,q);
  return ret;
}

double F_dCs_dt(double P, double Fc, double Gs, double TAUc)
{
  double ret;
  ret  = P - Fc * Gs - TAUc;
  return ret;
}

double F_dCr_dt(double Fc, double Gr, double TAUc)
{
  double ret;
  ret  = TAUc - Fc * Gr;
  return ret;
}

double F_dNs_dt(double Fn, double Gs, double TAUn)
{
  double ret;
  ret  = TAUn - Fn * Gs;
  return ret;
}

double F_dNr_dt(double Un, double Fn, double Gr, double TAUn)
{
  double ret;
  ret  = Un - Fn * Gr - TAUn;
  return ret;
}

void ThornTimeARU( int *steps, double *initials, 
            double const *TAIR, double const *TSOIL, double const *M, double const *N, 
            double const *FIRE, double const *A,
            double *Kl, double *gs, double *gr, double *KM, double *A0, double *N0,
            double *KA, double *Jc, double *Jn, double *q, double *RHOc, double *RHOn,
            double *Fc, double *Fn, 
            //  to be fitted environmental dependencies
            double *ma1,double *ma2,
            double *tn1,double *tn2, 
            double *mn1,double *mn2, double *mn3, double *mn4,
            double *tg1,double *tg2, double *tg3,double *tg4,
            double *mg1,double *mg2,
            double *tr1,double *tr2,
            double *f1, double *f2,
            // uncertainty parameters
            double *uMs, double *uMr, double *uCs, double *uCr,
            double *uNs, double *uNr,
            // return
            double *ABUTIME
                     )
{
	int		i;
	int		index;
	double	Ms_ ;
	double	Mr_ ;
	double	Cs_ ;
	double	Cr_ ;
	double	Ns_ ;
	double	Nr_ ;
	double  A0_E;
	double  N0_E;
	double  gs_E;
	double  gr_E;
	double  RsC;
	double  RrC;
	double  RsN;
	double  RrN;
	double  P;
	double  Un;
	double  Gs;
	double  Gr;
	double  TAUc;
	double  TAUn;
	double  sumb;
	double  LOSS;
    double  GF;
	
    sumb=0;

    Ms_ = *initials; 
    Mr_ = *initials; 
    Cs_ = Ms_ * 0.05;
    Cr_ = Mr_ * 0.05;
    Ns_ = Ms_ * 0.01;
    Nr_ = Ms_ * 0.01;

    for(i = 0; i < *steps; i ++) 
    { 
            index = i;  
            
            A0_E = *A0* A[index]*  trap1(M[index],*ma1,*ma2) ;
            
            N0_E = *N0* trap1(TSOIL[index],*tn1,*tn2) *
                        trap2(M[index],*mn1,*mn2,*mn3,*mn4) ; 
            
            GF = trap2(TSOIL[index],*tg1,*tg2,*tg3,*tg4)*trap1(M[index],*mg1,*mg2);             

            gs_E = *gs* GF;
            gr_E = *gr* GF; 
            
            RsC = F_RsC(*RHOc,Ms_,*q);
            RrC = F_RrC(*RHOc,Mr_,*q);
            RsN = F_RsN(*RHOn,Ms_,*q);
            RrN = F_RrN(*RHOn,Mr_,*q);
            
            P  = F_P(A0_E,Ms_,*KA,Cs_,*Jc);
            Un = F_Un(N0_E,Mr_,*KA,Nr_,*Jn);
            Gs = F_Gs(gs_E,Ms_,Cs_,Ns_);
            Gr = F_Gr(gr_E,Mr_,Cr_,Nr_);
            
            TAUc = F_TAUc(Cs_,Cr_,Ms_,Mr_,RsC,RrC);
            TAUn = F_TAUn(Ns_,Nr_,Ms_,Mr_,RsN,RrN);
            
            LOSS = *Kl*0.5 + *Kl*0.5*trap1(TAIR[index],*tr1,*tr2);
            
            Ms_ = MyMax(0.0, Ms_ + F_dMs_dt(Gs,LOSS,*KM,Ms_) + rnorm(0.0, *uMs) );
            Mr_ = MyMax(0.0, Mr_ + F_dMr_dt(Gr,LOSS,*KM,Mr_) + rnorm(0.0, *uMr) );
            Cs_ = MyMax(0.0, Cs_ + F_dCs_dt(P,*Fc,Gs,TAUc)   + rnorm(0.0, *uCs) );
            Cr_ = MyMax(0.0, Cr_ + F_dCr_dt(*Fc,Gr,TAUc)     + rnorm(0.0, *uCr) );
            Ns_ = MyMax(0.0, Ns_ + F_dNs_dt(*Fn,Gs,TAUn)     + rnorm(0.0, *uNs) );
            Nr_ = MyMax(0.0, Nr_ + F_dNr_dt(Un,*Fn,Gr,TAUn)  + rnorm(0.0, *uNr) );
            
            Ms_ = ( 1.0 - trap1(FIRE[index],*f1,*f2) ) * Ms_; 
            
            sumb = MyMax(0.0, Ms_ + Mr_);
            ABUTIME[index]=sumb;

    }
}



