#include <fftw3.h>
#include <string.h>
#include <random>
#include <complex.h>

#define REAL double
#define COMPLEX double _Complex

REAL a;
int seed;
int flag;

int NX, NY;
REAL DX, DY;
REAL LX, LY;
int NR;
int NK; 

REAL bp;      // magnitude of the partial Burger's vector
REAL diss;    // dissociation distance
REAL phi;     // character angle in degrees phi=0 (screw), phi=90 (edge)
REAL theta_l; // character of the leading partial
REAL theta_t; // character of the trailing partial

int E_per;    // flag to indicate if the generated force field satisfy periodicity of energy in the y direction


// normalized auto-correlation function in Fourier space for sigma_yz
REAL SigmaE_k(REAL kx, REAL ky)
{
    REAL k2 = kx*kx + ky*ky;
    return 60.*a*a*a*sqrt(M_PI)/(LX*LY)*ky*ky*(-a*sqrt(M_PI)*exp(-a*a*k2) + M_PI*(1.+2.*a*a*k2)*erfc(a*sqrt(k2))/(2.*sqrt(k2)));
}

// normalized auto-correlation function in Fourier space for sigma_xz
REAL SigmaS_k(REAL kx, REAL ky)
{
    REAL k2 = kx*kx + ky*ky;
    return 60.*a*a*a*sqrt(M_PI)/(LX*LY)*kx*kx*(-a*sqrt(M_PI)*exp(-a*a*k2) + M_PI*(1.+2.*a*a*k2)*erfc(a*sqrt(k2))/(2.*sqrt(k2)));
}

// normalized cross-correlation function between sigma_xz and sigma_yz in Fourier space
REAL SigmaCross_k(REAL kx, REAL ky)
{
    REAL k2 = kx*kx + ky*ky;
    return 60.*a*a*a*sqrt(M_PI)/(LX*LY)*kx*ky*(-a*sqrt(M_PI)*exp(-a*a*k2) + M_PI*(1.+2.*a*a*k2)*erfc(a*sqrt(k2))/(2.*sqrt(k2)));
}

// total force acting on the dislocation
REAL Force_k_diss(REAL kx, REAL ky)
{
    REAL bel=bp*sin(theta_l);
    REAL bsl=bp*cos(theta_l);
    REAL bet=bp*sin(theta_t);
    REAL bst=bp*cos(theta_t);
    REAL res = 0;
    res += bst*bst*SigmaS_k(kx,ky); 
    res += bet*bet*SigmaE_k(kx,ky);
    res += bsl*bsl*SigmaS_k(kx,ky); 
    res += bel*bel*SigmaE_k(kx,ky); 
    res += 2*bst*bsl*SigmaS_k(kx,ky)*cos(ky*diss);      
    res += 2*bet*bel*SigmaE_k(kx,ky)*cos(ky*diss);     
    res += 2*bst*bet*SigmaCross_k(kx,ky);
    res += 2*bsl*bel*SigmaCross_k(kx,ky);
    res += 2*bst*bel*SigmaCross_k(kx,ky)*cos(ky*diss);  
    res += 2*bsl*bet*SigmaCross_k(kx,ky)*cos(ky*diss);
    
//    if((kx==0) and (ky==0))
//    {
//        printf("%f %f %f %f %f %f\n", bel, bsl, bet, bst, creal(res), cimag(res));
//    }
    
    return res;
}


// total force acting on the dislocation for d=0
REAL Force_k_undiss(REAL kx, REAL ky)
{
    REAL be=sin(phi);
    REAL bs=cos(phi);

    REAL res = 0;
    res += bs*bs*SigmaS_k(kx,ky); 
    res += be*be*SigmaE_k(kx,ky);
    res += 2*be*bs*SigmaCross_k(kx,ky);
    
//    if((kx==0) and (ky==0))
//    {
//        printf("%f %f %f %f %f %f\n", bel, bsl, bet, bst, creal(res), cimag(res));
//    }
    
    return res;
}


float ReverseFloat( const float inFloat )
{
   float retVal;
   char *FloatToConvert = ( char* ) & inFloat;
   char *returnFloat = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = FloatToConvert[3];
   returnFloat[1] = FloatToConvert[2];
   returnFloat[2] = FloatToConvert[1];
   returnFloat[3] = FloatToConvert[0];
   
   return retVal;
}


void Write_field_vtk(REAL *F, const char *fname)
{        
	FILE *OutFile=fopen(fname,"w");
	
	fprintf(OutFile,"# vtk DataFile Version 2.0\n");
	fprintf(OutFile,"iter %d\n",0);
	fprintf(OutFile,"BINARY\n");
	fprintf(OutFile,"DATASET STRUCTURED_POINTS\n");
	fprintf(OutFile,"ORIGIN \t %f %f %f\n",0.,0.,0.);
	fprintf(OutFile,"SPACING \t %f %f %f\n", DX, DY, 1.0);
	fprintf(OutFile,"DIMENSIONS \t %d %d %d\n", NX, NY, 1);
	fprintf(OutFile,"POINT_DATA \t %d\n",NR);
	fprintf(OutFile,"SCALARS \t volume_scalars float 1\n");
	fprintf(OutFile,"LOOKUP_TABLE \t default\n");
	
	float temp;
	int ind;
	
    for(int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            ind = NY*i + j;
            temp=ReverseFloat(float(F[ind]));
            fwrite(&temp, sizeof(float), 1, OutFile);
        }
    }
	
	fclose(OutFile);
}

void Write_field_txt(REAL *F, const char *fname)
{        
	FILE *OutFile=fopen(fname,"w");
	
	fprintf(OutFile,"#SPACING \t %f %f %f\n", DX, DY, 1.0);
	fprintf(OutFile,"#DIMENSIONS \t %d %d %d\n", NX, NY, 1);
	fprintf(OutFile,"#POINT_DATA \t %d\n",NR);
	fprintf(OutFile,"#ORDERING \t in column (Ny*i+j)\n");
	
	float temp;
	int ind;
	
    for(int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            ind = NY*i + j;
            fprintf(OutFile, "%f\n", float(F[ind]));
        }
    }
	
	fclose(OutFile);
}


int main(int argc, char **argv)
{
    int ind = 0;
    REAL kx,ky;

    // read input command line

    NX = atoi(argv[1]);     // dimension along x
    NY = atoi(argv[2]);     // dimension along y

    DX = atof(argv[3]);     // grid spacing [AA]
    DY = atof(argv[4]);     // grid spacing [AA]
    
    a = atof(argv[5]);      // spreading parameter for the core structure [AA]
    diss = atof(argv[6]);  // dissociation distance [AA]
    phi = M_PI/180.*atof(argv[7]);    // dislocation charater converted to [rad]
    
    E_per = atof(argv[8]); // flag for periodic energy profiles
    
    seed = atof(argv[9]);   // random seed
    char* fname_out = argv[10];  // name of the output vtk 

    LX = NX*DX;
    LY = NY*DY;
    NR = NX*NY;
    NK = NX*(NY/2+1);
    
    theta_l = phi+M_PI/6.;
    theta_t = phi-M_PI/6.;

    bp=1.0/sqrt(3.);   // we consider btot=1.0
        
    REAL *Rr;          // pointer of the real space data
    COMPLEX *Rk;       // pointer of the fourier space data

    fftw_plan plan_R_r2c, plan_R_c2r;  // fft plan

    // allocate
    Rr = (REAL*) fftw_malloc(sizeof(REAL)*NR);
    Rk = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*NK);
        
    // prepare plans    
    plan_R_r2c = fftw_plan_dft_r2c_2d(NX, NY, Rr, reinterpret_cast<fftw_complex*>(Rk), FFTW_ESTIMATE);
    plan_R_c2r = fftw_plan_dft_c2r_2d(NX, NY, reinterpret_cast<fftw_complex*>(Rk), Rr, FFTW_ESTIMATE);

    // define the noise in Fourier space based on the auto-correlation function [Geslin et al. JMPS 2021]
    REAL N_k, M_k;
        
    std::default_random_engine generator(seed);
    std::normal_distribution<REAL> distribution(0.0,1.0);    
        
        
    ///////////////////////////////////
    // generate correlated component //
    ///////////////////////////////////
    for(int i=0; i<NX; i++)
    {
        for(int j=0; j<(NY/2+1); j++)
        {
             ind = (NY/2+1)*i + j;
                    
             kx = 2.*M_PI/LX*REAL(i);
             if(i>NX/2-1)
             {
                 kx = 2.*M_PI/LX*REAL(i-NX);
             }
             ky = 2*M_PI/LY*REAL(j);
                    
             // random numbers
             N_k = distribution(generator);
             M_k = distribution(generator);
            
             if(kx*ky>=0)
             {
                 N_k = N_k;
                 M_k = M_k;
             }
             else
             {
                 N_k = -N_k;
                 M_k = -M_k;
             }

             if(diss>0) //positive dissociation
             {
                 // /!\ special case for j=0 and j==NY/2 because of folding of C2R Fourier transform
                 if(j==0) 
                 {
                     if(E_per=1)
                     {
                        Rk[ind] = 0;
                     }
                     else
                     {
                        Rk[ind] = sqrt(Force_k_diss(kx,ky))*(N_k + M_k*_Complex_I);
                     }
                 }
                 else if(j==NY/2)
                 {
                     Rk[ind] = sqrt(Force_k_diss(kx,ky))*(N_k + M_k*_Complex_I);
                 }
                 else
                 {
                     Rk[ind] = sqrt(Force_k_diss(kx,ky)/2.)*(N_k + M_k*_Complex_I);
                 }
             }
             else if(diss==0.0) // null dissociation
             {
                 // /!\ special case for j=0 and j==NY/2 because of folding of C2R Fourier transform
                 if(j==0) 
                 {
                     if(E_per=1)
                     {
                        Rk[ind] = 0;
                     }
                     else
                     {
                        Rk[ind] = sqrt(Force_k_undiss(kx,ky))*(N_k + M_k*_Complex_I);
                     }
                 }
                 else if(j==NY/2)
                 {
                     Rk[ind] = sqrt(Force_k_undiss(kx,ky))*(N_k + M_k*_Complex_I);
                 }
                 else
                 {
                     Rk[ind] = sqrt(Force_k_undiss(kx,ky)/2.)*(N_k + M_k*_Complex_I);
                 }
             }
         }
     }
     Rk[0] = 0;
        
        
    // FFT back to real space    
    fftw_execute(plan_R_c2r);
            
    // print out variance of the fields
//    REAL var = 0;
//    for(int ind=0; ind<NR; ind++)
//    {
//        var += Rr[ind]*Rr[ind];
//    }
//    var = var/REAL(NR); 
//    printf("%f\n",float(var));
       
    
    Write_field_vtk(Rr, fname_out);
//    Write_field_txt(Rr, fname_out);

    return 0;
}
