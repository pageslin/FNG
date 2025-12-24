#include <fftw3.h>
#include <string.h>
#include <random>
#include <complex.h>

#define REAL double
#define COMPLEX double _Complex

int seed;
int flag;

int NX, NY;
REAL DX, DY;
REAL LX, LY;
int NR;
int NK; 

REAL sigma;
int E_per;


// normalized auto-correlation function in Fourier space
REAL S_uncorr_k(REAL kx, REAL ky)
{
    return (DX*DY)/(LX*LY);
}    


// normalized auto-correlation function in Fourier space with Gaussian convolution
REAL S_uncorr_conv_k(REAL kx, REAL ky)
{
    REAL k2 = kx*kx + ky*ky;
    return (2*M_PI)*sigma*sigma/(LX*LY)*exp(-k2*sigma*sigma/2.);
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

double ReverseDouble( const double inDouble )
{
   double retVal;
   char *DoubleToConvert = ( char* ) & inDouble;
   char *returnDouble = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnDouble[0] = DoubleToConvert[7];
   returnDouble[1] = DoubleToConvert[6];
   returnDouble[2] = DoubleToConvert[5];
   returnDouble[3] = DoubleToConvert[4];
   returnDouble[4] = DoubleToConvert[3];
   returnDouble[5] = DoubleToConvert[2];
   returnDouble[6] = DoubleToConvert[1];
   returnDouble[7] = DoubleToConvert[0];
      
   return retVal;
}


void Write_field(REAL *F, const char *fname)
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


int main(int argc, char **argv)
{
    int ind = 0;
    REAL kx,ky;

    // read input command line

    NX = atoi(argv[1]);     // dimension along x 
    NY = atoi(argv[2]);     // dimension along y

    DX = atof(argv[3]);     // grid spacing [AA]
    DY = atof(argv[4]);     // grid spacing [AA]
    
    sigma = atof(argv[5]);  // spreading length for Gaussian spreading [AA]
    
    E_per = atoi(argv[6]);  // E_per=1 if periodic energy landscape in y, E_per=0 otherwise 
    
    seed = atof(argv[7]);  // random seed
    char* fname_out = argv[8];  // name of the output vtk 

    
    LX = NX*DX;
    LY = NY*DY;

    NR = NX*NY;
    NK = NX*(NY/2+1);
        
    REAL *Rr;          // pointer of the real space data
    COMPLEX *Rk;       // pointer of the fourier space data

    fftw_plan plan_R_c2r;  // fft plan

    // allocate
    Rr = (REAL*) fftw_malloc(sizeof(REAL)*NR);
    Rk = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*NK);
        
        
    // prepare FFT plan   
    plan_R_c2r = fftw_plan_dft_c2r_2d(NX, NY, reinterpret_cast<fftw_complex*>(Rk), Rr, FFTW_ESTIMATE);


    // define the noise in Fourier space based on the auto-correlation function
    REAL Nk, Mk;
        
    std::default_random_engine generator(seed);
    std::normal_distribution<REAL> distribution(0.0,1.0);    
        
    ///////////////////////////////////////////////////
    // generate random noise without any convolution //
    ///////////////////////////////////////////////////
    if(sigma<=0)
    {
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
                Nk = distribution(generator);
                Mk = distribution(generator);
                        
                // /!\ special case for j=0 and j==NY/2 because of folding of C2R Fourier transform
                if(j==0) 
                {
                    if(E_per==1)
                    {
                        Rk[0]=0.0 ; // compatible with periodic energy landscape

                    }
                    else
                    {
                        Rk[ind] = sqrt(S_uncorr_k(kx,ky))*(Nk+Mk*_Complex_I);                        
                    }
                }
                else if(j==NY/2)
                {
                    Rk[ind] = sqrt(S_uncorr_k(kx,ky))*(Nk+Mk*_Complex_I);
                }
                else
                {
                    Rk[ind] = sqrt(S_uncorr_k(kx,ky)/2.)*(Nk+Mk*_Complex_I);
                }
            }
        }
    }
    else
    {
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
                Nk = distribution(generator);
                Mk = distribution(generator);
                        
                // /!\ special case for j=0 and j==NY/2 because of folding of C2R Fourier transform
                if(j==0) 
                {
                    if(E_per==1)
                    {
                        Rk[0]=0.0 ; // compatible with periodic energy landscape
                    }
                    else
                    {
                        Rk[ind] = sqrt(S_uncorr_conv_k(kx,ky))*(Nk+Mk*_Complex_I);
                    }
                }
                else if(j==NY/2)
                {
                    Rk[ind] = sqrt(S_uncorr_conv_k(kx,ky))*(Nk+Mk*_Complex_I);
                }
                else
                {
                    Rk[ind] = sqrt(S_uncorr_conv_k(kx,ky)/2.)*(Nk+Mk*_Complex_I);
                }
            }
        }
    }
    
    
        
    // FFT back to real space    
    fftw_execute(plan_R_c2r);
          
    // print out variance of the fields
    REAL var = 0;
    for(int ind=0; ind<NR; ind++)
    {
        var += Rr[ind]*Rr[ind];
    }
    var = var/REAL(NR); 

    // ouput vtk files
    Write_field(Rr, fname_out);
        
    printf("%f\n",var);

    return 0;

}
