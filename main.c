#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//Basic constants (SI)
#define PI (M_PI) //
#define h (6.626e-34) //J*s Plank
#define h_ (h/2/PI) //J*s Plank_
#define m0 (9.109e-31) //Free electron mass, kg
#define e (-1.602e-19) //Unit charge, Cl
#define kb (1.38e-23) //Boltzmann constant, J*K
#define eps0 (8.854187817e-12) //Vacuum permittivity, F/m
#define VxcA (11.4) //Capasso p46 parameter A for excahge-correlation potential
#define VxcB (0.6213) //Capasso p46 parameter B for excahge-correlation potential
#define aB (5.29e-11) //Bohr radius, m
//Parameters and intermediates

double T; //Temperature, K
double kT;
double epsB, epsW; //dielectric permittivity in barrier and well
double mB,mW; //Effective electron masses in barrier and in well
double mtB,mtW,mlB,mlW; //Transverse and longitudinal valley masses in Barrier and well
double a;//well width
double aBeW; //effective Bohr radius in well
double ns;//free carrier concentration, m-2
double Gamma; //Broadening for Delta
//Size Parameters
#define WFpointsMAX 100000
#define EigenCountMAX 30

// Size indicators and arrays

int EigenCount; //Number of s-q states in QW
double EigenValues [EigenCountMAX+1]; //s-q states in QW (E-Ef) - num from 1
int WFpointsCount;
double WFTable[WFpointsMAX][EigenCountMAX+1]; //[][0] - z in A;  wave functions in A^-1/2.
double ASTable[EigenCountMAX+1][EigenCountMAX+1]; // Tables of Collective Absorption Shift in eV
double MatrixElements [EigenCountMAX+1][EigenCountMAX+1];
//column 0 - z in Angstroms columns 1-* - WF (Normalized?)

void InitializeParameters() // can be redefined
{
    epsB=12.6;
    epsW=11.7;
    mtW=0.19; //m parallel (JAP) Si stressed
    mlW=0.92; //m perpendicular (JAP) Si stressed
    mtB=0.19; //SiGe not stressed. Proportion between Si0.8, Ge0.2
    mlB=0.92; //Si Ge not stressed.
    // Si: mt=0.19, ml=0.98; Ge: mt=0.0815, ml=1.59;
    mW=2*pow(mtW*mtW*mlW,1.0/3.0); // Density of states mass m* - why not like a transport?
    mB=2*pow(mtB*mtB*mlB,1.0/3.0);// Non-deformed, 6 subbands
    T=4;
    kT=kb*T;
    ns=0;
    a=200*1e-10; //Well width 200A in meters
    Gamma=0.005;// 5 meV
    aBeW=epsW/mW*aB;// effective Bohr radius
    EigenCount=0;
}

double Distribution(double Energy)
{
    return 1/(exp(fabs(Energy*e)/kT)+1); //EF=0 all over the model. Turn from eV.
}

void ReadStatus(char * StatusFile) //material parameters should be in accordance
{
    char s [101], s1[100], s2[100], s3[100];
    int i;
    double d1, V, ns_add;
    FILE *stream;
    stream = fopen(StatusFile, "r");
    if (stream==NULL)
    {
        printf("Error! Status file %s not found.\n",StatusFile);
        return;
    }
    ns=0; a=0; //free el concentr and well width
    while (fgets(s,100,stream)!=NULL)
    {
        // 95.0 Ang	SiStressed	ns= 8.111e+010 cm-2,	ps= 9.500e-017 cm-2
        if (sscanf(s,"%lf %s %s %s %le",&d1,s1,s2,s3,&ns_add)!=EOF)
            if (strcmp(s2,"SiStressed")==0)
            {
                ns+=ns_add; // in cm-2 three times added
                a+=d1;
            }

        //Temperature =  10.0 K
        if (sscanf(s,"%s = %lf",s1,&V)!=EOF)
            if (strcmp(s1,"Temperature")==0)
            {
                T=V;
                kT=kb*T;
            }
        // "Electron eigenvalue 6 = 3.072603e-002 eV"
        if (sscanf(s,"%s %s %d = %le",s1,s2,&i,&V)==EOF) continue;
        if (strcmp(s1,"Electron")==0)
        {
            EigenCount=i;
            EigenValues[i]=V;
        }
    };
    printf("Temperature set to %lf K. \n",T);
    printf("%d EigenValues found. \n",EigenCount);
    printf("Free carriers concentration = %le cm^-2. \n",ns);
    printf("Well width = %le A. \n",a);
    ns*=10000; //cm-2->m-2
    a*=1e-10;//A->m
    fclose(stream);
}

void WriteStatus()
{
    int i;
    for(i=1;i<=EigenCount;i++)
        printf("SQ Position [%d]= %le eV\n",i,EigenValues[i]);
}

void ReadWF(char * WFFile)
{
    char s1 [1000];
    int i=0;
    double V;
    FILE *stream;
    stream = fopen(WFFile, "r");
    if (stream==NULL)
    {
        printf("Error! Wave functions file %s not found.\n",WFFile);
        return;
    }
    //Read header
    //Y (ang)	el Psi 1
    fscanf(stream,"Y (ang) %s Psi",s1);
    while ((!feof(stream))&&(strcmp(s1,"el")==0))
    {
        fscanf(stream,"%d %s Psi",&i,s1);
    }
    if (i!=EigenCount)
    {
        printf("Error! Number of subbands %d does not match expected value %d.",i,EigenCount);
        return;
    }
    fgets(s1,1000,stream); //rest of the header
    WFpointsCount=0;
    fscanf(stream,"%le", &V); i=0; //row count
    while (!feof(stream))
    {
        WFTable[i][0]=V;
        int j;// WF count
        for (j=1;j<=EigenCount;j++)
        {
            fscanf(stream,"%le", &WFTable[i][j]);
            //printf("WF[%d][%d]=%le \n",i,j,WFTable[i][j]);
        }
        fgets(s1,1000,stream); //holes out
        if (feof(stream))
        {
            i++;
            continue;
        }
        fscanf(stream,"%le", &V);
        i++;
    }
    if (i>WFpointsMAX)
    {
        printf("Error! Too many points for WF.");
        return;
    }
    WFpointsCount=i; //0..i-1
    printf("%d WFs %d points each found.\n",EigenCount,WFpointsCount);
}



double diffz(int j) //effective differential dz at the j-th point
{
    double dz;
    if (j==0) dz=WFTable[1][0]-WFTable[0][0]; else
                if (j==WFpointsCount) dz=WFTable[WFpointsCount-1][0]-WFTable[WFpointsCount-2][0];
                    else dz=(WFTable[j+1][0]-WFTable[j-1][0]);
    return dz/2;
}
double diffWF(int j, int i) //effective differential for the i-th WF at the j-th point
{
    double dWF;
    if (j==0) dWF=WFTable[1][i]-WFTable[0][i]; else
                if (j==WFpointsCount) dWF=WFTable[WFpointsCount-1][i]-WFTable[WFpointsCount-2][i];
                    else dWF=(WFTable[j+1][i]-WFTable[j-1][i]);
    return dWF/2;
}

void WFNormalize()
{
    int i,j;
    double s;
    printf("WF normalization.\n");

    for (i=1;i<=EigenCount;i++) //foreach subband
    {
        s=0;
        for (j=0;j<WFpointsCount;j++)
        {
            s+=pow(WFTable[j][i],2)*diffz(j);
        }
        s/=WFTable[WFpointsCount-1][0]-WFTable[0][0];
        for (j=0;j<WFpointsCount;j++)
            WFTable[j][i]/=sqrt(s);
    }
}

void CheckNormalization() //OK
{
    int i,j;
    double s;
    for (i=1;i<=EigenCount;i++) //foreach subband
    {
        s=0;
        for (j=0;j<WFpointsCount;j++)
        {
            s+=pow(WFTable[j][i],2)*diffz(j);
        }
        printf("Norm(WF[%d])=%le\n",i,s);
    }
}

double MatrixElement(int i1, int i2) //Vorobjev metod p.45
{
    double s; int j;
    s=0;
    if (i1==i2) return 0;
    for (j=0;j<WFpointsCount;j++)
        s+=WFTable[j][i1]*diffWF(j,i2);
    return -h_*s;
    //  /a
}

void CalcMatrixElements ()
{
    int i,j;
    for (i=1; i<=EigenCount;i++)
        for (j=1; j<=EigenCount;j++)
            if (i>j) MatrixElements[i][j]=-MatrixElements[j][i];
            else MatrixElements[i][j]=MatrixElement(i,j);
}
double OscStrength(int i1, int i2) //Oscillator strength Vorobjev metod p.45 3.19
{
    double d;
    d=2/mW/m0*pow(fabs(MatrixElements[i1][i2]),2)/(EigenValues[i2]-EigenValues[i1])/e;
    return d;
}

double SigmaIJ(int i1, int i2, double En, double Gamma) //energies should be in eV (or like EigenValues)
//absorption coefficient addition at EN energy from i1 to i2
// Lorenz broadening with Gamma
{
    double DeltaF; //delta Fermi
    double sigma, DeltaE, Lorenzian;
    if (i2<i1) {printf("ERROR! Target level %d should be above source one %d in absorption!\n",i2,i1);};
    if (i2>EigenCount) {printf("ERROR! Too many EigenLevels in absorption!");};
    DeltaF=Distribution(EigenValues[i1])-Distribution(EigenValues[i2]);
    DeltaE=EigenValues[i2]-EigenValues[i1];
    Lorenzian=1/PI*(Gamma/(pow(DeltaE-En,2)+Gamma*Gamma));
    sigma=OscStrength(i1,i2)*DeltaF*Lorenzian; //davies 8.94, 8.92 Coeff-?
    //sigma=sigma*T; // comes from concentrations when they are obtained by integrating distribution functions
    // correct for Bolzmann statistics only
    sigma*=ns;
    return sigma;
}

double Absorption(double En) //All intersubband, Energy in eV, Gamma static
{
    int i,j;
    double s;
    s=0;
    for (i=1;i<=EigenCount-1;i++)
        for (j=i+1;j<=EigenCount;j++)
            s=s+SigmaIJ(i,j,En,Gamma);
    return s;
}

double AbsorptionWithShift(double En) //All intersubband, Energy in eV, Gamma, allowed for Collective shift
{
    int i,j;
    double s;
    s=0;
    for (i=1;i<=EigenCount-1;i++)
        for (j=i+1;j<=EigenCount;j++)
            s=s+SigmaIJ(i,j,En-ASTable[i][j],Gamma);
    return s;
}

double rs(double nz) //Capasso p46
{
    return pow(4*PI/3*pow(aBeW,3)*nz,-1./3.);
}

double Vxc(double nz) //Exchange-correlation energy
{
    double RS;
    RS=rs(nz);
    return -pow(9*PI/4,1./3.)*2/PI/RS*(1+VxcB/VxcA*RS*log(1+VxcA/RS))*e*e/(8*PI*epsW*eps0*aBeW);
}

double dVpodnz(double z) //should be negative
{
    double nz,dnz;
    nz=ns/a; //!!!!!????? = const !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dnz=nz/10000;
    return (Vxc(nz+dnz)-Vxc(nz))/dnz;
}

double Beta(int i1, int i2)// Capasso 59 - exciton shift . units OK
{
    double s; int i;
    //jntegral here
    s=0;
    for (i=0;i<WFpointsCount;i++) //Capasso 58
    {
        s+=pow(WFTable[i][i1]*WFTable[i][i2],2)*diffz(i)*dVpodnz(WFTable[i][0]);
    }
    return 2*ns/(EigenValues[i2]-EigenValues[i1])/e*s; //eV -> J
}

double AbsorptionShift(int i1, int i2) //Absorption maximum shift due to collective effects, eV
{
    double alpha, beta, E21, AS, S11;
    double sum1, sum2, dz;
    int i;
    E21=fabs(EigenValues[i2]-EigenValues[i1]);

    sum1=0;sum2=0;
    for (i=0;i<WFpointsCount;i++) //Capasso 58
    {
        dz=diffz(i);
        sum2+=WFTable[i][i1]*WFTable[i][i2]*dz;
        sum1+=sum2*sum2*dz;
    }
    S11=sum1*1e-10; //Angstrom -> m

    alpha=-2*e*ns/epsW/eps0*S11/E21; //Capasso 58, E12 in eV

    beta=Beta(i1,i2);

    AS=E21*(pow(1+alpha-beta,0.5)-1); // Capasso 57
    return AS;
}

void FillASTable()
//Fill Collective Absorption Shift Table.
//symmetric
{
    int i,j;
    for (i=1;i<=EigenCount;i++) //by rows - from
        for (j=1;j<=EigenCount;j++) //by columns - to
            if (i==j) ASTable[i][j]=0;
            else ASTable[i][j]=AbsorptionShift(i,j);
}

double OrthoCheck(int i1, int i2)
{
    double s; int j;
    s=0;
    for (j=0;j<WFpointsCount;j++)
        s+=WFTable[j][i1]*WFTable[j][i2]*diffz(j);
    return s;
}

void OSOut(char* FileName)
//output oscillator strengths as a table separated by tabs
//row - from; column - to
{
    double r, OSums [EigenCountMAX+1];
    double podgon=-3.5e-5/a/a;
    FILE* stream;
    stream=fopen(FileName,"w");
    if (stream==NULL)
    {
        printf("Error! Can not open %s for writing.\n",FileName);
        return;
    }
    int i,j;
    for (i=1;i<=EigenCount;i++) OSums[i]=0;
    fprintf(stream,"From\\To\t"); //header
    for (i=1;i<=EigenCount;i++) fprintf(stream,"Sub%d\t",i);
    fprintf(stream,"\n"); //header
    for (i=1;i<=EigenCount;i++) //by rows - from
    {
        fprintf(stream,"Sub%d\t",i);
        for (j=1;j<=EigenCount;j++) //by columns - to
            if (i==j) fprintf(stream,"0\t");
            else
            {
                r=OscStrength(i,j);
                OSums[j]+=r;
                fprintf(stream,"%6.2lf\t",r/podgon);
            }
        fprintf(stream,"\n");
    }
    fprintf(stream,"\nSUM%d\t",i);
    for (j=1;j<=EigenCount;j++) //by columns - to
        fprintf(stream,"%lf\t",OSums[j]);

    fprintf(stream,"\n");
    fclose (stream);
}
void DeltaEOut(char* FileName)
//output delta energies in meV as a table separated by tabs
//row - from; column - to
{
    double r;
    FILE* stream;
    stream=fopen(FileName,"w");
    if (stream==NULL)
    {
        printf("Error! Can not open %s for writing.\n",FileName);
        return;
    }
    int i,j;
    fprintf(stream,"From\\To\t"); //header
    for (i=1;i<=EigenCount;i++) fprintf(stream,"Sub%d\t",i);
    fprintf(stream,"\n"); //header
    for (i=1;i<=EigenCount;i++) //by rows - from
    {
        fprintf(stream,"Sub%d\t",i);
        for (j=1;j<=EigenCount;j++) //by columns - to
            if (i==j) fprintf(stream,"0\t");
            else
            {
                r=(EigenValues[j]-EigenValues[i])*1000; //to meV
                fprintf(stream,"%6.1lf\t",r);
            }
        fprintf(stream,"\n");
    }

    fclose (stream);
}

void DeltaFOut(char* FileName)
//output delta Fermi function for all transitions as a table separated by tabs
//row - from; column - to
{
    double r;
    FILE* stream;
    stream=fopen(FileName,"w");
    if (stream==NULL)
    {
        printf("Error! Can not open %s for writing.\n",FileName);
        return;
    }
    int i,j;
    fprintf(stream,"From\\To\t"); //header
    for (i=1;i<=EigenCount;i++) fprintf(stream,"Sub%d\t",i);
    fprintf(stream,"\n"); //header
    for (i=1;i<=EigenCount;i++) //by rows - from
    {
        fprintf(stream,"Sub%d\t",i);
        for (j=1;j<=EigenCount;j++) //by columns - to
            if (i==j) fprintf(stream,"0\t");
            else
            {
                r=(Distribution(EigenValues[j])-Distribution(EigenValues[i]));
                fprintf(stream,"%lf\t",r);
            }
        fprintf(stream,"\n");
    }

    fclose (stream);
}

void AbsShiftOut(char* FileName)
//output collective absorption shift in meV as a table separated by tabs
//row - from; column - to
{
    double r;
    FILE* stream;
    stream=fopen(FileName,"w");
    if (stream==NULL)
    {
        printf("Error! Can not open %s for writing.\n",FileName);
        return;
    }
    int i,j;
    fprintf(stream,"From\\To\t"); //header
    for (i=1;i<=EigenCount;i++) fprintf(stream,"Sub%d\t",i);
    fprintf(stream,"\n"); //header
    for (i=1;i<=EigenCount;i++) //by rows - from
    {
        fprintf(stream,"Sub%d\t",i);
        for (j=1;j<=EigenCount;j++) //by columns - to
            if (i==j) fprintf(stream,"0\t");
            else
            {
                r=ASTable[i][j]*1000; //to meV
                fprintf(stream,"%6.1lf\t",r);
            }
        fprintf(stream,"\n");
    }

    fclose (stream);
}

void FullAbsOut(char* FileName, double En1, double En2)
//output full (summarized) absorption coefficient vs energy in meV
//arbitrary units
{
    double podgon=-1e9;
    int steps=10000;
    double de=(En2-En1)/steps;
    double En;
    FILE* stream;
    stream=fopen(FileName,"w");
    if (stream==NULL)
    {
        printf("Error! Can not open %s for writing.\n",FileName);
        return;
    }
    En=En1;
    fprintf(stream,"%10s %10s\n","E","Sigma"); //header
    while (En<=En2)
        {
        fprintf(stream,"%10le %10le\n",En*1000,podgon*Absorption(En));
        En+=de;
        }
    fclose (stream);
}

void FullShiftedAbsOut(char* FileName, double En1, double En2)
//output full (summarized) absorption coefficient vs energy in meV allowed for collective shift
//arbitrary units
{
    double podgon=-1e22;
    int steps=10000;
    double de=(En2-En1)/steps;
    double En;
    FILE* stream;
    stream=fopen(FileName,"w");
    if (stream==NULL)
    {
        printf("Error! Can not open %s for writing.\n",FileName);
        return;
    }
    En=En1;
    fprintf(stream,"%10s %10s\n","E","Sigma"); //header
    while (En<=En2)
        {
        fprintf(stream,"%10le %10le\n",En*1000,podgon*AbsorptionWithShift(En));
        En+=de;
        }
    fclose (stream);
}

int main()
{
    //int i1, i2;
    InitializeParameters();
    ReadStatus("myqw_St_Status.txt");
    //WriteStatus();
    ReadWF("myqw_St_Wave.txt");
    //CheckNormalization();

    CalcMatrixElements();
    OSOut("OS.txt");
    DeltaEOut("DeltaE.txt");
 //   DeltaFOut("DeltaF.txt");
    FillASTable();
    AbsShiftOut("AbsShift.txt");
    FullAbsOut("absorp.txt",0.0001,0.100);
    FullShiftedAbsOut("AbsPlus.txt",0.0001,0.100);
    return 0;
}
