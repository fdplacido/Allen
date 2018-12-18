#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

/*
  Author: Pierre Billoir
  Date: 01/2018
  
 */
using namespace std;

#define sqr(x) ((x)*(x))

static constexpr int NBINXMAX = 100;
static constexpr int  NBINYMAX = 100;
static constexpr int NSTEP = 50;
static constexpr int TWODIPOLES = 0;


class Coef 
{
 public:
  int Degx1,Degx2,Degy1,Degy2;
  vector<double> x00,x10,x01,tx00,tx10,tx01,y00,y10,y01,ty00,ty10,ty01;
  Coef() {};
  ~Coef() {};
  int Read(FILE *inp, int degx1, int degx2, int degy1, int degy2);
};
Coef operator+(Coef a, Coef b);

Coef operator-(Coef a, Coef b);

Coef operator*(Coef a, double p);

struct parameters{ 
  double ZINI,ZFIN,PMIN,BEND,Txmax,Tymax,Xmax,Ymax,Dtxy,step;
  int Nbinx,Nbiny,XGridOption,YGridOption,QuadraticInterpolation,DEGX1,DEGX2,DEGY1,DEGY2;
  Coef C[NBINXMAX][NBINYMAX];
};

/*************************************************************************/



/*---------------------------------------------------------------------------------------------*/
void ReadCoef(char *name, parameters& params);

/*-------------------------------------------------------------------------------------------------------------------------*/
int extrap(const double zi,const double zf,const float xi,const float yi,const float txi,const float tyi,const float qop,const double bend,const int quad_interp,float& xf,float& yf,float& txf,float& tyf, const parameters params);


