#include "Parametrization.h"
#include <fstream>     // std::cout, std::ios

/*
  Author: Pierre Billoir
  Date: 01/2018
  
 */


int Coef::Read(FILE *inp, int degx1, int degx2, int degy1, int degy2) {
    double a;
    Degx1 = degx1; Degx2 = degx2; Degy1 = degy1; Degy2 = degy2;
    for(int i=0; i<degx2; i++) { fscanf(inp,"%lf",&a); x00.push_back(a); }
    for(int i=0; i<degx2; i++) { fscanf(inp,"%lf",&a); tx00.push_back(1e-3*a); }
    for(int i=0; i<degx1; i++) { fscanf(inp,"%lf",&a); x10.push_back(a); }
    for(int i=0; i<degx1; i++) { fscanf(inp,"%lf",&a); x01.push_back(a); }
    for(int i=0; i<degx1; i++) { fscanf(inp,"%lf",&a); tx10.push_back(1e-3*a); }
    for(int i=0; i<degx1; i++) { fscanf(inp,"%lf",&a); tx01.push_back(1e-3*a); }

    for(int i=0; i<degy2; i++) { fscanf(inp,"%lf",&a); y00.push_back(a); }
    for(int i=0; i<degy2; i++) { fscanf(inp,"%lf",&a); ty00.push_back(1e-3*a); }
    for(int i=0; i<degy1; i++) { fscanf(inp,"%lf",&a); y10.push_back(a); }
    for(int i=0; i<degy1; i++) { fscanf(inp,"%lf",&a); y01.push_back(a); }
    for(int i=0; i<degy1; i++) { fscanf(inp,"%lf",&a); ty10.push_back(1e-3*a); }
    for(int i=0; i<degy1; i++) { fscanf(inp,"%lf",&a); ty01.push_back(1e-3*a); }
    if(feof(inp)) return 0;
    return 1;
  }

Coef operator+(Coef a, Coef b)
{
  Coef c = a;
  //c.Degx1 = a.Degx1; c.Degx2 = a.Degx2; c.Degy1 = a.Degy1; c.Degy2 = a.Degy2;
  for(int i=0; i<c.Degx2; i++) { c.x00[i] += b.x00[i]; c.tx00[i] += b.tx00[i]; }
  for(int i=0; i<c.Degx1; i++) { c.x10[i] += b.x10[i]; c.x01[i] += b.x01[i]; c.tx10[i] += b.tx10[i]; c.tx01[i] += b.tx01[i]; }
  for(int i=0; i<c.Degy2; i++) { c.y00[i] += b.y00[i]; c.ty00[i] += b.ty00[i]; }
  for(int i=0; i<c.Degy1; i++) { c.y10[i] += b.y10[i]; c.y01[i] += b.y01[i]; c.ty10[i] += b.ty10[i]; c.ty01[i] += b.ty01[i]; }
  return c;
}

Coef operator-(Coef a, Coef b)
{
  Coef c = a;
  //c.Degx1 = a.Degx1; c.Degx2 = a.Degx2; c.Degy1 = a.Degy1; c.Degy2 = a.Degy2;
  for(int i=0; i<c.Degx2; i++) { c.x00[i] -= b.x00[i]; c.tx00[i] -= b.tx00[i]; }
  for(int i=0; i<c.Degx1; i++) { c.x10[i] -= b.x10[i]; c.x01[i] -= b.x01[i]; c.tx10[i] -= b.tx10[i]; c.tx01[i] -= b.tx01[i]; }
  for(int i=0; i<c.Degy2; i++) { c.y00[i] -= b.y00[i]; c.ty00[i] -= b.ty00[i]; }
  for(int i=0; i<c.Degy1; i++) { c.y10[i] -= b.y10[i]; c.y01[i] -= b.y01[i]; c.ty10[i] -= b.ty10[i]; c.ty01[i] -= b.ty01[i]; }
  return c;
}

Coef operator*(Coef a, double p)
{
  Coef c = a;
  //c.Degx1 = a.Degx1; c.Degx2 = a.Degx2; c.Degy1 = a.Degy1; c.Degy2 = a.Degy2;
  for(int i=0; i<c.Degx2; i++) { c.x00[i] *= p; c.tx00[i] *= p; }
  for(int i=0; i<c.Degx1; i++) { c.x10[i] *= p; c.x01[i] *= p; c.tx10[i] *= p; c.tx01[i] *= p; }
  for(int i=0; i<c.Degy2; i++) { c.y00[i] *= p; c.ty00[i] *= p; }
  for(int i=0; i<c.Degy1; i++) { c.y10[i] *= p; c.y01[i] *= p; c.ty10[i] *= p; c.ty01[i] *= p; }
  return c;
}

void ReadCoef(char *name, parameters& params)
{
  FILE *coef;
  coef = fopen(name,"r");
  
  fscanf(coef,"%lf %lf %lf %lf %lf %lf %lf",&params.ZINI,&params.ZFIN,&params.PMIN,&params.BEND,&params.Txmax,&params.Tymax,&params.Dtxy);
  fscanf(coef,"%d %d %d %d %d %d %d %d",&params.Nbinx,&params.Nbiny,&params.XGridOption,&params.YGridOption,&params.DEGX1,&params.DEGX2,&params.DEGY1,&params.DEGY2);
  printf("ZINI,ZFIN %5.0f %5.0f\n",params.ZINI,params.ZFIN);
  printf("PMIN %5.0f\n",params.PMIN);
  printf("Nbinx,Nbiny %d %d\n",params.Nbinx,params.Nbiny);
  printf("GridOptions %d %d\n",params.XGridOption,params.YGridOption);
  printf("DEGX1,DEGX2 %d %d   DEGY1,DEGY2 %d %d\n",params.DEGX1,params.DEGX2,params.DEGY1,params.DEGY2);
  printf("Dtxy %f\n",params.Dtxy);
  printf("bending at ZINI  %5.1f rad.MeV\n",params.BEND);

  for(int ix=0; ix<params.Nbinx; ix++) for(int iy=0; iy<params.Nbiny; iy++) params.C[ix][iy].Read(coef,params.DEGX1,params.DEGX2,params.DEGY1,params.DEGY2);
}

int extrap(const double zi,const double zf,const float xi,const float yi,const float txi,const float tyi,const float qop,const double bend,const int quad_interp,float& xf,float& yf,float& txf,float& tyf, const parameters params)
// extrapolation from plane zi to plane zf, from initial state (xi,yi,txi,tyi,qop) to final state (xf,yf,txf,tyf)
// the bending from origin to zi is approximated by adding bend*qop to xi/zi.
// quad_inperp (logical): if true, the quadratic interpolation is used (better, with a little bit more computations)
// XGridO[tion and YGridOption describe the choice of xy grid. By default, it is 1 (equally spaced values)   
{
  float xx,yy,dx,dy,ux,uy;
  int ix,iy;
  if(fabs(xi)>params.Xmax||fabs(yi)>params.Ymax) return 0;
  switch(params.XGridOption) {
    case 1: xx = xi/params.Xmax; break;
    case 2: xx = sqr(xi/params.Xmax); if(xi<0) xx = -xx; break;
    case 3: xx = xi/params.Xmax; xx = xx*xx*xx; break;
    case 4: xx = asin(xi/params.Xmax)*2/M_PI; break;
  }
  switch(params.YGridOption) {
    case 1: yy = yi/params.Ymax; break;
    case 2: yy = sqr(yi/params.Ymax); if(yi<0) yy = -yy; break;
    case 3: yy = yi/params.Ymax; yy = yy*yy*yy; break;
    case 4: yy = asin(yi/params.Ymax)*2/M_PI; break;
  }
  dx = params.Nbinx*(xx+1)/2; ix = dx; dx -= ix;
  dy = params.Nbiny*(yy+1)/2; iy = dy; dy -= iy;

  ux = (txi-xi/zi-bend*qop)/params.Dtxy; uy = (tyi-yi/zi)/params.Dtxy;
  if(fabs(ux)>2||fabs(uy)>2) return 0;

  Coef c;

  if(quad_interp) {
    float gx,gy;
    gx = dx-.5; gy = dy-.5;
    //if(gx*gx+gy*gy>.01) return 0;
    if(ix==0) { ix = 1; gx -= 1.; }
    if(ix==params.Nbinx-1) { ix = params.Nbinx-2; gx += 1.; }
    if(iy==0) { iy = 1; gy -= 1.; }
    if(iy==params.Nbiny-1) { iy = params.Nbiny-2; gy += 1.; }

    int rx,ry,sx,sy;
    rx = (gx>=0); sx = 2*rx-1; ry = (gy>=0); sy = 2*ry-1;
    Coef c00,cp0,c0p,cn0,c0n,cadd;
    c00 = params.C[ix][iy]; cp0 = params.C[ix+1][iy]; c0p = params.C[ix][iy+1]; c0n = params.C[ix][iy-1]; cn0 = params.C[ix-1][iy];
    cadd = params.C[ix+sx][iy+sy];
    float gxy = gx*gy, gx2 = gx*gx, gy2 = gy*gy, g2 = gx*gx+gy*gy;
    c = c00*(1-g2) + (cp0*(gx2+gx) + cn0*(gx2-gx) + c0p*(gy2+gy) + c0n*(gy2-gy))*.5
      + ((c00+cadd)*sx*sy - cp0*rx*sy + cn0*(!rx)*sy - c0p*ry*sx + c0n*(!ry)*sx)*gxy;
  }
  else {
    float ex,fx,ey,fy;
    int jx,jy;
    if(dx<.5) { jx = ix-1; ex = .5+dx; fx = .5-dx; } else { jx = ix+1; ex = 1.5-dx; fx = dx-.5; }
    if(dy<.5) { jy = iy-1; ey = .5+dy; fy = .5-dy; } else { jy = iy+1; ey = 1.5-dy; fy = dy-.5; }
    if(ix<0||ix>=params.Nbinx||iy<0||iy>=params.Nbiny || jx<0||jx>=params.Nbinx||jy<0||jy>=params.Nbiny) return 0;
    Coef c_ii = params.C[ix][iy], c_ij = params.C[ix][jy], c_ji = params.C[jx][iy], c_jj = params.C[jx][jy];
    //printf("x00 %f %f %f %f\n",c_ii.x00[0],c_ij.x00[0],c_ji.x00[0],c_jj.x00[0]);
    c = c_ii*ex*ey + c_ij*ex*fy + c_ji*fx*ey + c_jj*fx*fy;
  }
  // straight line
  xf = xi+txi*(zf-zi);
  yf = yi+tyi*(zf-zi);
  txf = txi;
  tyf = tyi;
  // corrections to straight line
  float fq = qop*params.PMIN;
  float ff = 1.f;
  for(int deg=0; deg<c.Degx2; deg++) {
    ff *= fq;
    xf  += c.x00[deg]*ff;
    txf += c.tx00[deg]*ff;
    if(deg>=c.Degx1) continue;
    xf  += ( c.x10[deg]*ux +  c.x01[deg]*uy)*ff;
    txf += (c.tx10[deg]*ux + c.tx01[deg]*uy)*ff;
  }
  ff = 1.f;
  for(int deg=0; deg<c.Degy2; deg++) {
    ff *= fq;
    yf  +=  c.y00[deg]*ff;
    tyf += c.ty00[deg]*ff;
    if(deg>=c.Degy1) continue;
    yf  += ( c.y10[deg]*ux +  c.y01[deg]*uy)*ff;
    tyf += (c.ty10[deg]*ux + c.ty01[deg]*uy)*ff;
  }
  return 1;
}

/**************************************************************************************/

/************************************************************************************************************************/

// #include "Field_RK_eloss.cc"

// int main(int argc, char *argv[])
// {
//   FILE *req, *inp, *res;
//   char name_coef[200], name_part[200], name_res[200];

//   req = fopen("use_param.req","r");
//   fscanf(req,"%*s %s %s %s",name_coef,name_part,name_res);
//   printf("coef_file %s   part_file %s\n",name_coef,name_part);
//   fscanf(req,"%*s %d",&QuadraticInterpolation);
//   printf("QuadraticInterpolation %d\n",QuadraticInterpolation);

//   ReadCoef(name_coef);
//   inp = fopen(name_part,"r");
//   res = fopen(name_res,"w");
//   ReadField();

//   // reading the particle file 
//   Txmax = Tymax = .25; Xmax = ZINI*Txmax; Ymax = ZINI*Tymax;
//   int code,signe;
//   float px,py,pz,vx,vy,vz,tx,ty,qop,qopi,qopf;

//   for(int itry=0; itry<100000; itry++) {
//     fscanf(inp,"%d  %lf %lf %lf  %lf %lf %lf  %d",&code,&px,&py,&pz,&vx,&vy,&vz,&signe);
//     if(feof(inp)) break;
//     qop = signe/sqrt(px*px+py*py+pz*pz);
//     if(vz>ZINI||1/fabs(qop)<PMIN) continue;
//     tx = px/pz; ty = py/pz;
//     if(fabs(tx)>Txmax||fabs(ty)>Tymax) continue;

//     float xi,yi,txi,tyi,xr,yr,txr,tyr,xf,yf,txf,tyf;
//     // propagation par RK4 de vz a ZINI
//     nstep_RK4(20,vz,(ZINI-vz)/20,vx,vy,tx,ty,qop,0,xi,yi,txi,tyi,qopi);
//     // propagation par RK4 et extrap() de ZINI a ZFIN
//     float eloss = 2./3e5;
//     nstep_RK4(NSTEP,ZINI,(ZFIN-ZINI)/NSTEP,xi,yi,txi,tyi,qopi,eloss,xr,yr,txr,tyr,qopf);
//     if(extrap(ZINI,ZFIN,xi,yi,txi,tyi,qop,BEND,QuadraticInterpolation,xf,yf,txf,tyf))
//       fprintf(res,"%7.4f %7.4f %7.4f   %7.4f %7.4f %7.4f %7.4f\n",1e3*qop,tx,ty,xf-xr,yf-yr,1e3*(txf-txr),1e3*(tyf-tyr));
//   }

//   return 1;
// }
