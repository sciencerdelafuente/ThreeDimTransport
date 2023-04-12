/**************************************************************************************/  
//  Computes the two-layer map of particles and the three-dimensional singular values 
//  from the analytical ABC flow with an added settling velocity.
//  Related work and results: https://link.aps.org/doi/10.1103/PhysRevE.104.065111
/**************************************************************************************/  
//  Code written by Rebeca de la Fuente :  2020                                 
/**************************************************************************************/  
#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<cstdlib>
#include <cmath>
# define PI  3.141592653589793238462643383279502884L 


using namespace std;

long double RandomFloat(long double a, long double b); 

 
void jacobi (int n, long double **a, long double *d, long double **v);

int main(void)
{



     srand48(time(NULL));
     srand(time(NULL));



long double A,B,C;
long double k1,k2,k3,k4,k,l1,l2,l3,l4,l,s1,s2,s3,s4,s,xn,yn,zn,tn,h;
long double x11,y11,z11,x21,y21,z21,x12,y12,z12,x22,y22,z22,x13,y13,z13,x23,y23,z23;
long double T;
bool perturbation1,perturbation2,perturbation3,perturbation4;
long double min_x,min_y,min_z,max_x,max_y,max_z,angulo;
long double i_x,i_y,i_z,f_x,f_y,f_z;
long double dx,dy,dz;
int ksx,ksy,ksz;
long double long_x,long_y,long_z;
int box;
long double xn0,yn0,zn0,tn0;
long double det,traza,L1,L2,lll;
long double Lyapunov,Lyapunov_aux;
long double var_xx,var_yx,var_zx,var_xy,var_yy,var_zy,var_xz,var_yz,var_zz;
long double k5,k6,k7,k8,k9,k10,k11,k12;
long double l5,l6,l7,l8,l9,l10,l11,l12;
long double s5,s6,s7,s8,s9,s10,s11,s12;
long double f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12;
long double r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12;
long double RR[3][3],Tensor[3][3];
long double ent0,ent1,ent2,p,q,phy,ll,sigma,llM,llm;
long double L3;
long double vec[3],normal[3],vec2[3],ll2,vecM[3],vecm[3];
B=sqrt(2);
C=sqrt(3);
A=1;
min_x=0;
min_y=0;
min_z=0;
max_x=2*PI;
max_y=2*PI;
max_z=2*PI;
long_x=max_x-min_x;
long_y=max_y-min_y;
long_z=max_z-min_z;
tn=0;
    int n=3;
   long double **a, *d, **v;
 
   a = new long double* [n], d = new long double [n], v = new long double* [n];
 
   for(int j=0; j<n; j++){
 
      a[j] = new long double [n], v[j] = new long double [n];
    
   }  
long double tau0;
long double longitud;
int Tmaximum;
long double layer1,layer2;
long double c_aux;
long double xf,yf,zf;
int ksx1,ksz1,ksx2,ksz2,box2,ksy1,ksy2;
long double anterior_xn,anterior_yn,anterior_zn,posterior_xn,posterior_yn,posterior_zn;
long double lambda;
long double dx_boxes,kout,entropia;
long double tau1,tau2,tau3,tau4;
long double boxes=50;
long double number_boxes=boxes*boxes;
dx_boxes=long_x/boxes;
long double particles_per_box=70;
dx=dx_boxes/particles_per_box;
long double constante=-3.15;
long double vx,vy,vz,wx,wy,wz,sigma2_s,sigma2,medio;

 layer1=11;
 layer2=0;
 h=0.1;


long double distanciaa,distanciamax;
long double distancia_layers=fabs(layer2-layer1);
long double distancia_aux=layer2-max_y;
long double koutx,kouty,koutz;
long double auxiliar,auxiliar2,entropia2,Lyapunov_d,lll2;
long double tau_anterior,tau_posterior;
long double distancia_esfera,x_m,z_m;
long double vector_r,vector_s,absr,abss,angulov,Proj_s;


        long double moduloM,modulom,modulo,productoV;
        int total_particles=particles_per_box*particles_per_box;
        long double *cubex;
        cubex=(long double *)malloc(total_particles*sizeof(long double*));
        long double *cubey;
        cubey=(long double *)malloc(total_particles*sizeof(long double*));
        long double *cubez;
        cubez=(long double *)malloc(total_particles*sizeof(long double*));
        long double numero_mx0,numero_mx1,numero_my0,numero_my1,numero_mz0,numero_mz1;
        long double cube_mx0,cube_mx1,cube_my0,cube_my1,cube_mz0,cube_mz1;

        long double posterior_var_xx,posterior_var_yx,posterior_var_zx,posterior_var_xy,posterior_var_yy;
        long double posterior_var_zy,posterior_var_xz,posterior_var_yz,posterior_var_zz;
        long double anterior_var_xx,anterior_var_yx,anterior_var_zx,anterior_var_xy,anterior_var_yy;
        long double anterior_var_zy,anterior_var_xz,anterior_var_yz,anterior_var_zz;
        long double distancia_total,distancia_menor;

        long double VFxM,VFyM,VFzM,VFxm,VFym,VFzm,mx,my,mz,px,py,pz,qx,qy,qz;
        long double distanciacapa=fabs(layer2-layer1);
        long double xnBC,ynBC,znBC,elipse_menor;

        long double factor_numerador,factor_denominador,factor,LSM,LSm;
        long double mayor,projectionm,projectionAngulo;
    long double auxiliar3,projectionM;
    long double modulo1,modulo2,area;
    long double Tau=6;
    long double LyapunovM=0;
    long double Lyapunovn=0;
    long double Lyapunovm=0;


    long double ProyeccionMayorMedia,ProyeccionMenorMedia,ProyeccionMayor,ProyeccionMenor,AreaMedia;
    long double ZZ[2][2];
    long double xf_dx,yf_dx,xf_dy,yf_dy,zf_dx,zf_dy;
    long double perturbation=0.000001;


     double *L;
    L=(double *)malloc(number_boxes*sizeof(double*));

    double *Lx;
    Lx=(double *)malloc(boxes*sizeof(double*));

    double *Lz;
    Lz=(double *)malloc(boxes*sizeof(double*));

   long double ElipseMenorMedia,Koutx,Koutz,Lyapunovaux;

          int contando=0;
          long double radio,cota;

          radio=dx_boxes;
   long double proyeccion[2][2];
   long double Taumedia;
   long double dx1=(2*PI)/particles_per_box;
   long double dx2=radio/particles_per_box;
   long double projectiont,proyt,Lyapunovpeq,Lyapunovaux2;
   long double normal2[3];

   long double radio1,radio2,longitud1,longitud2;



    long double Px=1;
    long double Py=1;

   for(long double Px=min_x;Px<max_x;Px=Px+0.1)
  {
   for(long double Py=min_x;Py<max_x;Py=Py+0.1)
   {

    long double Pz=layer1;



   T=20;
  // for(int tiemposs=0;tiemposs<800;tiemposs++)
  // {

//    T=T+0.025;



                xn0=Px;
                yn0=Py;
                zn0=Pz;


 
         var_xx=1;
         var_yx=0;
         var_zx=0;
         var_xy=0;
         var_yy=1; 
         var_zy=0;
         var_xz=0;
         var_yz=0;
         var_zz=1;


         xn=xn0;
         yn=yn0;
         zn=zn0;
         tau0=0;
         for(long double tiempo=0;tiempo<T;tiempo=tiempo+h)
         {

            k1=A*sin(zn)+C*cos(yn);
            k2=B*sin(xn)+A*cos(zn);
            k3=C*sin(yn)+B*cos(xn)+constante;
            k4=-C*sin(yn)*var_xy+A*cos(zn)*var_xz;
            k5=-C*sin(yn)*var_yy+A*cos(zn)*var_yz;
            k6=-C*sin(yn)*var_zy+A*cos(zn)*var_zz;
            k7=B*cos(xn)*var_xx-A*sin(zn)*var_xz;
            k8=B*cos(xn)*var_yx-A*sin(zn)*var_yz;
            k9=B*cos(xn)*var_zx-A*sin(zn)*var_zz;
            k10=-B*sin(xn)*var_xx+C*cos(yn)*var_xy;
            k11=-B*sin(xn)*var_yx+C*cos(yn)*var_yy;
            k12=-B*sin(xn)*var_zx+C*cos(yn)*var_zy;

            l1=A*sin(zn+0.5*h*k3)+C*cos(yn+0.5*h*k2);
            l2=B*sin(xn+0.5*h*k1)+A*cos(zn+0.5*h*k3);
            l3=C*sin(yn+0.5*h*k2)+B*cos(xn+0.5*h*k1)+constante;
            l4=-C*sin(yn+0.5*h*k2)*(var_xy+0.5*h*k7)+A*cos(zn+0.5*h*k3)*(var_xz+0.5*h*k10);
            l5=-C*sin(yn+0.5*h*k2)*(var_yy+0.5*h*k8)+A*cos(zn+0.5*h*k3)*(var_yz+0.5*h*k11);
            l6=-C*sin(yn+0.5*h*k2)*(var_zy+0.5*h*k9)+A*cos(zn+0.5*h*k3)*(var_zz+0.5*h*k12);
            l7=B*cos(xn+0.5*h*k1)*(var_xx+0.5*h*k4)-A*sin(zn+0.5*h*k3)*(var_xz+0.5*h*k10);
            l8=B*cos(xn+0.5*h*k1)*(var_yx+0.5*h*k5)-A*sin(zn+0.5*h*k3)*(var_yz+0.5*h*k11);
            l9=B*cos(xn+0.5*h*k1)*(var_zx+0.5*h*k6)-A*sin(zn+0.5*h*k3)*(var_zz+0.5*h*k12);
            l10=-B*sin(xn+0.5*h*k1)*(var_xx+0.5*h*k4)+C*cos(yn+0.5*h*k2)*(var_xy+0.5*h*k7);
            l11=-B*sin(xn+0.5*h*k1)*(var_yx+0.5*h*k5)+C*cos(yn+0.5*h*k2)*(var_yy+0.5*h*k8);
            l12=-B*sin(xn+0.5*h*k1)*(var_zx+0.5*h*k6)+C*cos(yn+0.5*h*k2)*(var_zy+0.5*h*k9);

            s1=A*sin(zn+0.5*h*l3)+C*cos(yn+0.5*h*l2);
            s2=B*sin(xn+0.5*h*l1)+A*cos(zn+0.5*h*l3);
            s3=C*sin(yn+0.5*h*l2)+B*cos(xn+0.5*h*l1)+constante;
            s4=-C*sin(yn+0.5*h*l2)*(var_xy+0.5*h*l7)+A*cos(zn+0.5*h*l3)*(var_xz+0.5*h*l10);
            s5=-C*sin(yn+0.5*h*l2)*(var_yy+0.5*h*l8)+A*cos(zn+0.5*h*l3)*(var_yz+0.5*h*l11);
            s6=-C*sin(yn+0.5*h*l2)*(var_zy+0.5*h*l9)+A*cos(zn+0.5*h*l3)*(var_zz+0.5*h*l12);
            s7=B*cos(xn+0.5*h*l1)*(var_xx+0.5*h*l4)-A*sin(zn+0.5*h*l3)*(var_xz+0.5*h*l10);
            s8=B*cos(xn+0.5*h*l1)*(var_yx+0.5*h*l5)-A*sin(zn+0.5*h*l3)*(var_yz+0.5*h*l11);
            s9=B*cos(xn+0.5*h*l1)*(var_zx+0.5*h*l6)-A*sin(zn+0.5*h*l3)*(var_zz+0.5*h*l12);
            s10=-B*sin(xn+0.5*h*l1)*(var_xx+0.5*h*l4)+C*cos(yn+0.5*h*l2)*(var_xy+0.5*h*l7);
            s11=-B*sin(xn+0.5*h*l1)*(var_yx+0.5*h*l5)+C*cos(yn+0.5*h*l2)*(var_yy+0.5*h*l8);
            s12=-B*sin(xn+0.5*h*l1)*(var_zx+0.5*h*l6)+C*cos(yn+0.5*h*l2)*(var_zy+0.5*h*l9);


            f1=A*sin(zn+h*s3)+C*cos(yn+h*s2);
            f2=B*sin(xn+h*s1)+A*cos(zn+h*s3);
            f3=C*sin(yn+h*s2)+B*cos(xn+h*s1)+constante;
            f4=-C*sin(yn+h*s2)*(var_xy+h*s7)+A*cos(zn+h*s3)*(var_xz+h*s10);
            f5=-C*sin(yn+h*s2)*(var_yy+h*s8)+A*cos(zn+h*s3)*(var_yz+h*s11);
            f6=-C*sin(yn+h*s2)*(var_zy+h*s9)+A*cos(zn+h*s3)*(var_zz+h*s12);
            f7=B*cos(xn+h*s1)*(var_xx+h*s4)-A*sin(zn+h*s3)*(var_xz+h*s10);
            f8=B*cos(xn+h*s1)*(var_yx+h*s5)-A*sin(zn+h*s3)*(var_yz+h*s11);
            f9=B*cos(xn+h*s1)*(var_zx+h*s6)-A*sin(zn+h*s3)*(var_zz+h*s12);
            f10=-B*sin(xn+h*s1)*(var_xx+h*s4)+C*cos(yn+h*s2)*(var_xy+h*s7);
            f11=-B*sin(xn+h*s1)*(var_yx+h*s5)+C*cos(yn+h*s2)*(var_yy+h*s8);
            f12=-B*sin(xn+h*s1)*(var_zx+h*s6)+C*cos(yn+h*s2)*(var_zy+h*s9);


            r1=(k1+2*l1+2*s1+f1)/6;
            r2=(k2+2*l2+2*s2+f2)/6;
            r3=(k3+2*l3+2*s3+f3)/6;
            r4=(k4+2*l4+2*s4+f4)/6;
            r5=(k5+2*l5+2*s5+f5)/6;
            r6=(k6+2*l6+2*s6+f6)/6;
            r7=(k7+2*l7+2*s7+f7)/6;
            r8=(k8+2*l8+2*s8+f8)/6;
            r9=(k9+2*l9+2*s9+f9)/6;
            r10=(k10+2*l10+2*s10+f10)/6;
            r11=(k11+2*l11+2*s11+f11)/6;
            r12=(k12+2*l12+2*s12+f12)/6;

            xn=xn+h*r1;
            yn=yn+h*r2;
            zn=zn+h*r3;
            var_xx=var_xx+h*r4;
            var_yx=var_yx+h*r5;
            var_zx=var_zx+h*r6;
            var_xy=var_xy+h*r7;
            var_yy=var_yy+h*r8;
            var_zy=var_zy+h*r9;
            var_xz=var_xz+h*r10;
            var_yz=var_yz+h*r11;
            var_zz=var_zz+h*r12;
            tau0=tau0+h;



         }



          xf=xn;
          yf=yn;
          zf=zn;


            vx=A*sin(zn)+C*cos(yn);
            vy=B*sin(xn)+A*cos(zn);
            vz=C*sin(yn)+B*cos(xn)+constante;

               RR[0][0]=var_xx;
               RR[1][0]=var_yx;
               RR[2][0]=var_zx;

               RR[0][1]=var_xy;
               RR[1][1]=var_yy;
               RR[2][1]=var_zy;

               RR[0][2]=var_xz;
               RR[1][2]=var_yz;
               RR[2][2]=var_zz;


               //Green-Strain tensor
               for(int kk1=0;kk1<3;kk1++)
               {
                 for(int kk2=0;kk2<3;kk2++)
                 {
                    a[kk1][kk2]=(RR[0][kk1]*RR[0][kk2])+(RR[1][kk1]*RR[1][kk2])+(RR[2][kk1]*RR[2][kk2]);
                 }
               }
               jacobi(n, a, d, v); 


  if(d[0]>d[1] && d[0]>d[2]){llM=d[0];}
  if(d[1]>d[0] && d[1]>d[2]){llM=d[1];}
  if(d[2]>d[1] && d[2]>d[0]){llM=d[2];}
  if((d[0]>d[1] && d[0]<d[2]) || (d[0]<d[1] && d[0]>d[2])){ll=d[0];}
  if((d[1]>d[0] && d[1]<d[2]) || (d[1]<d[0] && d[1]>d[2])){ll=d[1];}
  if((d[2]>d[1] && d[2]<d[0]) || (d[2]<d[1] && d[2]>d[0])){ll=d[2];}
  if(d[0]<d[1] && d[0]<d[2]){llm=d[0];}
  if(d[1]<d[0] && d[1]<d[2]){llm=d[1];}
  if(d[2]<d[1] && d[2]<d[0]){llm=d[2];}


       llM=sqrt(llM);
       ll=sqrt(ll);
       llm=sqrt(llm);



         xn=xn0+perturbation;
         yn=yn0;
         zn=zn0;
         for(long double tiempo=0;tiempo<T;tiempo=tiempo+h)
         {

            k1=A*sin(zn)+C*cos(yn);
            k2=B*sin(xn)+A*cos(zn);
            k3=C*sin(yn)+B*cos(xn)+constante;

            l1=A*sin(zn+0.5*h*k3)+C*cos(yn+0.5*h*k2);
            l2=B*sin(xn+0.5*h*k1)+A*cos(zn+0.5*h*k3);
            l3=C*sin(yn+0.5*h*k2)+B*cos(xn+0.5*h*k1)+constante;

            s1=A*sin(zn+0.5*h*l3)+C*cos(yn+0.5*h*l2);
            s2=B*sin(xn+0.5*h*l1)+A*cos(zn+0.5*h*l3);
            s3=C*sin(yn+0.5*h*l2)+B*cos(xn+0.5*h*l1)+constante;

            f1=A*sin(zn+h*s3)+C*cos(yn+h*s2);
            f2=B*sin(xn+h*s1)+A*cos(zn+h*s3);
            f3=C*sin(yn+h*s2)+B*cos(xn+h*s1)+constante;


            r1=(k1+2*l1+2*s1+f1)/6;
            r2=(k2+2*l2+2*s2+f2)/6;
            r3=(k3+2*l3+2*s3+f3)/6;

            xn=xn+h*r1;
            yn=yn+h*r2;
            zn=zn+h*r3;
            tau0=tau0+h;
          }

          xf_dx=xn;
          yf_dx=yn;
          zf_dx=zn;

         xn=xn0;
         yn=yn0+perturbation;
         zn=zn0;
         tau0=0;
         for(long double tiempo=0;tiempo<T;tiempo=tiempo+h)
         {

            k1=A*sin(zn)+C*cos(yn);
            k2=B*sin(xn)+A*cos(zn);
            k3=C*sin(yn)+B*cos(xn)+constante;

            l1=A*sin(zn+0.5*h*k3)+C*cos(yn+0.5*h*k2);
            l2=B*sin(xn+0.5*h*k1)+A*cos(zn+0.5*h*k3);
            l3=C*sin(yn+0.5*h*k2)+B*cos(xn+0.5*h*k1)+constante;

            s1=A*sin(zn+0.5*h*l3)+C*cos(yn+0.5*h*l2);
            s2=B*sin(xn+0.5*h*l1)+A*cos(zn+0.5*h*l3);
            s3=C*sin(yn+0.5*h*l2)+B*cos(xn+0.5*h*l1)+constante;

            f1=A*sin(zn+h*s3)+C*cos(yn+h*s2);
            f2=B*sin(xn+h*s1)+A*cos(zn+h*s3);
            f3=C*sin(yn+h*s2)+B*cos(xn+h*s1)+constante;


            r1=(k1+2*l1+2*s1+f1)/6;
            r2=(k2+2*l2+2*s2+f2)/6;
            r3=(k3+2*l3+2*s3+f3)/6;

            xn=xn+h*r1;
            yn=yn+h*r2;
            zn=zn+h*r3;
            tau0=tau0+h;
          }


          xf_dy=xn;
          yf_dy=yn;
          zf_dy=zn;
       


       vecM[0]=((xf_dx-xf)/perturbation);
       vecM[1]=((yf_dx-yf)/perturbation);
       vecM[2]=((zf_dx-zf)/perturbation);

       vec[0]=((xf_dy-xf)/perturbation);
       vec[1]=((yf_dy-yf)/perturbation);
       vec[2]=((zf_dy-zf)/perturbation);


               RR[0][0]=vecM[0];
               RR[1][0]=vecM[1];
               RR[2][0]=vecM[2];

               RR[0][1]=vec[0];
               RR[1][1]=vec[1];
               RR[2][1]=vec[2];


               a[0][0]=RR[0][0]*RR[0][0]+RR[1][0]*RR[1][0]+RR[2][0]*RR[2][0];
               a[0][1]=RR[0][0]*RR[0][1]+RR[1][0]*RR[1][1]+RR[2][0]*RR[2][1];
               a[1][0]=a[0][1];
               a[1][1]=RR[1][1]*RR[1][1]+RR[0][1]*RR[0][1]+RR[2][1]*RR[2][1];


               traza=a[0][0]+a[1][1];
               det=a[0][0]*a[1][1]-a[1][0]*a[0][1];
               L1=(traza+sqrt(traza*traza-4*det))*0.5;
               L2=(traza-sqrt(traza*traza-4*det))*0.5;
               L1=(L1);
               L2=(L2);
               LSM=L1;
               LSm=L2;
               if(L2>L1)
               {
                  LSM=L2;
                  LSm=L1;
               }


  LSM=sqrt(LSM);
  LSm=sqrt(LSm);




       LSM=log(LSM)/T;
       LSm=log(LSm)/T;
       llM=log(llM)/T;
       ll=log(ll)/T;
       llm=log(llm)/T;



  
          cout<<T<<" "<<llM<<" "<<ll<<" "<<llm<<"        "<<LSM<<" "<<LSm<<endl; //"     "<<llM+ll+llm<<"   "<<LSM+LSm<<endl;
 // cout<<LSM+LSm<<"   "<<llM+ll<<endl;

       cout<<LSM<<"  "<<llM<<endl;
  
         
       // }


    }}
}























void jacobi (int n, long double **a, long double *d, long double **v) {
 
   // Define las variables de tipo entero
 
   int i, j, ip, iq, nrot;
 
// Define las variables de tipo doble
 
   long double *b, *z;
 
   b = new long double [n];
 
   z = new long double [n];
 
   b = new long double [n]; z = new long double [n];
 
   long double c, g, h, s, sm, t, tau, theta, tresh;
 
// Inicializa a la matriz identidad
 
   for (ip = 0; ip < n; ip++) { 
 
      for (iq = 0; iq < n; iq++) {
 
         v[ip][iq] = 0;
 
      }
     
      v[ip][ip] = 1;
 
   }
 
// Inicializa b y d a la diagonal de a
 
   for (ip = 0; ip < n; ip++) {   
 
      b[ip] = a[ip][ip];
 
      d[ip] = b[ip];
 
      z[ip] = 0;
 
   }
 
   nrot = 0;
 
   for (i = 0; i < 50; i++) {
 
      sm = 0;
     
      for (ip = 0; ip < n - 1; ip++) {
 
         for (iq = ip + 1; iq < n; iq++) {
 
            sm +=fabs(a[ip][iq]);
 
         }
 
      } 
 
      if (sm == 0) break;
     
      if (i < 4)
     
         tresh = 0.2*sm/(n*n);
 
      else
 
         tresh = 0.0;
 
      for (ip =0; ip < n -1; ip++) {
 
         for (iq = ip + 1; iq < n; iq++) {
 
            g = 100.0*fabs(a[ip][iq]);
 
            if(i>4 && (long double)(fabs(d[ip])+g) == (long double)fabs(d[ip])
 
               && (long double)(fabs(d[iq])+g) == (long double)fabs(d[iq]))
 
               a[ip][iq] = 0.0;
 
            else if (fabs(a[ip][iq]) > tresh) {
 
               h = d[iq] - d[ip];
 
               if ((long double)(fabs(h)+g) == (long double)fabs(h))
 
                  t = (a[ip][iq])/h;   // t = 1/(2theta)
 
               else {
 
                  theta = 0.5*h/(a[ip][iq]);
 
                  t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
 
                  if(theta < 0.0) t = -t;
 
               }
 
                  c = 1.0/sqrt(1+t*t);
 
                  s = t*c;
 
                  tau = s/(1.0+c);
 
                  h = t*a[ip][iq];
 
                  z[ip] -=h;
 
                  z[iq] +=h;
 
                  d[ip] -=h;
 
                  d[iq] +=h;
 
                  a[ip][iq] = 0.0;
 
// Varía desde 0 hasta  ip - 1
 
               for (j =0; j < ip; j++) { 
 
                  g = a[j][ip];
 
                  h = a[j][iq];
 
                  a[j][ip] = g - s*(h+g*tau);
 
                  a[j][iq] = h + s*(g-h*tau);
 
               } 
 
// Varía desde ip+1 hasta  iq - 1
 
               for (j =ip+1; j < iq; j++) { 
 
                  g = a[ip][j];
 
                  h = a[j][iq];
 
                  a[ip][j] = g - s*(h+g*tau);
 
                  a[j][iq] = h + s*(g-h*tau);
 
               } 
 
               for (j =iq+1; j < n; j++) { 
 
                  g = a[ip][j];
 
                  h = a[iq][j];
 
                  a[ip][j] = g - s*(h+g*tau);
 
                  a[iq][j] = h + s*(g-h*tau);
 
               } 
 
               
               for (j =0; j < n; j++) { 
 
                  g = v[j][ip];
 
                  h = v[j][iq];
 
                  v[j][ip] = g - s*(h+g*tau);
 
                  v[j][iq] = h + s*(g-h*tau);
 
 
                  } 
 
               ++(nrot);
 
            }
 
         }
 
      }
 
         for (ip = 0; ip < n; ip++) {
 
            b[ip] = b[ip]+z[ip];
 
            d[ip] = b[ip];
 
            z[ip] = 0.0;
 
         }
       
   }   
 
   delete [] b, z;
 
}


long double RandomFloat(long double a, long double b) 
{
    long double random = ((long double) rand()) / (long double) RAND_MAX;
    long double diff = b - a;
    long double r = random * diff;
    return a + r;
}
