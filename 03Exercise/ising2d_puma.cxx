//
//  main.cpp
//  Ising
//
//  Created by Paolo Umari on 13/12/19.
//  Copyright (c) 2019 unipd. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;


double energia(double Jconst, int nx, int ny, double* spins)
{
    double ene=0.;
    int icount=0;
    for(int j=0;j<ny;j++){
        for( int i=0; i<nx;i++){
            
            int pv_x,pv_y,pv_count;
            
            if(i<nx-1){
                pv_x=i+1;
            }else{
                pv_x=0;
            }
            pv_y=j;
            pv_count=pv_y*nx+pv_x;
            ene+=-0.5*Jconst*spins[icount]*spins[pv_count];
            
            if(i>0){
                pv_x=i-1;
            }else{
                pv_x=nx-1;
            }
            pv_y=j;
            pv_count=pv_y*nx+pv_x;
            ene+=-0.5*Jconst*spins[icount]*spins[pv_count];
            
            
            if(j<ny-1){
                pv_y=j+1;
            }else{
                pv_y=0;
            }
            pv_x=i;
            pv_count=pv_y*nx+pv_x;
            ene+=-0.5*Jconst*spins[icount]*spins[pv_count];
            
            if(j>0){
                pv_y=j-1;
            }else{
                pv_y=ny-1;
            }
            pv_x=i;
            pv_count=pv_y*nx+pv_x;
            ene+=-0.5*Jconst*spins[icount]*spins[pv_count];

            
           
            icount++;
       }
        
    }
    
    return ene;
}


double numero_random( void){
    const long int a=16807;
    const long int c=0;
    const long m=2147483647;
    static long int x0=1;
    long int x1;
    double r;
    
    x1=(a*x0+c)%m;
    r=((double) x1)/((double) m);
    x0=x1;
    
    return r;
    
}

void genera_configurazione(int nx, int ny, double *spins)
{
    double r;
    int icount=0;
    for(int j=0;j<ny;j++){
        for( int i=0; i<nx;i++){
            r =numero_random();
            if(r<0.5){
                spins[icount]=-1;
            }else{
                spins[icount]=1;
            }
            icount++;
        }
    }
    
}


double cambia_configurazione(double Jconst,int nx, int ny, double *spins, double* spins1)
{
    for(int i=0;i< nx*ny;i++){
        spins1[i]=spins[i];
    }
    double r;
    int i,j;
    int icount;
    double ene=0.;
    
    r=numero_random();
    icount=floor(r*(nx*ny));
  //  std::cout << " Cambia" << icount << '\n' ;
    spins1[icount]*=-1.;
    j=floor((double)icount/(double)nx);
    i=icount-j*nx;
    
    int pv_x,pv_y,pv_count;
    
    if(i<nx-1){
        pv_x=i+1;
    }else{
        pv_x=0;
    }
    pv_y=j;
    pv_count=pv_y*nx+pv_x;
    ene+=-Jconst*spins1[icount]*spins1[pv_count];
    ene-=-Jconst*spins[icount]*spins[pv_count];
    
    if(i>0){
        pv_x=i-1;
    }else{
        pv_x=nx-1;
    }
    pv_y=j;
    pv_count=pv_y*nx+pv_x;
    ene+=-Jconst*spins1[icount]*spins1[pv_count];
    ene-=-Jconst*spins[icount]*spins[pv_count];
    
    
    if(j<ny-1){
        pv_y=j+1;
    }else{
        pv_y=0;
    }
    pv_x=i;
    pv_count=pv_y*nx+pv_x;
    ene+=-Jconst*spins1[icount]*spins1[pv_count];
    ene-=-Jconst*spins[icount]*spins[pv_count];
    
    if(j>0){
        pv_y=j-1;
    }else{
        pv_y=ny-1;
    }
    pv_x=i;
    pv_count=pv_y*nx+pv_x;
    ene+=-Jconst*spins1[icount]*spins1[pv_count];
    ene-=-Jconst*spins[icount]*spins[pv_count];
    

    return ene;
    
}

double magnetizzazione(int nx, int ny, double *spins)
{   int icount=0;
    double ma=0;
    for(int j=0;j<ny;j++){
        for( int i=0; i<nx;i++){
            ma+=spins[icount];
            icount++;
        }
    }
    return ma;
}

int main(int argc, const char * argv[])
{
    int nx,ny;
    long int nsteps;
    double J,kbT;
    double ene_ave=0., ene_var=0.;
    double ma_ave=0., ma_var=0.;
    double masq_ave=0., enesq_ave=0.;
    
    double ene0,ene1;
    double ratio,r;
    int metodo;
    int icount=0;
    double diff_ene;
    
    std::cout << "Nx :\n";
    std::cin >> nx;
    std::cout << "Nx :\n";
    std::cin >> ny;
    std::cout << "J :\n";
    std::cin >> J;
    std::cout << "kbT :\n";
    std::cin >> kbT;
    std::cout << "N passi :\n";
    std::cin >> nsteps;
    std::cout << "1-un spin alla volta 2-tuttu i :\n";
    std::cin >> metodo;

    
    std::cout.precision(10);
    
    
    double* ene = new double[nsteps];
    double* ma = new double[nsteps];
    double* spins = new double[nx*ny];
    double* spins1 = new double[nx*ny];
    
    genera_configurazione( nx,  ny, spins);
    ene0=energia( J,  nx,  ny, spins);
    ene[0]=ene0;
    ma[0]=magnetizzazione( nx,  ny, spins);
    
    for(long int n=0; n <nsteps; n++){
        if(n%1000==0) std::cout << "Passo :" << n <<'\n';
        
        if(metodo==1){
             diff_ene=cambia_configurazione( J,nx,  ny, spins,spins1);
        }else{
              genera_configurazione( nx,  ny, spins1);
            diff_ene=energia( J,  nx,  ny, spins1)-ene0;
        }
           // ene1=energia( J,  nx,  ny, spins1);
          //  ene0=energia( J,  nx,  ny, spins);
           
         //  std::cout << "Diff ene"  <<ene0-ene1 <<'\n';
            ratio=exp((-diff_ene)/kbT);
            if(ratio>1.0){
                ene0+=diff_ene;;
                icount++;
                ene[icount]=ene0;
                ma[icount]=magnetizzazione( nx,  ny, spins);
               // std::cout << "Energia > 1  :" << ene0 <<'\n';
                for(int i=0;i< nx*ny;i++){
                    spins[i]=spins1[i];
                }
            }else{
                r=numero_random( );
              //  std::cout << "r ratio :" << r << ' ' << ratio <<'\n';
                if(r<ratio){
               //     std::cout << "r ratio :" << r << ' ' << ratio <<'\n';
                    ene0+=diff_ene;
                    icount++;
                    ene[icount]=ene0;
                    ma[icount]=magnetizzazione( nx,  ny, spins);
                //    std::cout << "Energia  < 1  :" << ene0 <<'\n';
                    for(int i=0;i< nx*ny;i++){
                        spins[i]=spins1[i];
                    }
                    
                }else{
                    icount++;
                    ene[icount]=ene0;
                    ma[icount]=magnetizzazione( nx,  ny, spins);
                    
                }
            
            
            
            }
    }
    
    //calcola medie e varianza
    ofstream filemc;
    filemc.open("montecarlo.dat",ios::out);
    filemc.precision(10);

    
    for(long int n=0; n <nsteps; n++){
        ene_ave+=ene[n];
        ma_ave+=ma[n];
        enesq_ave+=pow(ene[n],2);
        masq_ave+=pow(ma[n],2);
        if(n%100==0) filemc << n << "  " << ene[n]/(nx*ny) << "  " << ma[n]/(nx*ny) << '\n';
    }
                                 
    filemc.close();
    
    ene_ave/=(nsteps*nx*ny);
    enesq_ave/=(nsteps*nx*ny);
    ma_ave/=(nsteps*nx*ny);
    masq_ave/=(nsteps*nx*ny);
    
    std::cout << "Energia per sito : "<< ene_ave << "  +/-  " << sqrt((enesq_ave-pow(ene_ave,2))/nsteps) << '\n';
    std::cout << "Magnetizzazione totale : "<< ma_ave << "  +/-  " << sqrt((masq_ave-pow(ma_ave,2))/nsteps) << '\n';
    
    delete( ene);
    delete(ma);
    delete (spins);
    delete (spins1);
   
    
    return 0;
}


