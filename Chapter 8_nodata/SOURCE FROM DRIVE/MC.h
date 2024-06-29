void copycoordinates(long int idp)
{
     parts[idp].ox = parts[idp].x;
     parts[idp].oy = parts[idp].y;
     parts[idp].oz = parts[idp].z;

}

void rollback(long int idp)
{
     parts[idp].x = parts[idp].ox;
     parts[idp].y = parts[idp].oy;
     parts[idp].z = parts[idp].oz;

}


void translation()
{
    //random
    // (ran3(&mySys.seed);
}



void do_MC_sweep()
{
 
}

void do_MC(){

    char dumpname[100];
    char restartname[100];

    FILE* f = fopen("energy.dat", "a");
    FILE* g = fopen("acceptance.dat", "a"); 

    sprintf(restartname,"restartpoint.dat");
    
    for(mySys.step=0; mySys.step < mySys.NSteps; mySys.step++){
         
        do_MC_sweep();
        //if(mySys.step % 1000 == 0)  WriteConf(restartname);
     
        if(mySys.step % mySys.NPrint == 0){ 
            //printf("dumping...\n");
        }
    }
   
    fclose(f); fclose(g);
   
}

