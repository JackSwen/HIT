#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void count(double *dat, int *zone, int dat_row, int zxindex, int zyindex, int zzindex,double Xori, double Yori, double Zori, double cubesize){
for (int i = 0; i<dat_row; i++){
 double x = *(dat+i*3)-(Xori);
 double y = *(dat+i*3+1)-(Yori);
 double z = *(dat+i*3+2)-(Zori);
 int zx = x/(double)cubesize;
 int zy = y/(double)cubesize;
 int zz = z/(double)cubesize;
 (*(zone+(zx)*zyindex*zzindex+(zy)*zzindex+(zz)))++;//using the index to represent the little cube
}}
// counting the ions number appeared in different zone;


void summary(int *zone,int *zone_summ,int zxindex,int zyindex,int zzindex){
//after counting, the matrix should be summarized as a standard form
//the standard form is n*4 matrix, in which, the n means the total number of the cube.

 for (int x = 0; x < zxindex; x++){
   for (int y = 0; y < zyindex; y++){
    for (int z = 0; z < zzindex; z++){
 // printf("%d %d %d %d\n",x,y,z,*(zone+x*zyindex*zzindex+y*zzindex+z));
  *(zone_summ+(x*zyindex*zzindex+y*zzindex+z)*4+3)=*(zone+x*zyindex*zzindex+y*zzindex+z);
  *(zone_summ+(x*zyindex*zzindex+y*zzindex+z)*4)=x;
  *(zone_summ+(x*zyindex*zzindex+y*zzindex+z)*4+1)=y;
  *(zone_summ+(x*zyindex*zzindex+y*zzindex+z)*4+2)=z;
   }
  } 
 }
}


void sort (int *sort, int row){
//this is the bubble sorting method. the input should be have 4 columns,
//the first 3 columns mean the x, y, z index of the cube.
//the 4th column means the total ions number in the cube.
//the 4th column is the one used to sort. 
  for (int i = 0; i < row; i++){
     double index = *(sort+i*4+3);
     double x = *(sort+i*4);
     double y = *(sort+i*4+1);
     double z = *(sort+i*4+2);
     for (int j = i; j < row ; j++) {
      if( index < *(sort+j*4+3))
      {*(sort+i*4+3) = *(sort+j*4+3);
      *(sort+i*4)= *(sort+j*4);
      *(sort+i*4+1) = *(sort+j*4+1);
      *(sort+i*4+2) = *(sort+j*4+2);

      *(sort+j*4+3) = index;
      *(sort+j*4) = x;
      *(sort+j*4+1) = y;
      *(sort+j*4+2) = z;

       index = *(sort+i*4+3);
       x = *(sort+i*4);
       y = *(sort+i*4+1);
       z = *(sort+i*4+2);

       }
    }
  }
}

void sort_d (double *sort, int row){
//this is the bubble sorting method. the input should be have 4 columns,
//the first 3 columns mean the x, y, z index of the cube.
//the 4th column means the total ions number in the cube.
//the 4th column is the one used to sort. 
  for (int i = 0; i < row; i++){
     double index = *(sort+i*4+3);
     double x = *(sort+i*4);
     double y = *(sort+i*4+1);
     double z = *(sort+i*4+2);
     for (int j = i; j < row ; j++) {
      if( index < *(sort+j*4+3))
      {*(sort+i*4+3) = *(sort+j*4+3);
      *(sort+i*4)= *(sort+j*4);
      *(sort+i*4+1) = *(sort+j*4+1);
      *(sort+i*4+2) = *(sort+j*4+2);

      *(sort+j*4+3) = index;
      *(sort+j*4) = x;
      *(sort+j*4+1) = y;
      *(sort+j*4+2) = z;

       index = *(sort+i*4+3);
       x = *(sort+i*4);
       y = *(sort+i*4+1);
       z = *(sort+i*4+2);

       }
    }
  }
}
void cut_num(int *cut_row, int *zone_summ, double ave, int zxindex, int zyindex, int zzindex){
//the clustering should be stopped by the average of the number of all ions.
//using an index  cut_row to stop the classification.
   *cut_row = 0;
  for (int i = 0; i<zxindex*zyindex*zzindex; i++){
  (*cut_row)++;
  int index = *(zone_summ+4*i+3);
  if(index < ave)
 break;     
}
}



void cluster_ini(int *ion, int* ion_index, int *zone_summ, int sum_cube, int ave){
//initial clustering is to find the top N cube as the initial cluster
//the rule is to make sure these cube do not contact to each other.

//initialize the orinial matrix for further comparison and  cluster;

 *(ion)=*zone_summ;
 *(ion+1)=*(zone_summ+1);
 *(ion+2)=*(zone_summ+2);
 *(ion+3)=*(zone_summ+3);
//first cube donnot need to be compare;
//we got one cluster;
int cou = 0;
 for (int j = 1; j < sum_cube; j++){
  int k = 0;
 for (int jj = 0; jj < sum_cube; jj++){
 int Dist_sqr = pow(*(zone_summ+j*4+0) - *(ion+jj*4+0),2) + pow(*(zone_summ+j*4+1)-*(ion+jj*4+1),2) + pow(*(zone_summ+j*4+2)-*(ion+jj*4+2),2);
 if (Dist_sqr > 3 )
k++;}
  if (k > (sum_cube-1)){
 cou++;
 *(ion+cou*4+0) = *(zone_summ+j*4+0);
 *(ion+cou*4+1) = *(zone_summ+j*4+1);
 *(ion+cou*4+2) = *(zone_summ+j*4+2);
 *(ion+cou*4+3) = *(zone_summ+j*4+3);
  }
 if (*(zone_summ+j*4+3) < ave)
 break;
}

(*ion_index) = (cou+1);
//the total cluster number

}

void cluster_classical(int *ion, double *new_ion, int *zone_summ, int *cut_row, int swen,double Xori,double Yori, double Zori, double cubesize){

 for (int i = 0; i<*cut_row;i++) {
   for (int j = 0; j<swen;j++){
  int dis = pow(*(zone_summ+4*i)-*(ion+j*4),2)+ pow(*(zone_summ+4*i+1)-*(ion+j*4+1),2)+ pow(*(zone_summ+4*i+2)-*(ion+j*4+2),2);
  if (dis < 4)
 {
   *(new_ion+4*j+3) += *(zone_summ+4*i+3);
//   printf("%d %f %d\n",j,*(new_ion+4*j+3),*(zone_summ+4*i+3));
   break;
 }    
}
}
//first loop to calculate the total ions number of the whole cluster.usually 3*3*3;

 for (int ii = 0; ii<*cut_row;ii++) {
   for (int jj = 0; jj<swen; jj++){
  int dis = pow(*(zone_summ+4*ii)-*(ion+jj*4),2) + pow(*(zone_summ+4*ii+1)-*(ion+jj*4+1),2) + pow(*(zone_summ+4*ii+2)-*(ion+jj*4+2),2);
  if (dis < 4)
   {
   double ratio = (*(zone_summ+4*ii+3))/(*(new_ion+4*jj+3));  
//   printf("%f\t",ratio);
   *(new_ion+4*jj) += ((*(zone_summ+4*ii)+0.5)*cubesize+Xori)*ratio;
   *(new_ion+4*jj+1) += ((*(zone_summ+4*ii+1)+0.5)*cubesize+Yori)*ratio;
   *(new_ion+4*jj+2) += ((*(zone_summ+4*ii+2)+0.5)*cubesize+Zori)*ratio;
//   printf("%f\t%f\t%f\n",*(new_ion+4*jj),*(new_ion+4*jj+1),*(new_ion+4*jj+2));
   break;
   } 
  
}
}
}

int main(int argc, char *argv[]) {


FILE *ffra = fopen("frame.number","r");
double wc[3];
for (int i = 0 ; i < 3; i++){
fscanf(ffra,"%lf",wc+i);
}
double FRAMENUM = wc[0];

fclose(ffra);

double cubesize = atof(argv[1]);


unsigned int BUFFER_SIZE = 1024;
   FILE *fpread = fopen("SSS.SSS", "r");
   char line[BUFFER_SIZE];
int row = 0;
int col = 3;
   while (fgets(line, BUFFER_SIZE, fpread) != NULL) {
       row++;
   }
rewind(fpread);
 double* dat = (double*) malloc(3*row*sizeof( double ) );
//the pointer dat is to store all ions postion which is n*3;
 memset(dat,0.0,3*row*sizeof(double));
for (int i = 0; i<row*3; i++){
fscanf(fpread,"%lf",dat+i);}
//input the data from file to prepared matrix dat
fclose (fpread);

double Xori=*(dat),Xfin=*(dat),Yori=*(dat+1),Yfin=*(dat+1),Zori=*(dat+2),Zfin=*(dat+2);

for (int i = 0; i<row; i++){
if (*(dat+i*3) < Xori){
Xori = *(dat+i*3);
}
else if(*(dat+i*3)>Xfin){
Xfin = *(dat+i*3);
}
if (*(dat+i*3+1) < Yori){
Yori = *(dat+i*3+1);
}
else if(*(dat+i*3+1)>Yfin){
Yfin = *(dat+i*3+1);
}
if (*(dat+i*3+2)<Zori){
Zori = *(dat+i*3+2);
}
else if(*(dat+i*3+2)>Zfin){
Zfin = *(dat+i*3+2);
}}
const double X = Xfin - Xori;
const double Y = Yfin - Yori;
const double Z = Zfin - Zori;
int zxindex = (X/cubesize + 1);
int zyindex = (Y/cubesize + 1);
int zzindex = (Z/cubesize + 1);
int row_tot = zxindex*zyindex*zzindex;




double ave;
int sum_cube = zxindex*zyindex*zzindex;
//the number of cubes;
ave = row/(double)sum_cube;
//the average of ions in all cubes;
int *zone = (int*) malloc(zxindex*zyindex*zzindex*sizeof(int));
//pointer "zone" is the solvate box (zxindex * zyindex *zzindex);
memset (zone, 0, zxindex*zyindex*zzindex*sizeof(int));
count(dat,zone,row,zxindex,zyindex,zzindex,Xori,Yori,Zori, cubesize);
//run the function count to calculate the ions number 
//from matrix "dat" to "zone"

int *zone_summ = (int*) malloc(4*zxindex*zyindex*zzindex*sizeof(int));
memset(zone_summ,0,4*zxindex*zyindex*zzindex*sizeof(int));
summary(zone,zone_summ,zxindex,zyindex,zzindex);
//summary the solvate box from 3D box to a 2D matrix;
//data from zone to zone_summ;

sort(zone_summ,sum_cube);
//ranking all numbers; no data transfer;
//bubble sorting;

int *cut_row;
int jack = 0;
cut_row = &jack;
//we need to declare a variable and then the pointer can be used.
//the cut_row (jack) is used for classification

cut_num(cut_row, zone_summ,ave,zxindex,zyindex,zzindex);
//use the threshold (ave) to calculate the how many rows should be 
//taken into classification.

int* ionnumpred = (int*) malloc(4*sum_cube*sizeof(int));
memset (ionnumpred,0,4*sum_cube*sizeof(int));
for(int io = 0; io < 4*sum_cube; io++){ *(ionnumpred+io) = -1000;}
//initialize the room to store the first cube of every cluster;
int swen = 0;
int* ionnumpred_index = &swen;
//set the number to log how many cluster do we have;

cluster_ini(ionnumpred,ionnumpred_index,zone_summ,sum_cube,ave);

double* new_ion = (double*) malloc(4*swen*sizeof(double));
memset(new_ion,0.0,4*swen*sizeof(double));
cluster_classical(ionnumpred,new_ion,zone_summ,cut_row,swen,Xori,Yori,Zori,cubesize);

sort_d(new_ion,swen);
//for (int sss = 0; sss < 10; sss++){
//printf("%s %d %d %d %d\n","the sss cluster is x, y, z, counting",*(ionnumpred+sss*4), *(ionnumpred+sss*4+1), *(ionnumpred+sss*4+2), *(ionnumpred+sss*4+3));
//}
for (int HH = 0; HH < swen; HH++){
printf("%f %f %f %f %f\n",*(new_ion+HH*4), *(new_ion+HH*4+1), *(new_ion+HH*4+2), *(new_ion+4*HH+3),*(new_ion+HH*4+3)/(double)FRAMENUM);
//printf("%d %d %d %d\n", ion[HH][0], ion[HH][1], ion[HH][2],ion[HH][3]);
}
printf("%s\t%d\n","the number of cubes is",sum_cube);
printf("%s\t%d\n","the number of total ions is",row);
printf("%s\t%f\n","the average of ions in cube is",ave);
printf("%s\t%d\n","the stop row is",jack);
printf("%s\t%d\n","the number of clusters is",swen);

free(dat);
free(zone);
free(zone_summ);
free(ionnumpred);
free(new_ion);
return 0;
}


