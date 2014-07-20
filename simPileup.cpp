/*
  simulate pileupdata, 
  1) simulating genotypes
  2) simulate reads and qscores according to the genotypes

 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <cmath>
#include <cassert>

float bfreq[4]={0.25,0.5,0.75,1};
char bases[5] = "ACGT";
double minfreq =1e-4, myConst;

template<typename T>
void norm(T*ary,size_t x){
  assert(x>0);
  T ssum=ary[0];
  for(int i=1;i<x;i++)     ssum +=ary[i];
  fprintf(stderr,"sum:%f\n",ssum);
  for(int i=0;i<x;i++) ary[i] /=ssum;
}




double simfreq() {
  return exp(drand48()*myConst-myConst);
}

int rpois(int mean){
  double limit = exp(-mean);
  double x = drand48();
  int ret =0;
  while (x > limit) {
    ret++;
    x *= drand48();
  }
  return ret;
}

int pickInterval(int n,float points[4]){
  double r=drand48();
  // fprintf(stderr,"r:%f\n",r);
  if(r<points[1]){
    // fprintf(stderr,"a,c\n");
    if(r<points[0])
      return 0;
    else
      return 1;
  }else{
    //    fprintf(stderr,"g,t\n");
    if(r<points[2])
      return 2;
    else
      return 3;
  }

}


void per_ind_hap(char genos,int dep,float err){
  //  fprintf(stderr,"[%s] genos:%d dep:%d err:%f\n",__FUNCTION__,genos,dep,err);
  for(int d=0;d<dep;d++){
  
    float r = drand48();
    //fprintf(stderr,"d:%d r:%f err:%f\n",d,r,err);  
    if(r>err)
      fprintf(stdout,"%c" , bases[genos]);
    else{
      char b;
      while(1){
	int i = pickInterval(4,bfreq);
	fprintf(stderr,"i:%d\n",i);
	if(i!=genos){
	  b=i;
	  break;
	}
      }
      fprintf(stdout,"%c" , bases[b]);
    }
  }
  fprintf(stdout,"\t");

}


int main(int argc,char **argv){

  myConst = -log(minfreq);
  
  int mean =4;
  int ploidy =1;
  int nind =5;
  int nsites=1e5;
  char *outname=NULL;
  int n;
  float var=0.6;
  float err = 0.00;
  if (argc==1) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: simPileup  [options] \n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -p <int>            list of positions or regions\n");
    fprintf(stderr, "   -i <int>            list of input BAM filenames, one per line [null]\n");
    fprintf(stderr, "   -s <int>            minQLen\n");
    fprintf(stderr, "   -o <filename>       base quality threshold\n");
    fprintf(stderr, "   -v <float>       variability\n");
    fprintf(stderr, "   -m <int>       mean seq depth\n");
    fprintf(stderr, "   -e <int>      errorrateh\n");
    fprintf(stderr, "\n");
    return 1;
  }

  double *truefreq=NULL;


  while ((n = getopt(argc, argv, "p:i:s:o:v:m:e:")) >= 0) {
    switch (n) {
    case 'p':ploidy = atoi(optarg); break; 
    case 'i': nind = atoi(optarg); break;  
    case 's': nsites = atoi(optarg); break;
    case 'o': outname = strdup(optarg); break;
    case 'm': mean = atoi(optarg); break;
    case 'v': var = atof(optarg); break;
    case 'e': err = atof(optarg); break;
    }
  }
  fprintf(stderr,"p:%d i:%d s:%d o:%s m:%d v:%f\n",ploidy,nind,nsites,outname,mean,var);
  size_t dim=ploidy==1?nind+1:2*nind+1;
  fprintf(stderr,"dim:%zu\n",dim);
  truefreq=ploidy==1?new double[dim]:new double[dim];
  for(int i=0;i<dim;i++) truefreq[i] =0;
  
  assert(ploidy==1);
  
  for(int s=1;s<=nsites;s++){
    fprintf(stdout,"chr1\t%d\t",s);
    char ref=pickInterval(4,bfreq);//<- 0,1,2,3
    char alt=ref;
    if(drand48()<=var){
      char b = pickInterval(4,bfreq);
      //  fprintf(stderr,"b:%d\n",b);
      alt = b==ref?(b+1)%4:b;//<- 0,1,2,3
    }
    fprintf(stdout,"%c\t%c\t",bases[ref],bases[alt]);
   
    float sfreq = simfreq();
    fprintf(stdout,"%f\t",sfreq);
    int aCount =0;
    for(int n=0;n<nind;n++){
      char b=drand48()>sfreq?ref:alt;
      if(b!=ref) aCount++;
      int d = rpois(mean);
      //fprintf(stdout,"%c\t",bases[b]);
      //exit(0);
      //fprintf(stderr,"-----b:%d\n",b);
      per_ind_hap(b,d,err);
      //fprintf(stderr,"-----\n");
      //exit(0);
    }
    truefreq[aCount]++;
  }
  fprintf(stdout,"\n");
  norm(truefreq,dim);
  fprintf(stderr,"truefreq:\n");
  for(int i=0;i<dim;i++)
    fprintf(stderr,"%f\t",truefreq[i]);
  fprintf(stderr,"\n");

  return 0;
}
