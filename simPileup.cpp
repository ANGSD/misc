

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
#include <vector>

const char * PLP =".pileup";
const char * GENO =".geno";
const char * TRUEFREQ =".truefreq";

std::vector<char *> ofiles;

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

FILE *getFP(char *fname,const char *post){
  fprintf(stderr,"fnaem:%s post:%s\n",fname,post);
  unsigned len = strlen(fname)+strlen(post)+1;
  char *tmp =(char*) malloc(strlen(fname)+strlen(post)+1);
  snprintf(tmp,len,"%s%s",fname,post);
  
  FILE *fp = NULL;
  if(!((fp=fopen(tmp,"wb")))){
    fprintf(stderr,"Problem opening file: \'%s\'\n",tmp);
    exit(0);
  }
  ofiles.push_back(tmp);
  return fp;
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


void per_ind_hap(char genos,int dep,float err,FILE *fp){
  //  fprintf(stderr,"[%s] genos:%d dep:%d err:%f\n",__FUNCTION__,genos,dep,err);
  for(int d=0;d<dep;d++){
  
    float r = drand48();
    //fprintf(stderr,"d:%d r:%f err:%f\n",d,r,err);  
    if(r>err)
      fprintf(fp,"%c" , bases[genos]);
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
      fprintf(fp,"%c" , bases[b]);
    }
  }
  fprintf(fp,"\t");

}

void print_info(FILE *fp){
   fprintf(stderr, "\n");
    fprintf(stderr, "Usage: simPileup  [options] \n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -p <int>            ploidy\n");
    fprintf(stderr, "   -i <int>            number of individuals\n");
    fprintf(stderr, "   -s <int>            number of sites\n");
    fprintf(stderr, "   -o <filename>       outputfilename\n");
    fprintf(stderr, "   -v <float>          variability\n");
    fprintf(stderr, "   -m <int>            mean seq depth\n");
    fprintf(stderr, "   -e <int>            errorrateh\n");
    fprintf(stderr, "\n");
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
    print_info(stderr);
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
  if(!outname){
    fprintf(stderr,"Must supply -o\n");
    print_info(stderr);
    return 1;
  }
  FILE *fp_pileup = getFP(outname,PLP);
  FILE *fp_geno = getFP(outname,GENO);
  fprintf(stderr,"p:%d i:%d s:%d o:%s m:%d v:%f\n",ploidy,nind,nsites,outname,mean,var);
  size_t dim=ploidy==1?nind+1:2*nind+1;
  fprintf(stderr,"dim:%zu\n",dim);
  truefreq=ploidy==1?new double[dim]:new double[dim];
  for(int i=0;i<dim;i++) truefreq[i] =0;
  
  assert(ploidy==1);
  
  for(int s=1;s<=nsites;s++){//over sites
    fprintf(fp_pileup,"chr1\t%d\t",s);
    char ref=pickInterval(4,bfreq);//<- 0,1,2,3
    char alt=ref;

    if(drand48()<=var){
      char b = pickInterval(4,bfreq);
      alt = b==ref?(b+1)%4:b;//<- 0,1,2,3
    }

    fprintf(fp_pileup,"%c\t%c\t",bases[ref],bases[alt]);
   
    float sfreq = simfreq();
    // fprintf(stdout,"%f\t",sfreq);
    int aCount =0;
    for(int n=0;n<nind;n++){
      char b=drand48()>sfreq?ref:alt;
      fprintf(fp_geno,"%c ",bases[b]);
      if(b!=ref) aCount++;
      int d = rpois(mean);
      per_ind_hap(b,d,err,fp_pileup);
    }
    
    fprintf(fp_pileup,"\n");
    fprintf(fp_geno,"\n");
    truefreq[aCount]++;
  }
  fprintf(stdout,"\n");
  norm(truefreq,dim);
  fprintf(stderr,"truefreq:\n");
  for(int i=0;i<dim;i++)
    fprintf(stderr,"%f\t",truefreq[i]);
  fprintf(stderr,"\n");
  #if 1
  fprintf(stderr,"Dumping files\n");
  for(unsigned i=0;i<ofiles.size();i++)
    fprintf(stderr,"\t%s\n",ofiles[i]);

  #endif
  return 0;
}
