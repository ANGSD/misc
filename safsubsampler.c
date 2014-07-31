#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <sys/stat.h>

char tmp[1024];

size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}


int main(int argc,char**argv){
  long seed = time(NULL);
  fprintf(stderr,"seed:%ld\n",seed);
  if(argc!=4){
    fprintf(stderr,"Program to subsample saf files, supply ./safsubsamler file.saf nchr outname prefix\n");
    return 0;
  }
  FILE *fp = NULL;
  if(((fp=fopen(argv[1],"rb")))==NULL){
    fprintf(stderr,"Problem opening file: \'%s\'\n",argv[1]);
    return 0;
  }
  int nChr=0;
  if(((nChr=atoi(argv[2])))<=0){
    fprintf(stderr,"need to supply a positive value for the second parameter (the number of chrs)\n");
    return 0;
  }
  char *prefix = strdup(argv[3]);
  

  snprintf(tmp,1024,"%s.%ld.saf",prefix,seed);
  FILE *fpout= NULL;
  if(((fpout=fopen(tmp,"w")))==NULL){
    fprintf(stderr,"Problem opening file: \'%s\'\n",tmp);
    return 0;
  }
  fprintf(stderr,"ARGS: fname:%s nchr:%d out.saf:%s out.bin:",argv[1],nChr,tmp);
  snprintf(tmp,1024,"%s.%d.bin.gz",prefix,seed);
  gzFile gz= Z_NULL;
  if(((gz=gzopen(tmp,"wb")))==Z_NULL){
    fprintf(stderr,"Problem opening file: \'%s\'\n",tmp);
    return 0;
  }
  fprintf(stderr,"%s\n",tmp);
  
  if(fsize(argv[1]) % (nChr+1)*sizeof(double) ){
    fprintf(stderr,"Problem with saffile looks truncated\n");
    return 0;
  }


  //Program can either sample based on previously calculated offsets or recompute offsets
  int nSites = fsize(argv[1])/(nChr+1)/sizeof(double);
  fprintf(stderr,"file contains: %d sites\n",nSites);

  int *offs = malloc(nSites*sizeof(int));
  

  if(!isatty(fileno(stdin)) ){
    fprintf(stderr,"will load offsets\n");
    
    if(nSites != fread(offs,sizeof(int),nSites,stdin)){
      fprintf(stderr,"problening reading full chunk from uncompressed piped data .bin file\n");
      return 0;
    }
  }else{
    fprintf(stderr,"will compute offsets\n");
    for(int i=0;i<nSites;i++){
      offs[i] = lrand48() %nSites;
      //    fprintf(stdout,"%d\n",offs[i]);
    }
  }
  //  return 0;
  if(sizeof(int)*nSites!=gzwrite(gz,offs,sizeof(int)*nSites)){
    fprintf(stderr,"Problem writing offsets into compressed file\n");
    return 0;
      
  }
  
  double *d=malloc(fsize(argv[1]));
  fread(d,nSites*(nChr+1),sizeof(double),fp);
  for(int i=0;i<nSites;i++){
    //    fprintf(stderr,"i:%d offs:%d offB:%lu\n",i,offs[i],offs[i]*sizeof(double)*(nChr+1));
    if(nChr+1!=fwrite(d+offs[i]*(nChr+1),sizeof(double),nChr+1,fpout)){
      fprintf(stderr,"Problem writing data to file\n");
      return 0;
    }
  }
  gzclose(gz);fclose(fpout);fclose(fp);
  free(d);free(offs);
  return 0;
}
