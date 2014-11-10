#include "faidx.h"
#include <pthread.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <ctype.h>
typedef struct{
  char *fastaname;
  pthread_mutex_t aMut;//mutex will be locked unlocked whenever we do something
  faidx_t *fai;//contains the faidx structure
  char *seqs;//contains the reference for the current chr;
  int curChr;//the exact chromosome name for the seqs above
  int chrLen;//length of chromosome
}perFasta;



//this will initialize our data
perFasta *init( char *fname){

  fprintf(stderr,"\t-> Reading fasta: %s\n",fname);

  //check that fa hasn't been updated
  char *strtsk=NULL;
  strtsk = (char*)calloc(strlen(fname) + 5, 1);
  sprintf(strtsk, "%s.fai", fname);
  free(strtsk);
    
  perFasta *r= new perFasta;
  r->fastaname=NULL;
  r->fai = NULL;
  r->seqs = NULL;
  r->curChr = -1;
  if(pthread_mutex_init(&r->aMut,NULL)){
    fprintf(stderr,"[%s:%s] error initializing mutex \n",__FILE__,__FUNCTION__);
    exit(0);
  }

  if(NULL==(r->fai = fai_load(fname))){
    fprintf(stderr,"[%s:%s] error reading fai file:%s\n",__FILE__,__FUNCTION__,fname);
    exit(0);
  }
  r->fastaname=strdup(fname);
  return r;
}


char *loadChr(perFasta *f, char*chrName,int chrId){
  fprintf(stderr,"[%s] \t->loading chr:%s from faidx=%p curchr=%d\n",__FUNCTION__,chrName,f,f->curChr);
  free(f->seqs);
  f->seqs=NULL;
  //fprintf(stderr,"f->curChr=%d chrId=%d\n",f->curChr,chrId);
  f->curChr = chrId;
  f->seqs = faidx_fetch_seq(f->fai, chrName, 0, 0x7fffffff, &f->chrLen);
  if(f->seqs==NULL){
    fprintf(stderr,"\n[%s] Error loading fasta info from chr:\'%s\' \n",__FUNCTION__,chrName);
    f->chrLen=0;
  }

  return f->seqs;
}




int main(int argc,char **argv){
  perFasta **fs = new perFasta*[argc-1];
  for(int i=1;i<argc;i++){
    fs[i-1] = init(argv[i]);
  }
  int nit=0;
  char **chr=getnam(fs[0]->fai,&nit);

  fprintf(stdout,"chr\tpos");
  for(int i=1;i<argc;i++)
    fprintf(stdout,"\t%s",basename(argv[i]));
  fprintf(stdout,"\n");

  for(int i=0;i<nit;i++){
    for(int j=0;j<argc-1;j++)
      loadChr(fs[j],chr[i],i);
    for(int j=1;j<argc-1;j++)
      assert(fs[0]->chrLen==fs[j]->chrLen);
    for(int j=0;j<fs[0]->chrLen;j++){
      char b=toupper(fs[0]->seqs[j]);
      if(b=='N')
	continue;
      int skip=1;
      for(int f=1;f<argc-1;f++)
	if(b!=toupper(fs[f]->seqs[j])){
	  skip=0;
	  break;
	}
      if(skip==0){
	fprintf(stdout,"%s\t%d",chr[i],j);
	for(int f=0;f<argc-1;f++)
	  fprintf(stdout,"\t%c",toupper(fs[f]->seqs[j]));
	fprintf(stdout,"\n");
      }
    }
  }
    
  
  return 0;
} 
