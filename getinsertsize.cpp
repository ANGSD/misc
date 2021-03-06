#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <cassert>
#include <cmath>
#include <map>
#define LENS 100000


typedef std::map<int,size_t> aMap;

size_t *rlen = new size_t[LENS];
size_t *span = new size_t[LENS];

aMap rlen_map;
aMap span_map;


void myIns(int i,aMap &asso){
  aMap::iterator it=asso.find(i);
  if(it!=asso.end())
    it->second++;
  else
    asso.insert(std::pair<int,size_t>(i,1));
}

double mymean(size_t *ary,aMap & asso,int maxbound){
  double nsum = 0;
  double nobs =0;
  for(size_t i=0;i<LENS;i++){
    if((maxbound!=-1)&& (i>maxbound))
      continue;
    if(ary[i]){
      nsum += i*ary[i];
      nobs += ary[i];
    }
    
  }

  for(aMap::iterator it=asso.begin();it!=asso.end();++it){
    if(maxbound!=-1&&it->first>maxbound)
      continue;
    
    nsum += (it->first)*(it->second);
    nobs += it->second;
  }
  return nsum/nobs;
}

int whichmax(size_t *ary,aMap &asso){
  int id =0;
  for(int i=1;i<LENS;i++)
    if(ary[i]>ary[id])
      id=i;
  for(aMap::iterator it=asso.begin();it!=asso.end();++it)
    if(it->second>ary[id])
      id=it->first;
  
  return id;
}

double myvar(size_t *ary,aMap &asso,double mean,int maxbound){
  // fprintf(stderr,"maxbound:%d\n",maxbound);
  double nsum = 0;
  double nobs =0;
  for(size_t i=0;i<LENS;i++){
    if((maxbound!=-1)&&(i>maxbound))
      continue;
    if(ary[i]){
      nsum += (i-mean)*(i-mean) *ary[i];
      nobs += ary[i];
    }
    
  }

  for(aMap::iterator it=asso.begin();it!=asso.end();++it){
    if(maxbound!=-1&&it->first>maxbound)
      continue;
    nsum += (it->first-mean)*(it->first-mean)*(it->second);
    nobs += it->second;
  }
  return sqrt(nsum/(nobs-1));
}

void myprint(size_t *ary,aMap &asso,FILE *fp){
  for(int i=0;i<LENS;i++)
    if(ary[i])
      fprintf(fp,"%d\t%zu\n",i,ary[i]);
  for(aMap::iterator it=asso.begin();it!=asso.end();++it)
    fprintf(fp,"%d\t%zu\n",it->first,it->second);
	
  

}

size_t stat[5]={0,0,0,0,0};

#define unmapped 0
#define mapped 1
#define orphan 2
#define paired 3
#define used 4

int main() {
  if(isatty(fileno(stdin))){
    fprintf(stderr,"pipe uncompressed samfile (no header) into program\n");
    return 0;
  }

  memset(rlen,0,sizeof(size_t)*LENS);
  memset(span,0,sizeof(size_t)*LENS);
  char *buf = new char[LENS];
  char **toks = new char*[12];
  size_t nreads=0;
  while(fgets(buf,LENS,stdin)){
    nreads++;
    //    fprintf(stderr,"buf:%s\n",buf);
    toks[0] = strtok(buf,"\n\t ");
    for(int i=1;i<12;i++){
      toks[i] = strtok(NULL,"\n\t ");
    }
    
    int flag = atoi(toks[1]);
    if(flag &4){
      //      fprintf(stderr,"Read unmapped \'%s\'\n",toks[1]);
      stat[unmapped]++;
      continue;
    }else
      stat[mapped]++;


    if(!(flag&1)){//check if odd
      fprintf(stderr,"Alignment record must be paired \'%s\'\n",toks[1]);
      return 0;
    }

    #if 0
    //disregard flag column, which is only meanning full if isize has been specified in bwa mapping
    if(!(flag & 2)){
      stat[orphan]++;
      continue;
    }else
      stat[paired]++;
    #endif
    
    if(flag&8){
      stat[orphan]++;
      continue;
    }else
      stat[paired]++;

    

    char *cig = toks[5];
    const char *tmp = "IDNSHP=X";
    int cont=1;
    for(unsigned i=0;cont&&i<strlen(tmp);i++)
      if(strchr(cig,tmp[i]))
	cont =0;
    if(cont==0){
      //   fprintf(stderr,"Unperfect match: \'%s\'\n",cig);
      continue;
    }
    if(toks[6][0]!='=')
      continue;
    stat[used]++;

    //    fprintf(stderr,"Perfect Match:%s\n",cig);
    char *tok = strchr(cig,'M');
    assert(tok);
    *tok = '\0';
    int len = atoi(cig);
    if(len<LENS-10)
      rlen[len]++;
    else
      myIns(len,rlen_map);

    if(toks[6][0]!='=')
      continue;
    //    fprintf(stderr,"tok6:%s\n",toks[6]);

    if(((len = atoi(toks[8])))<=0)
      continue;
    
    // fprintf(stderr,"spanlen:%d\n",len);
    if(len<LENS-10)
      span[len]++;
    else
      myIns(len,span_map);
    
  }
  FILE *fp = stdout;
  fprintf(fp,"\"nReads(total number):\"%zu\n",nreads);
  fprintf(fp,"\"unmapped(flag not 4):\"%zu\n",stat[0]);
  fprintf(fp,"\"mapped(flag 4):\"%zu\n",stat[1]);
  fprintf(fp,"\"orphan(flag 8):\"%zu\n",stat[2]);
  fprintf(fp,"\"paired(flag not 8):\"%zu\n",stat[3]);
  fprintf(fp,"\"used(same chr and mate mapped):\"%zu\n",stat[4]);
#if 1
  fprintf(stdout,"#insertsizes len count");
  myprint(span,span_map,stdout);
  fprintf(stdout,"#rlen len count");
  myprint(rlen,rlen_map,stdout);
#endif
  double mea = mymean(rlen,rlen_map,-1);
  double var = myvar(rlen,rlen_map,mea,-1);
  //  fprintf(stdout,"Read length mean:%f\t sqrt(var):%f\n",mea,var);

  int pivot = whichmax(span,span_map);
  mea = mymean(span,span_map,3*pivot);
  var = myvar(span,span_map,mea,3*pivot);
  //fprintf(stdout,"Read span mean:%f\t sqrt(var):%f\n",mea,var);
  
  //cleanup
  delete [] rlen;
  delete [] span;
  delete [] buf;
  delete [] toks;
  return 0;
}
