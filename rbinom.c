/*
  BINV from 'binomial random variate generation, voratas kachitvichyanukul, bruce w. schmeiser, feb 1988 vol 31, nr 2
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
int binv(int n,double p){
  double q=1-p;
  double s=p/q;
  double a=(n+1)*s;
  long double r=pow(q,n);
  double u=drand48();
  int x=0;

  while(u>r){
    u-=r;x++;
    r=((a/(double)x)-s)*r;
  }
  return x;
}

#ifdef __WITH_MAIN__
int main(){
  for(int i=0;i<1e6;i++){
    fprintf(stdout,"%d\n",binv(100,0.6));


  }

}
#endif
