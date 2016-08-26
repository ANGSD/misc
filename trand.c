#include <stdlib.h>

typedef struct {
  struct drand48_data drand_buf;
  double res;
}myStruct;

#define set_seed(tr,seed) srand48_r(seed,&tr.drand_buf)
#define uni(tr) drand48_r(&tr.drand_buf,&tr.res)?tr.res:tr.res

int main(){
  myStruct ms;
  set_seed(ms,10);
  uni(ms);

}
