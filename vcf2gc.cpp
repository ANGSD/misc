#include <cstdio>
#include <cstring>

#define LENS 40960


int print(char *tmp){
  if(tmp[0]=='.'&&tmp[2]=='.')
    return -1;
  else {
    int en = tmp[0]-'0';
    int to = tmp[2]-'0';
    //    fprintf(stderr,"%d %d\n",en,to);
    return en+to;
  }

}


int main(){
  char *buf = new char[LENS];
  
  while(fgets(buf,LENS,stdin)){
    if(buf[0]=='#')
      continue;
    //fprintf(stderr,"%s",buf);
    fprintf(stdout,"%s\t",strtok(buf,"\n\t "));
    fprintf(stdout,"%s\t",strtok(NULL,"\n\t "));
    strtok(NULL,"\n\t ");
    fprintf(stdout,"%s\t",strtok(NULL,"\n\t "));
    fprintf(stdout,"%s",strtok(NULL,"\n\t "));
    strtok(NULL,"\n\t ");
    strtok(NULL,"\n\t ");
    strtok(NULL,"\n\t ");
    strtok(NULL,"\n\t ");
    char *tmp;
    while(((tmp=strtok(NULL,"\n\t ")))) {
      //      fprintf(stderr,"%s",tmp);
      int v = print(tmp);
      if(v==-1)
	fprintf(stdout,"\tNA");
      else
	fprintf(stdout,"\t%d",v);
    }
    fprintf(stdout,"\n");
  }

  return 0;
}
