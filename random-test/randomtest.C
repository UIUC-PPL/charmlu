#include <iostream>
#include <cassert>

#include <sys/time.h>

using namespace std;

int random_seed;


/// Return a random integer between 0 and num-1 inclusive
unsigned int randInt(unsigned int num, const char* name, int seed=0){
  assert(num < 1000);

  unsigned long hash = 0;
  unsigned int c;
  unsigned char * str = (unsigned char*)name;

  while (c = *str++){
    unsigned int c2 = (c+64)%128;
    unsigned int c3 = (c2*5953)%127;
    hash = c3 + (hash << 6) + (hash << 16) - hash;
  }
  
  unsigned long shuffled1 = (hash*2083)%7907;
  unsigned long shuffled2 = (seed*4297)%2017;

  unsigned long shuffled3 = (random_seed*4799)%7919;

  unsigned int namehash = shuffled3 ^ shuffled1 ^ shuffled2;

  unsigned int result = ((namehash * 6029) % 1117) % num;


  //  cout << " result = " << result << " seed = " << seed << " random_seed = " << random_seed  << endl;

  assert(result >=0 && result < num);
  return result;
}




int main(){
  int range = 2;

  struct timeval tp;
  gettimeofday(& tp, NULL);
  random_seed = (((int)tp.tv_usec*7537)%4751) ^ (((int)tp.tv_sec * 2719)%6343);

  int match = 0;
  int nomatch = 0;

  for(int i=0;i<10000;i++){

    if(randInt(4, "block size", i) == randInt(2,"which mapping",i) )
      match++;
    else
      nomatch++;
  }

  cout << "match: " << match << endl;
  cout << "nomatch: " << nomatch << endl;

  cout << "match/nomatch = " << (double)match/(double)nomatch << " expected "  << 1.0 / 4.0  <<  endl;
  

  return 0;
}
