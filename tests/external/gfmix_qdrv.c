#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "gfmixff.h"
/* mode "P": read 20 freqs, 2 gammas, t -> print reference P(t) (column-major)
   mode "rates": read alpha and K     -> print K category rates and weights  */
int main(int argc, char **argv){
  int nchar=20, nca=2, i, j;
  int cia[20]={0,0,1,2,2,2,2,0,2,1,2,1,1,1,0,2,2,2,1,2};
  double fr[20], ca[3], S[400], W[400], Lam[400], M[400], mu, t;
  char mode = (argc>1 && argv[1][0]=='r') ? 'r' : 'P';
  if(mode=='r'){
    double alpha, wt[64], rated[64]; int K;
    if(scanf("%lf %d",&alpha,&K)!=2) return 2;
    DiscreteGamma(wt,rated,alpha,alpha,K,0);
    for(i=0;i<K;i++) printf("%.17g %.17g\n", rated[i], wt[i]);
    return 0;
  }
  for(i=0;i<nchar;i++) if(scanf("%lf",&fr[i])!=1) return 2;
  for(i=0;i<nca;i++)   if(scanf("%lf",&ca[i])!=1) return 2;
  ca[nca]=1.0;
  if(scanf("%lf",&t)!=1) return 2;
  LG_exchangeability(S);
  mu=frca(nchar,ca,cia,fr);
  Sfr2WLam(nchar,S,fr,W,Lam);
  Mij_nod(nchar,t,fr,W,Lam,M);
  for(i=0;i<nchar;i++){
    for(j=0;j<nchar;j++) printf("%.17g ", M[i+j*nchar]);
    printf("\n");
  }
  return 0;
}
