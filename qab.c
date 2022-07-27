//(a,b) optimization testing.


//INCLUSIONS

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef METHDECL
#define METHDECL
#include "qabdecl.h"
#endif

//STATIC DECLARATIONS
static GEN alg_count_Q_worker(GEN founddef, GEN foundindef, long N1, long N2, long prec);

//FUNCTIONS

//BASIC METHODS

//Given (a, b), returns the discriminant of (a,b/Q). a and b can be factored, and we can also pass in [a, b] for a instead.
GEN ab_disc(GEN a, GEN b){
  pari_sp top=avma;
  if(typ(a)==t_VEC){b=gel(a, 2);a=gel(a, 1);}
  GEN afact, bfact;
  if(typ(a)==t_INT) afact=Z_factor(a);
  else{afact=a;a=factorback(afact);}
  if(typ(b)==t_INT) bfact=Z_factor(b);
  else{bfact=b;b=factorback(bfact);}
  GEN plist=ZV_sort_uniq(shallowconcat(gel(afact, 1), gel(bfact, 1)));//May be missing 2, and may include -1
  GEN disc=gen_1;
  if(hilbert(a, b, gen_2)==-1) disc=gen_2;
  long ind=1;
  while(lg(plist)>ind && cmpis(gel(plist, ind), 2)<=0) ind++;//Skipping -1 and 2
  while(ind<lg(plist)){
	GEN p=gel(plist, ind);
	if(hilbert(a, b, p)==-1) disc=mulii(disc, p);
    ind++;	
  }
  return gerepilecopy(top, disc);
}

//Given (a, b), returns the set of ramifying primes in (a,b/Q). a and b can be factored, and we can also pass in [a, b] for a instead.
GEN ab_ramprimes(GEN a, GEN b){
  pari_sp top=avma;
  if(typ(a)==t_VEC){b=gel(a, 2);a=gel(a, 1);}
  GEN afact, bfact;
  if(typ(a)==t_INT) afact=Z_factor(a);
  else{afact=a;a=factorback(afact);}
  if(typ(b)==t_INT) bfact=Z_factor(b);
  else{bfact=b;b=factorback(bfact);}
  GEN plist=ZV_sort_uniq(shallowconcat(gel(afact, 1), gel(bfact, 1)));//May be missing 2, and may include -1
  GEN rprimes=vectrunc_init(lg(plist)+1);
  if(hilbert(a, b, gen_2)==-1) vectrunc_append(rprimes, gen_2);
  long ind=1;
  while(lg(plist)>ind && cmpis(gel(plist, ind), 2)<=0) ind++;//Skipping -1 and 2
  while(ind<lg(plist)){
	GEN p=gel(plist, ind);
	if(hilbert(a, b, p)==-1) vectrunc_append(rprimes, p);
    ind++;	
  }
  return gerepilecopy(top, rprimes);
}

//Returns (a, b) such that A=(a,b/F) has |N_{F/Q}(discriminant)|=D and infinite ramification prescribed by infram (a length [F:Q] vector of 0's/1's), if it exists. If it does not, this returns 0. If F!=Q, then there may be multiple choices for the primes ramifying. This picks the first one possible (when we factorize p).
GEN algabfromram(GEN F, GEN D, GEN infram){
  pari_sp top=avma;
  if(typ(D)!=t_INT || signe(D)!=1) pari_err_TYPE("D should be a positive integer", D);//The absolute norm to Q should be a positive integer
  GEN pfac=Z_factor(D);
  long nfacsp1=lg(gel(pfac, 1));//# prime factors+1
  long nramplaces=nfacsp1-1;
  for(long i=1;i<lg(infram);i++) if(!gequal0(gel(infram, i))) nramplaces++;//Adding the number of oo ramified places
  if(nramplaces%2==1) return gc_const(top, gen_0);//Odd number of ramification places, BAD
  GEN pfacideals=zerovec(nfacsp1-1), hass=cgetg(nfacsp1, t_VEC), possideals;
  long expon;
  for(long i=1;i<nfacsp1;i++){
    expon=itos(gcoeff(pfac, i, 2));//Coefficient of p we desire.
    possideals=idealprimedec(F, gcoeff(pfac, i, 1));
    for(long j=1;j<lg(possideals);j++){
      if(pr_get_f(gel(possideals, j))==expon){
        gel(pfacideals, i)=gel(possideals, j);//We win!
        break;
      }
    }
    if(gequal0(gel(pfacideals, i))) return gc_const(top, gen_0);//Nope, return 0
    gel(hass, i)=gen_1;//The vector of 1's
  }
  GEN A=alginit(F, mkvec3(gen_2, mkvec2(pfacideals, hass), infram), -1, 1);
  return gerepileupto(top, alggetab(A));
}

//Returns [a, b] given an algebra.
GEN alggetab(GEN A){
  pari_sp top=avma;
  GEN L=alg_get_splittingfield(A), pol=rnf_get_pol(L);//L=K(sqrt(a))
  long Lvar=rnf_get_varn(L);
  GEN a=gneg(gsubst(pol, Lvar, gen_0));//Polynomial is of the form x^2-a, so plug in 0 and negate to get a
  GEN b=lift(alg_get_b(A));//This extracts a and b from the algebra.
  return gerepilecopy(top, mkvec2(a, b));
}

//Returns the norm to Q of the discriminant of A
GEN algnormdisc(GEN A){
  pari_sp top=avma;
  GEN nf=alg_get_center(A);//Field
  GEN rams=algramifiedplacesf(A);
  GEN algdisc=gen_1;
  for(long i=1;i<lg(rams);i++) algdisc=mulii(algdisc, idealnorm(nf, gel(rams, i)));//Norm to Q of the ramification
  return gerepileupto(top, algdisc);
}

//Returns the vector of finite ramified places of the algebra A.
GEN algramifiedplacesf(GEN A){
  pari_sp top=avma;
  GEN hass=alg_get_hasse_f(A);//Shallow
  long nhass=lg(gel(hass, 2));
  GEN rp=vectrunc_init(nhass);
  for(long i=1;i<nhass;i++){
    if(gel(hass, 2)[i]==0) continue;//Unramified
    vectrunc_append(rp, gmael(hass, 1, i));//Ramified
  }
  return gerepilecopy(top, rp);
}


//COUNTING ALGEBRAS


//Counts the number of new and distinct algebras (a,b/Q) with N1<=max(|a|,|b|)<=N2. Returns [founddef, foundindef, defcount, indefcount], where each entry is a length N2-N1+1 vecsmall, the ith entry representing the number of algebras with max(|a|,|b|)=N1+i-1. This is not stack clean, and we update the reference to olddiscs, so make sure that the passed olddiscs can be changed.
static GEN alg_count_Q_worker(GEN founddef, GEN foundindef, long N1, long N2, long prec){
  GEN Q=nfinit(pol_x(1), prec);//rational numbers
  pari_sp av=avma;//found(in)def tracks the algebras found
  long ldef=lg(founddef), lindef=lg(foundindef), lN=N2-N1+1;//Tracks the number of algebras found
  GEN countdef=vecsmall_ei(lN, 1), countindef=vecsmall_ei(lN, 1);//Number of new algebras for each N1<=n<=N2
  if(N1>1) countindef[1]=0;//If N1=1, we start with M(2, Q), else we start with no algebras.
  countdef[1]=0;//a=1 only corresponds to M_2(Q), so we start this count at 1.
  long aind=0;
  for(long a=N1;a<=N2;a++){//we do -a first, then a
    aind++;
    if(a%20==0) pari_printf("a=%d done\n", a);
    if(gc_needed(av, 1)) gerepileall(av, 4, &founddef, &foundindef, &countdef, &countindef);
    GEN aval=stoi(-a);
    for(long b=1;b<=a;b++){//We do -b first, then b
      GEN bval=stoi(-b);
      GEN disc=algnormdisc(alginit(Q, mkvec2(aval, bval), 0, 0));
      founddef=setunion_i(founddef, mkvec(disc));
      if(lg(founddef)>ldef){ldef++;countdef[aind]++;}
      if(uissquare(b)) continue;//b is square, no need to do it.
      bval=stoi(b);
      disc=algnormdisc(alginit(Q, mkvec2(aval, bval), 0, 0));
      foundindef=setunion_i(foundindef, mkvec(disc));
      if(lg(foundindef)>lindef){lindef++;countindef[aind]++;}
    }
    if(uissquare(a)) continue;//a is square, no need to do it.
    aval=stoi(a);
    for(long b=1;b<a;b++){//We do -b first, then b. We do a=b at the end, since we don't need to redo the cases (a, -a) and (-a, a)
      GEN bval=stoi(-b);
      GEN disc=algnormdisc(alginit(Q, mkvec2(aval, bval), 0, 0));
      foundindef=setunion_i(foundindef, mkvec(disc));
      if(lg(foundindef)>lindef){lindef++;countindef[aind]++;}
      if(uissquare(b)) continue;//b is square, no need to do it.
      bval=stoi(b);
      disc=algnormdisc(alginit(Q, mkvec2(aval, bval), 0, 0));
      foundindef=setunion_i(foundindef, mkvec(disc));
      if(lg(foundindef)>lindef){lindef++;countindef[aind]++;}
    }
    GEN disc=algnormdisc(alginit(Q, mkvec2(aval, aval), 0, 0));//a=b
    foundindef=setunion_i(foundindef, mkvec(disc));
    if(lg(foundindef)>lindef){lindef++;countindef[aind]++;}
  }
  return mkvec4(founddef, foundindef, countdef, countindef);
}

//Counts the number of distinct algebras (a,b/Q) with max(|a|,|b|)<=N. Returns [defcount, indefcount], where each entry is a length N vecsmall, the ith entry representing the number of algebras with max(|a|,|b|)=N.
GEN alg_count_Q(long N, long prec){
  pari_sp top=avma;
  GEN dat=alg_count_Q_worker(cgetg(1, t_VEC), mkvec(gen_1), 1, N, prec);
  return gerepilecopy(top, mkvec2(gel(dat, 3), gel(dat, 4)));
}

//Takes the partial results from alg_count_Q_tofile, and computes alg_count_Q with the new N using this already found input.
GEN alg_count_Q_append(long N, char *oldfname, char *newfname, long prec){
  pari_sp top=avma;
  char *fullfile1=stack_sprintf("data/%s.dat", oldfname);
  GEN olddata=gp_readvec_file(fullfile1);
  GEN newdata=alg_count_Q_worker(gel(olddata, 2), gel(olddata, 3), itos(gel(olddata, 1))+1, N, prec);
  GEN defcount=shallowconcat(gel(olddata, 4), gel(newdata, 3));
  GEN indefcount=shallowconcat(gel(olddata, 5), gel(newdata, 4));
  char *fullfile2=stack_sprintf("data/%s.dat", newfname);
  FILE *f=fopen(fullfile2, "w");
  pari_fprintf(f, "%d\n%Ps\n%Ps\n%Ps\n%Ps\n", N, gel(newdata, 1), gel(newdata, 2), defcount, indefcount);
  fclose(f);
  return gerepilecopy(top, mkvec4(gel(newdata, 1), gel(newdata, 2), defcount, indefcount));
}

//Does alg_count_Q, with also presenting the results to the file data/fname.dat.
GEN alg_count_Q_tofile(long N, char *fname, long prec){
  pari_sp top=avma;
  if(!pari_is_dir("data")){
    int s=system("mkdir -p data");
    if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY data");
  }
  char *fullfile=stack_sprintf("data/%s.dat", fname);
  FILE *f=fopen(fullfile, "w");//Now we have created the output file f.
  GEN dat=alg_count_Q_worker(cgetg(1, t_VEC), mkvec(gen_1), 1, N, prec);
  pari_fprintf(f, "%d\n%Ps\n%Ps\n%Ps\n%Ps\n", N, gel(dat, 1), gel(dat, 2), gel(dat, 3), gel(dat, 4));
  fclose(f);
  return gerepilecopy(top, dat);
}



//Returns a better pair (a, b).
//GEN algbetterab_Q(GEN ab){
  //pari_sp top=avma;
  //GEN a=gel(ab, 1), b=gel(ab, 2);
//}
