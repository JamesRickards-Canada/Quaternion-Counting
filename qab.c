//(a,b) optimization testing.


//INCLUSIONS

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#include "pari/paripriv.h"

#ifndef METHDECL
#define METHDECL
#include "qabdecl.h"
#endif

//STATIC DECLARATIONS
static ulong ab_discu(long a, long b, GEN apdivs, GEN bpdivs);
static int gomel(long t);
static GEN alg_count_alg_worker(GEN disc, ulong N1, ulong N2);
static GEN alg_count_Q_hash(GEN founddef, GEN foundindef, ulong N1, ulong N2);
static int gp_read_stream_buf(FILE *fi, Buffer *b);

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

//Returns the discriminant of (a,b/Q). apdivs and bpdivs are vecsmalls of the prime divisors
static ulong ab_discu(long a, long b, GEN apdivs, GEN bpdivs){
  ulong disc=1;
  long aind=1, bind=1, la=lg(apdivs), lb=lg(bpdivs);
  while(aind<la && bind<lb){//Loop through until we run of of primes from either a or b.
    ulong p1=apdivs[aind], p2=bpdivs[bind];
    ulong p=(p1>p2) ? p2:p1;//p=min(p1, p2)
    if(hilbertss(a, b, p)==-1) disc=disc*p;
    if(p==p1) aind++;
    if(p==p2) bind++;
  }
  while(aind<la){//More primes in a to check
    ulong p=apdivs[aind];
    if(hilbertss(a, b, p)==-1) disc=disc*p;
    aind++;
  }
  while(bind<lb){//More primes in b to check
    ulong p=bpdivs[bind];
    if(hilbertss(a, b, p)==-1) disc=disc*p;
    bind++;
  }
  if((a&1) && (b&1) && hilbertss(a, b, 2)==-1) disc=disc*2;//2 ramification when a and b are odd.
  return disc;//No garbage!!
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

/* t a t_INT, is t = 3,5 mod 8 ? */
static int gomel(long t){
  if(t==0) return 0;
  if(t<0) t=-t;
  switch(t&7){//t mod 8
    case 3:
    case 5: return 1;
    default: return 0;
  }
}

//Adapted from hilbertii.
long hilbertss(long x, long y, ulong p){
  pari_sp top=avma;
  if(x==0 || y==0) return 0;
  if(p==1) pari_err_PRIME("hilbertii", stoi(p));
  long xrem, yrem, z;
  long oddvx=odd(z_lvalrem(x, p, &xrem));
  long oddvy=odd(z_lvalrem(y, p, &yrem));
  if(p==2){/* x, y are p-units, compute hilbert(xrem * p^oddvx, yrem * p^oddvy, p) */
    z=((xrem&3)==3 && (yrem&3)==3)? -1: 1;
    if(oddvx && gomel(yrem)) z=-z;
    if(oddvy && gomel(xrem)) z=-z;
  }
  else{
    z=(oddvx && oddvy && (p&3)==3)? -1: 1;
    if(oddvx && kross(yrem, p)<0) z=-z;
    if(oddvy && kross(xrem, p)<0) z=-z;
  }
  return gc_long(top, z);
}

//vecsquarefreeu, but splits the return into [def, indef] (each entry a vecsmall).
GEN vecsquarefreeu_split(ulong a, ulong b){
  pari_sp top=avma;
  GEN all=vecfactorsquarefreeu(a, b);
  long lgall=lg(all);
  GEN def=vecsmalltrunc_init(lgall), indef=vecsmalltrunc_init(lgall);
  for(long i=1;i<lgall;i++){
	if(!gel(all, i)) continue;
	if(lg(gel(all, i))%2==0) vecsmalltrunc_append(def, a+i-1);//Odd number of prime divisors, definite
	else vecsmalltrunc_append(indef, a+i-1);
  }
  return gerepilecopy(top, mkvec2(def, indef));
}


//COUNTING ALGEBRAS

//Counts instances of the algebra with disc disc among N1<=|a|<=N2 and |b|<=|a|. We only count one of (a,b) and (b,a), and restrict to a,b being squarefree. We can also pass in a set of discriminants as a Vec/Vecsmall.
GEN alg_count_alg(GEN disc, ulong N1, ulong N2){
  pari_sp top=avma;
  if(N2==0){N2=N1;N1=1;}
  return gerepilecopy(top, alg_count_alg_worker(disc, N1, N2));
}

//alg_count_alg, with writing to a file.
void alg_count_alg_tofile(GEN disc, ulong N1, ulong N2, char *fname){
  pari_sp top=avma;
  if(!pari_is_dir("data")){
    int s=system("mkdir -p data");
    if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY data");
  }
  char *fullfile=stack_sprintf("data/%s.dat", fname);
  FILE *f=fopen(fullfile, "w");//Now we have created the output file f.
  GEN dat=alg_count_alg_worker(disc, N1, N2);
  if(typ(dat)==t_VECSMALL) pari_fprintf(f, "%Ps\n", dat);
  else for(long i=1;i<lg(dat);i++) pari_fprintf(f, "%Ps\n", gel(dat, i));
  fclose(f);
  set_avma(top);
}

//Counts instances of the algebra with disc disc among |a|,|b|<=N. We only count one of (a,b) and (b,a), and restrict to a,b being squarefree. We can also pass in a set of discriminants as a Vec/Vecsmall. NOT stack clean.
static GEN alg_count_alg_worker(GEN disc, ulong N1, ulong N2){
  int isvec=1;
  long t=typ(disc);
  if(t==t_INT){isvec=0;disc=mkvecsmall(itos(disc));}
  else if(t==t_VEC) disc=gtovecsmall(disc);
  ulong lgdisc=lg(disc), lN=N2-N1+1, aind=0;
  GEN countset=cgetg(lgdisc, t_VEC);
  for(long i=1;i<lgdisc;i++) gel(countset, i)=zero_zv(lN);
  GEN facts=vecfactorsquarefreeu(1, N2);//Find prime divisors of all numbers between 1 and N.
  for(long a=N1;a<=N2;a++){
    aind++;
    if(a%100==0) pari_printf("a=%d done\n", a);
    GEN apdivs=gel(facts, a);
    if(!apdivs) continue;//WLOG squarefree
    for(long b=1;b<=a;b++){
      GEN bpdivs=gel(facts, b);
      if(!bpdivs) continue;//WLOG squarefree
      ulong foundd=ab_discu(-a, -b, apdivs, bpdivs);
      for(long i=1;i<lgdisc;i++){if(foundd==disc[i]){gel(countset, i)[aind]++;break;}}
      foundd=ab_discu(a, -b, apdivs, bpdivs);
      for(long i=1;i<lgdisc;i++){if(foundd==disc[i]){gel(countset, i)[aind]++;break;}}
      foundd=ab_discu(a, b, apdivs, bpdivs);
      for(long i=1;i<lgdisc;i++){if(foundd==disc[i]){gel(countset, i)[aind]++;break;}}
      if(a==b) continue;
      foundd=ab_discu(-a, b, apdivs, bpdivs);
      for(long i=1;i<lgdisc;i++){if(foundd==disc[i]){gel(countset, i)[aind]++;break;}}
    }
  }
  if(!isvec) return gel(countset, 1);
  return countset;
}

//Counts the number of new and distinct algebras (a,b/Q) with N1<=max(|a|,|b|)<=N2. Returns [founddef, foundindef, defcount, indefcount], where the last two entries are length N2-N1+1 vecsmalls, the ith entry representing the number of algebras with max(|a|,|b|)=N1+i-1. This is not stack clean.
static GEN alg_count_Q_hash(GEN founddef, GEN foundindef, ulong N1, ulong N2){
  long lN=N2-N1+1;//Tracks the number of algebras found as Vecsmalls
  hashtable *hdef=hash_create_ulong(100, 1);//Make the hashtable using the stack. For definite discriminants
  hashtable *hindef=hash_create_ulong(100, 1);//For indefinite discriminants
  for(long i=1;i<lg(founddef);i++) hash_insert(hdef, (void *)founddef[i], NULL);//Inserting the found discriminants into our hash
  for(long i=1;i<lg(foundindef);i++) hash_insert(hindef, (void *)foundindef[i], NULL);
  GEN facts=vecfactorsquarefreeu(1, N2);//Find prime divisors of all numbers between 1 and N2.
  GEN countdef=zero_zv(lN), countindef=zero_zv(lN);//Number of new algebras for each N1<=n<=N2
  long aind=0;
  for(long a=N1;a<=N2;a++){//we do -a first, then a
    aind++;
    if(a%100==0) pari_printf("a=%d done\n", a);
    GEN apdivs=gel(facts, a);
    if(!apdivs) continue;//WLOG squarefree
    for(long b=1;b<=a;b++){//We do -b first, then b
      GEN bpdivs=gel(facts, b);
      if(!bpdivs) continue;//WLOG squarefree
      ulong disc=ab_discu(-a, -b, apdivs, bpdivs);
      if(!hash_search(hdef, (void *)disc)){//New disc!
        countdef[aind]++;
        hash_insert(hdef, (void *)disc, NULL);
      }
      disc=ab_discu(-a, b, apdivs, bpdivs);
      if(!hash_search(hindef, (void *)disc)){//New disc!
        countindef[aind]++;
        hash_insert(hindef, (void *)disc, NULL);
      }
      disc=ab_discu(a, b, apdivs, bpdivs);
      if(!hash_search(hindef, (void *)disc)){//New disc!
        countindef[aind]++;
        hash_insert(hindef, (void *)disc, NULL);
      }
      if(a==b) continue;//If a=b, (a,-b) and (-a,b) give the same algebra.
      disc=ab_discu(a, -b, apdivs, bpdivs);
      if(!hash_search(hindef, (void *)disc)){//New disc!
        countindef[aind]++;
        hash_insert(hindef, (void *)disc, NULL);
      }
    }
  }
  GEN defalg=hash_keys(hdef);
  GEN indefalg=hash_keys(hindef);
  hash_destroy(hdef);
  hash_destroy(hindef);
  return mkvec4(defalg, indefalg, countdef, countindef);
}

//Counts the number of distinct algebras (a,b/Q) with max(|a|,|b|)<=N. Returns [defcount, indefcount], where each entry is a length N vecsmall, the ith entry representing the number of algebras with max(|a|,|b|)=N.
GEN alg_count_Q(ulong N){
  pari_sp top=avma;
  GEN dat=alg_count_Q_hash(cgetg(1, t_VECSMALL), cgetg(1, t_VECSMALL), 1, N);
  return gerepilecopy(top, mkvec2(gel(dat, 3), gel(dat, 4)));
}

//Takes the partial results from alg_count_Q_tofile, and computes alg_count_Q with the new N using this already found input.
void alg_count_Q_append(ulong N, char *oldfname, char *newfname){
  pari_sp top=avma;
  char *fullfile1=stack_sprintf("data/%s_defalg.dat", oldfname);
  FILE *f=fopen(fullfile1, "r");
  GEN olddefalg=readvecsmall(f);
  fclose(f);
  char *fullfile2=stack_sprintf("data/%s_indefalg.dat", oldfname);
  f=fopen(fullfile2, "r");
  GEN oldindefalg=readvecsmall(f);
  fclose(f);
  char *fullfile3=stack_sprintf("data/%s_defcount.dat", oldfname);
  f=fopen(fullfile3, "r");
  GEN olddefcount=readvecsmall(f);
  fclose(f);
  char *fullfile4=stack_sprintf("data/%s_indefcount.dat", oldfname);
  f=fopen(fullfile4, "r");
  GEN oldindefcount=readvecsmall(f);
  fclose(f);
  GEN newdata=alg_count_Q_hash(olddefalg, oldindefalg, lg(olddefcount), N);
  GEN defcount=shallowconcat(olddefcount, gel(newdata, 3));
  GEN indefcount=shallowconcat(oldindefcount, gel(newdata, 4));
  
  char *fullfile5=stack_sprintf("data/%s_defalg.dat", newfname);
  f=fopen(fullfile5, "w");
  for(long i=1;i<lg(gel(newdata, 1));i++) pari_fprintf(f, "%d\n", gel(newdata, 1)[i]);//Write out the definite algebras
  fclose(f);
  char *fullfile6=stack_sprintf("data/%s_indefalg.dat", newfname);
  f=fopen(fullfile6, "w");
  for(long i=1;i<lg(gel(newdata, 2));i++) pari_fprintf(f, "%d\n", gel(newdata, 2)[i]);//Write out the indefinite algebras
  fclose(f);
  char *fullfile7=stack_sprintf("data/%s_defcount.dat", newfname);
  f=fopen(fullfile7, "w");
  for(long i=1;i<lg(defcount);i++) pari_fprintf(f, "%d\n", defcount[i]);//Write out the definite counts
  fclose(f);
  char *fullfile8=stack_sprintf("data/%s_indefcount.dat", newfname);
  f=fopen(fullfile8, "w");
  for(long i=1;i<lg(indefcount);i++) pari_fprintf(f, "%d\n", indefcount[i]);//Write out the indefinite counts
  fclose(f);
  set_avma(top);
}

//Does alg_count_Q, with also presenting the results to the file data/fname.dat.
void alg_count_Q_tofile(ulong N, char *fname){
  pari_sp top=avma;
  if(!pari_is_dir("data")){
    int s=system("mkdir -p data");
    if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY data");
  }
  GEN dat=alg_count_Q_hash(cgetg(1, t_VECSMALL), cgetg(1, t_VECSMALL), 1, N);
  char *fullfile1=stack_sprintf("data/%s_defalg.dat", fname);
  FILE *f=fopen(fullfile1, "w");//Now we have created the output file f.
  for(long i=1;i<lg(gel(dat, 1));i++) pari_fprintf(f, "%d\n", gel(dat, 1)[i]);//Write out the definite algebras
  fclose(f);
  char *fullfile2=stack_sprintf("data/%s_indefalg.dat", fname);
  f=fopen(fullfile2, "w");
  for(long i=1;i<lg(gel(dat, 2));i++) pari_fprintf(f, "%d\n", gel(dat, 2)[i]);//Write out the indefinite algebras
  fclose(f);
  char *fullfile3=stack_sprintf("data/%s_defcount.dat", fname);
  f=fopen(fullfile3, "w");
  for(long i=1;i<lg(gel(dat, 3));i++) pari_fprintf(f, "%d\n", gel(dat, 3)[i]);//Write out the definite counts
  fclose(f);
  char *fullfile4=stack_sprintf("data/%s_indefcount.dat", fname);
  f=fopen(fullfile4, "w");
  for(long i=1;i<lg(gel(dat, 4));i++) pari_fprintf(f, "%d\n", gel(dat, 4)[i]);//Write out the indefinite counts
  fclose(f);
  set_avma(top);
}





//OTHER METHODS

//Given v, returns a vector/vecsmall w (same type as v) where w[i]=v[1]+v[2]+...+v[i].
GEN veccumu(GEN v){
  long lg;
  GEN w=cgetg_copy(v, &lg);
  if(lg==1) return w;//Empty
  if(typ(v)==t_VEC){
    gel(w, 1)=gcopy(gel(v, 1));
    for(long i=2;i<lg;i++) gel(w, i)=gadd(gel(v, i), gel(w, i-1));
  }
  else{//Vecsmall
    w[1]=v[1];
    for(long i=2;i<lg;i++) w[i]=v[i]+w[i-1];
  }
  return w;
}

//Write v to data/fname.dat, one entry per line. If v is a vector of vectors/vecsmalls, assume they are all the same length AND type, and writes them in columns separated by spaces
void writevecs(GEN v, char *fname){
  long t=typ(v), lv=lg(v);
  if(lv==1) return;
  if(!pari_is_dir("data")){
    int s=system("mkdir -p data");
    if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY data");
  }
  char *fullfile=stack_sprintf("data/%s.dat", fname);
  FILE *f=fopen(fullfile, "w");//Now we have created the output file f.
  if(t==t_VECSMALL){
    for(long i=1;i<lv;i++) pari_fprintf(f, "%d\n", v[i]);
  }
  else{//Vector
    t=typ(gel(v, 1));
    if(t==t_VECSMALL){
      for(long i=1;i<lg(gel(v, 1));i++){
        for(long j=1;j<lv-1;j++) pari_fprintf(f, "%d ", gel(v, j)[i]);
        pari_fprintf(f, "%d\n", gel(v, lv-1)[i]);
      }
    }
    else if(t==t_VEC){
      for(long i=1;i<lg(gel(v, 1));i++){
        for(long j=1;j<lv-1;j++) pari_fprintf(f, "%Ps ", gmael(v, j, i));
        pari_fprintf(f, "%Ps\n", gmael(v, lv-1, i));
      }
    }
    else{
      for(long i=1;i<lv;i++) pari_fprintf(f, "%Ps\n", gel(v, i));
    }
  }
  fclose(f);
}

//Copied from es.c, needed for readvecsmall
static int gp_read_stream_buf(FILE *fi, Buffer *b){
  input_method IM;
  filtre_t F;
  init_filtre(&F, b);
  IM.file=(void*)fi;
  IM.myfgets=(fgets_t)&fgets;
  IM.getline=&file_input;
  IM.free=0;
  return input_loop(&F, &IM);
}

//Adapted from gp_readvec_stream in es.c
GEN readvecsmall(FILE *fi){
  pari_sp top=avma;
  Buffer *b=new_buffer();
  long i=1, n=16;
  GEN z=cgetg(n+1, t_VECSMALL);
  for(;;){
    if(!gp_read_stream_buf(fi, b)) break;
    if(!*(b->buf)) continue;
    if(i>n){
      n<<=1;
      z=vecsmall_lengthen(z, n);
    }
    z[i++] = itos(readseq(b->buf));
  }
  setlg(z, i);
  delete_buffer(b);
  return gerepilecopy(top, z);
}

