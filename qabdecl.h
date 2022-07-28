//BASIC METHODS
GEN ab_disc(GEN a, GEN b);
GEN ab_ramprimes(GEN a, GEN b);
GEN algabfromram(GEN F, GEN D, GEN infram);
GEN alggetab(GEN A);
GEN algnormdisc(GEN A);
GEN algramifiedplacesf(GEN A);
long hilbertss(long x, long y, ulong p);

//COUNTING ALGEBRAS

GEN alg_count_Q(ulong N, long prec);
GEN alg_count_Q_append(long N, char *oldfname, char *newfname, long prec);
GEN alg_count_Q_tofile(long N, char *fname, long prec);