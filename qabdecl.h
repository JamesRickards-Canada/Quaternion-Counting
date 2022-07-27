//BASIC METHODS
GEN ab_disc(GEN a, GEN b);
GEN ab_ramprimes(GEN a, GEN b);
GEN algabfromram(GEN F, GEN D, GEN infram);
GEN alggetab(GEN A);
GEN algnormdisc(GEN A);
GEN algramifiedplacesf(GEN A);

//COUNTING ALGEBRAS

GEN alg_count_Q(long N, long prec);
GEN alg_count_Q_append(long N, char *oldfname, char *newfname, long prec);
GEN alg_count_Q_tofile(long N, char *fname, long prec);