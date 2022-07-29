//BASIC METHODS
GEN ab_disc(GEN a, GEN b);
GEN ab_ramprimes(GEN a, GEN b);
GEN algabfromram(GEN F, GEN D, GEN infram);
GEN alggetab(GEN A);
GEN algnormdisc(GEN A);
GEN algramifiedplacesf(GEN A);
long hilbertss(long x, long y, ulong p);

//COUNTING ALGEBRAS
GEN alg_count_alg(GEN disc, ulong N1, ulong N2);
void alg_count_alg_tofile(GEN disc, ulong N1, ulong N2, char *fname);
GEN alg_count_Q(ulong N);
void alg_count_Q_append(ulong N, char *oldfname, char *newfname);
void alg_count_Q_tofile(ulong N, char *fname);

//OTHER METHODS
GEN veccumu(GEN v);