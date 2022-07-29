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



//visual.c

//DATA
GEN integerbin(GEN v, GEN binlen, GEN binstart);
GEN integerbin_cumu(GEN v, GEN binlen, GEN binstart);
GEN veccount(GEN v);
GEN vecsmallcount(GEN v);
long ZV_countnonpos(GEN v);

//HISTOGRAMS
void hist_autocompile(GEN minx, GEN maxx, char *imagename, char *plotoptions, int open);
GEN hist_make(GEN data, char *imagename, int compilenew, int open, char *plotoptions, long prec);
GEN hist_tobins(GEN data, GEN minx, GEN maxx, GEN nbins, int toscale, int compilenew, char *imagename, char *plotoptions, int open, long prec);
GEN hist_tobins_defaultbins(GEN data, GEN minx, GEN maxx, int toscale, int compilenew, char *imagename, char *plotoptions, int open, long prec);
GEN hist_rebin(GEN data, GEN histdata, GEN nbins, long prec);
GEN hist_rerange(GEN data, GEN histdata, GEN minx, GEN maxx, long prec);
GEN hist_rescale(GEN data, GEN histdata, int scale, long prec);

//REGRESSIONS
GEN OLS(GEN X, GEN y, int retrsqr);
GEN OLS_nointercept(GEN X, GEN y, int retrsqr);
GEN OLS_single(GEN x, GEN y, int retrsqr);
GEN rsquared(GEN X, GEN y, GEN fit);

//TEX
GEN tex_makecolours(int ncol);
void tex_compile(char *imagename, int open);
void tex_recompile(GEN data);