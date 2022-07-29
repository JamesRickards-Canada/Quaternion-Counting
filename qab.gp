\\print("\n\nType '?qab' for help.\n\n");
\\addhelp(qab, "This package is for testing optimizing (a, b) in a quaternion algebra.\n For each subtopic ``P (p)'', call ?p to access a basic description and list of methods. Subtopics:\n ");

\\BASIC METHODS
	install(ab_disc,"GD0,G,",,"./libqab.so");
	addhelp(ab_disc,"Inputs a, b.\n Returns the discriminant of the quaternion algebra (a,b/Q). We can input a as [a,b], and can pass in the factorizations of a and/or b instead.");
	install(ab_ramprimes,"GD0,G,",,"./libqab.so");
	addhelp(ab_ramprimes,"Inputs a, b.\n Returns the set of finite primes ramifying in the quaternion algebra (a,b/Q). We can input a as [a,b], and can pass in the factorizations of a and/or b instead.");
	install(algabfromram,"GGG",,"./libqab.so");
	addhelp(algabfromram,"Input F, D, infram.\n Returns a pair (a, b) such that (a, b/F) is a quaternion algebra with N_(F/Q)(disc)=D and infinite ramification specified by the length [F:Q] vector infram (0=unram, 1=ram), if it exists. Otherwise, returns 0. If there are multiple such algebras, we only return one of them.");
	install(alggetab,"G",,"./libqab.so");\
	addhelp(alggetab,"Input A, a quaternion algebra.\n Returns the pair (a, b) such that A=(a,b/F).");
	install(algnormdisc,"G","algnormdisc","./libqab.so");
	addhelp(algnormdisc,"Input A, an algebra.\n Returns the norm to Q of the discriminant of A.");
	install(algramifiedplacesf,"G","algramifiedplacesf","./libqab.so");
	addhelp(algramifiedplacesf,"Input A, an algebra.\n Returns the vector of finite places that ramify.");
	install(hilbertss,"lLLU",,"./libqab.so");
	addhelp(hilbertss,"Inputs x, y, p with p>0.\n Returns the Hilbert symbol (x, y)_p, where x and y are longs, and p is a ulong.");
	install(vecsquarefreeu_split,"UU",,"./libqab.so");
	addhelp(vecsquarefreeu_split,"Inputs a, b, positive integers.\n Returns [def, indef], vecsmalls of the squarefree integers between a and b, where the entries in def have odd numbers of prime factors, and indef have even numbers of prime factors.");

\\COUNTING ALGEBRAS
	install(alg_count_alg,"GUD0,U,",,"./libqab.so");
	addhelp(alg_count_alg,"Inputs disc, N1, N2.\n Counts instances of the algebra with disc disc among N1<=|a|<=N2 and 1<=|b|<=|a|. We only count one of (a,b) and (b,a), and restrict to a,b being squarefree. We can also pass in disc as a Vec/Vecsmall of discriminants, and N2=0 means the range [1, N1].");
	install(alg_count_alg_tofile,"vGUUs",,"./libqab.so");
	addhelp(alg_count_alg_tofile,"Input discs, N1, N2, fname.\n Does alg_count_alg, and writes the results to data/fname.dat.");
	install(alg_count_Q,"U",,"./libqab.so");
	addhelp(alg_count_Q,"Input N, a postive integer.\n Counts the number of distinct algebras (a,b/Q) with max(|a|,|b|)<=N. Returns [defcount, indefcount], where each entry is a length N vecsmall, the ith entry representing the number of algebras with max(|a|,|b|)=N.");
	install(alg_count_Q_append,"vUss",,"./libqab.so");
	addhelp(alg_count_Q_append,"Inputs N, oldfname, newfname.\n Does alg_count_Q, finding and writing [founddef, foundindef, countdef, countindef] to the file data/newfname.dat. We assume that the file data/oldfname.dat already includes partial computations, i.e. for a smaller N.");
	install(alg_count_Q_tofile,"vUs",,"./libqab.so");
	addhelp(alg_count_Q_tofile,"Inputs N, fname.\n Does alg_count_Q, finding and writing [founddef, foundindef, countdef, countindef] to the file data/fname.dat.");

\\OTHER METHODS
	install(veccumu,"G",,"./libqab.so");
	addhelp(veccumu,"Input v, a vector/vecsmall.\n Returns w, where w[i]=v[1]+v[2]+...+v[i].");
	install(writevecs,"vGs",,"./libqab.so");
	addhelp(writevecs,"Inputs v, fname.\n Writes v to data/fname.dat, one entry per row. v can also be a vector of vectors/vecsmalls, all of the same type and length. In this case, each row contains v[1][i], v[2][i], ..., v[n][i] separated by spaces.");


\\visual.c
	\\DATA
		install("integerbin","GGD0,G,",,"./libqab.so");
		addhelp(integerbin,"Inputs v, binlen, {binstart=0}.\n Assumes v is a sorted list of integers, and puts them into bins of length binlen, starting with binstart (assumed to be 0). Returns [binends, counts], with binends being the last number in the bin.");
		install("integerbin_cumu","GGD0,G,",,"./libqab.so");
		addhelp(integerbin_cumu,"Inputs v, binlen, {binstart=0}.\n Assumes v is a sorted list of integers, and puts them into bins of length binlen, starting with binstart (assumed to be 0). Returns [binends, counts], with binends being the last number in the bin. This is cumulative, so counts is increasing.");
		install("veccount","G",,"./libqab.so");
		addhelp(veccount,"Input v, a vector.\n Returns [uniq, count], where uniq is the sorted vector v with repeats removed, and count is the corresponding number of times they appear in v.");
		install("vecsmallcount","G",,"./libqab.so");
		addhelp(vecsmallcount,"Input v, a vecsmall.\n Returns [uniq, count], where uniq is the sorted vecsmall v with repeats removed, and count is the corresponding number of times they appear in v.");
		install("ZV_countnonpos","lG",,"./libqab.so");
		addhelp(ZV_countnonpos,"Input v, a sorted vector of integers.\n Returns the number of entries that are nonpositive.");
		
	\\HISTOGRAMS
		install("hist_make","GrD0,L,D0,L,Drp",,"./libqab.so");
		addhelp(hist_make,"Inputs data, imagename, {compilenew=0}, {open=0}, {plotoptions=NULL}: sorted list of real numbers data, name of the tikz picture, {compilenew=0, 1}, {open=0, 1}, {plotoptions=NULL, string}.\n Automatically bins the data, and creates a pdf of the histogram using tikz and externalize. The output is in the folder /images, with the build file (named imagename_build.tex) being in the folder /images/build. If compilenew=0 assumes the LaTeX document to compile the plot is pre-made, and otherwise this method automatically writes it. If additionally, plotoptions!=NULL, this string is inserted in between \\begin{axis} and \\end{axis} in the LaTeX document (allowing one to customize how the histogram looks). If open=1, the pdf is automatically opened (only works with Linux subsystem for Windows). The returned value is used to modify the histogram, e.g. changing the bins, scaling it, and changing the range.");
		install("hist_rebin","GGGp",,"./libqab.so");
		addhelp(hist_rebin,"Inputs data, histdata, nbins: the sorted data, the length 7 vector output of a hist_ method, the number of bins.\n Rebins the data according to the new number of bins, and updates histdata.");
		install("hist_rerange","GGGGp",,"./libqab.so");
		addhelp(hist_rerange,"Inputs data, histdata, minx, maxx: the sorted data, the length 7 vector output of a hist_ method, minimum x-value, maximum x-value.\n Rebins the data according to the new minimum and maximum value. This is useful when there are outliers that skew the look of the graph. Returns the updated histdata.");
		install("hist_rescale","GGLp",,"./libqab.so");
		addhelp(hist_rescale,"Inputs data, histdata, scale: the sorted data, the length 7 vector output of a hist_ method, and scale=0, 1.\n If scale=1 scales the data so the total area is 1, and if scale=0 uses the absolute count for the y-axis. Returns the updated histdata.");

	\\REGRESSIONS & PLOTS
		install("OLS","GGD1,L,",,"./libqab.so");
		addhelp(OLS,"Inputs X, y, {retrsqr=1}:  m*n matrix X with top row being all 1's, length n column vector y, retrsqr=0, 1.\n Performs ordinary least squares regression on the data, where the n inputs are the columns of X, and the outputs are the entries of y. We must include a constant term, hence why the first row of X must be all 1's. If retrsqr=1, returns [pararms, R^2], and otherwise returns params, where params is the length m column vector of best fit parameters.");
		install("OLS_nointercept","GGD1,L,",,"./libqab.so");
		addhelp(OLS_nointercept,"Inputs X, y, {retrsqr=1}: vector X, column vector y (of same length), retrsqr=0, 1.\n Performs ordinary least squares regression on the data assuming that y[i]=c*X[i], i.e. the y-intercept is 0. Returns c if retrsqr=0, or [c, R^2] otherwise.");
		install("OLS_single","GGD1,L,",,"./libqab.so");
		addhelp(OLS_single,"Inputs x, y, {retrsqr=1}: vector x, column vector y, retrsqr=0, 1. Performs linear regression for a single variable (essentially a macro for OLS with y=mx+b.");
		install("rsquared","GGG",,"./libqab.so");
		addhelp(rsquared,"Inputs X, y, fit: X and y data supplied to OLS, and fit the proposed fit (a column vector of parameters). This returns the R^2 value for this proposal.");

	\\TEX
		install("tex_recompile","vG",,"./libqab.so");
		addhelp(tex_recompile,"Input data, the output of a tex image creation call.\n Recompiles the image, returning nothing. This is used when you edit the LaTeX document by hand.");

	\\GENERAL HELP
		addhelp(visual,"This package deals with visualizing data. Subtopics:\n Data (data)\n Histograms (hist)\n Regressions/plots (reg)\n Tex (tex)");
		addhelp(data,"integerbin, veccount, vecsmallcount.");
		addhelp(hist,"Installed methods:\n hist_make, hist_rebin, hist_rerange, hist_rescale.");
		addhelp(reg,"OLS, OLS_nointercept, OLS_single, rsquared.");
		addhelp(tex,"tex_recompile");
	
default(parisize, "4096M");\\Must come at the end