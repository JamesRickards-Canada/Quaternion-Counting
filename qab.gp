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


\\COUNTING ALGEBRAS
	install(alg_count_Q,"Up",,"./libqab.so");
	addhelp(alg_count_Q,"Input N, a postive integer.\n Counts the number of distinct algebras (a,b/Q) with max(|a|,|b|)<=N. Returns [defcount, indefcount], where each entry is a length N vecsmall, the ith entry representing the number of algebras with max(|a|,|b|)=N.");
	install(alg_count_Q_append,"vUssp",,"./libqab.so");
	addhelp(alg_count_Q_append,"Inputs N, oldfname, newfname.\n Does alg_count_Q, finding and writing [founddef, foundindef, countdef, countindef] to the file data/newfname.dat. We assume that the file data/oldfname.dat already includes partial computations, i.e. for a smaller N.");
	install(alg_count_Q_tofile,"vUsp",,"./libqab.so");
	addhelp(alg_count_Q_tofile,"Inputs N, fname.\n Does alg_count_Q, finding and writing [founddef, foundindef, countdef, countindef] to the file data/fname.dat.");

\\VISUAL
	\\Files from another project, just comment them out if needed
	install(OLS,"GGD1,L,",,"../Apollonian/libapol.so");
	addhelp(OLS,"Inputs X, y, {retrsqr=1}:  m*n matrix X with top row being all 1's, length n column vector y, retrsqr=0, 1.\n Performs ordinary least squares regression on the data, where the n inputs are the columns of X, and the outputs are the entries of y. We must include a constant term, hence why the first row of X must be all 1's. If retrsqr=1, returns [pararms, R^2], and otherwise returns params, where params is the length m column vector of best fit parameters.");
	install(OLS_nointercept,"GGD1,L,",,"../Apollonian/libapol.so");
	addhelp(OLS_nointercept,"Inputs X, y, {retrsqr=1}: vector X, column vector y (of same length), retrsqr=0, 1.\n Performs ordinary least squares regression on the data assuming that y[i]=c*X[i], i.e. the y-intercept is 0. Returns c if retrsqr=0, or [c, R^2] otherwise.");
	
default(parisize, "4096M");\\Must come at the end