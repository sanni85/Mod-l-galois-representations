def lookforec(l,n,a,B):
        """
        Input:

        l is 2, 3 or 5.

        n is a positive integer, which is the Artin conductor of the
        representation.

        a is the list of traces in the format [p, a_p] for p prime.

        B is a bound set for computations

        Output: a list of labels of elliptic curves over QQ.  The mod
        l representation associated to each elliptic curve in the list
        is isomorphic to the representation whose traced are in the
        list a
        """
	F = GF(l)

        h = prod([p for (p,ap) in a if p<=B and (ap^2-(p+1)^2)%l == 0], 1)

	# this is the level raising condition, the bound B is a stopping parameter for the search
	# the code below sets the right power of l in the search for the conductor of the elliptic curves over QQ

        max_exp = 8 if l==2 else 5 if l==3 else 2
        h = lcm(h,l^max_exp)

        conductors = [n*i for i in h.divisors() if n*i< CremonaDatabase().largest_conductor() ]

	res=[] # will hold the output list of labels

	for E in cremona_optimal_curves(conductors):
		wt,cd = serre_invariants(E, l)
		if n%cd==0:
			m=E.conductor()/n
			avoid=prime_divisors(m)
			primes_to_check=[ p for p in prime_range(100) if p not in avoid]
			ap_list = [[p, F(E.ap(p))] for p in primes_to_check]
                	a_check = [i for i in a if i[0] in primes_to_check]

			# here is the check outside the ramified primes
			if ap_list == a_check:
				# here the check for ramified primes
				cont = 0
				for (p,ap) in a:
					if p in avoid:
                	                	Eap = F(E.ap(p))
                        	        	if F(Eap^2-ap*Eap+p):
							cont=1
                                        	break
				if cont==0:
					res.append(E);
	cds=set([serre_invariants(E, l)[1] for E in res])
	print "the smallest conductor is", min(cds);
	for i in cds:
		print "conductor ", i; 
		res1=[E.label()for E in res if serre_invariants(E, l)[1]==i]
		print res1;
	return [E.label()for E in res]
