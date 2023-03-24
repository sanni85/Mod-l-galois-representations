# This code is an implementation of the main theorems in "Détermination du poids et du conducteur associés aux représentations des points de p-torsion d'une courbe elliptique" by Alain Kraus in Rozprawy Matematycznetomnrw serii: 364 wydano: 1997. 
#Author: Samuele Anni




def case2(D):    
    #This is an auxiliary function that returns the weight for mod 2 representations associated to the elliptic curve over QQ which has additive reduction at 2. Reference: Théorème 3 of Kraus' paper.
    #Input: D the discriminant of the minimal model of the elliptic curve over QQ
    vd = valuation(D,2)
    krho = 2 if vd%2==0 else 4
    return krho          


def case3(D,j,c4,c6):
    #This is an auxiliary function that returns the weight for mod 3 representations associated to the elliptic curve over QQ which has additive reduction at 3. Reference: Théorème 2 of Kraus' paper. 
    #Input: D is the discriminant, j the j-invariant, c4 and c6 the C4 and C6 invariants of the minimal model of the elliptic curve over QQ
    vd=valuation(D,3);
    D1=(D*3^(-vd));
    if valuation(j,3) < 0 :
        if vd % 3 == 0 :
            if (D1 % 9 == 1) or (D1 % 9 == 8) :
                krho=2;
            else:
                krho=6;
        else:
            krho=8;
    else:
        if vd == 3 :
            if valuation(c6,3) == 3 :
                if (D1 % 9 == 1) or (D1 % 9 == 8) :
                    krho=2;
                else:
                    krho=6;
            else:
                krho=6;   
        elif vd == 4 : 
            krho=8;
        elif vd == 5 :
            if valuation(c6,3) == 3 :
                krho=8;
            else:
                krho=4;
        elif vd == 6 :
            if valuation(c6,3) == 3 :
                if (D1 % 9 == 1) or (D1 % 9 == 8) :
                    krho=2;
                else:
                    krho=6;
            else:
                krho=6;   
        elif vd == 7 :
            krho=8;   
        elif vd == 9 :
            krho=2; 
        elif vd == 10 : 
            krho=4;
        elif vd == 11 : 
            if valuation(c6,3) == 6 :
                krho=4;
            else:
                krho=8;
        elif vd == 12 :
            krho=2;
        else:
            krho=4;
    return krho;    


def casep(p,vj,vd,vc4,vc6):
    # This is an auxiliary function that returns the weight for % p representations associated to the elliptic curve over QQ which has additive reduction at p. Reference: Théorème 1 of Kraus' paper. 
    # p is a prime p>3, vj, vd, vc4, vc6 are the p-adic valuations of the j-invariant, the  discriminant, the C4 and C6 invariants of the minimal %el of the elliptic curve over QQ
    if vj < 0 :
        if vj % p == 0 :
            krho=(p^2+3)/2;
        else:   
            krho=(p+1)^2/2; 
    else:
        if vd == 2 :
            if p == 7 :
                if vc4 == 1 :
                    krho=14;
                else:
                    krho=2;
            elif p % 3 == 1 and not p == 7 :
                krho=(p^2+4*p+7)/6;
            else:
                krho=(p^2+6*p+5)/6;
        elif vd == 3 :
            if p == 5 :
                if vc6 == 2 :
                    krho=10;
                else:
                    krho=2;
            else:
                if p % 4 == 1 :
                    krho=(p^2+2*p+5)/4;
                else:
                    krho=(p^2+4*p+3)/4;
        elif vd == 4 :
            if p % 3 == 1 :
                krho=(p^2+p+4)/3;
            else:
                krho=(p^2+3*p+2)/3;
        elif vd == 6 :
            krho=(p^2+3)/2;
        elif vd == 8 :
            if p % 3 == 1 :
                krho=(p^2+4*p+1)/3;
            else:
                if vc4 == 3 :
                    krho=(p^2+3*p+2)/3;
                else:
                    krho=(p^2+5)/3;
        elif vd == 9 :
            if p % 4 == 1 :
                 krho=(p^2+6*p+1)/4;
            else:
                if vc6 == 5 :
                    krho=(p^2+4*p+3)/4;
                else:
                    krho=(p^2+7)/4;
        else:
            if p == 5 :
                krho=2;
            else:
                if p % 3 == 2 and not p == 5 :
                    if vc4 == 4 :
                        krho=(p^2+6*p+5)/6;
                    else:
                        krho=(p^2+11)/6;
                else:
                    krho=(p^2+10*p+1)/6;
    return krho;


def conductor(p, EC):
    # This is an auxiliary function that returns the conductor for % p representations associated to EC. Reference: Proposition part II.A of Kraus' paper. #
    # p is a prime, EC an elliptic curve over QQ

    E=EC.minimal_model();
    D=E.discriminant();
    N=E.conductor();
    FN=N.factor();
    cd=1;
    print(N, FN);
    for l in [i[0] for i in FN] :
        print(l);
        if not l == p :
            fl=valuation(N,l);
            print(fl);
            flp=0;
            if E.has_multiplicative_reduction(l):   
                if valuation(D,l) % p == 0 :
                    flp=1;
            else:
                if p > 3 :
                    flp=0;
                    print(flp); 
                elif p == 3 :
                    if E.kodaira_symbol(3) in [KodairaSymbol(i) for i in ["IV", "IV*"]] :
                        flp=1;
                    else:
                        flp=0;
                    print(flp);     
                else:
                    if E.kodaira_symbol(l) == KodairaSymbol("In") :
                        if valuation(D, l) % 2 == 0 and valuation(D, l) > 0 :
                            flp=2;
                        else:
                            flp=1;
                    elif E.kodaira_symbol(l) in [KodairaSymbol(i) for i in ["III", "III*"]] :
                        flp=1;
                    elif E.kodaira_symbol(l) in [KodairaSymbol(i) for i in ["II","II*","IV", "IV*"]] :
                        flp=0;
                    print(flp);     
            print(flp);            
            cd=cd*l^(fl-flp);
            print(cd);
    return cd;



def serre_invariants(EC,p):
    # This function returns the Serre's weight and Serre's conductor for the % p representations associated to EC
    # Input: EC an elliptic curve over QQ, p is a prime
    if not EC.base_field().degree() == 1 :
        print("only implemented for elliptic curves over the rationals");
    else:
        E=EC.minimal_model();
        D=E.discriminant();
        j=E.j_invariant();
        c4=E.c4();
        c6=E.c6(); 
	# Remark before part I.A#
        if E.has_good_reduction(p): 
            krho=2;
        elif E.has_multiplicative_reduction(p): 
            if valuation(j,p) > 0 :
                krho=2;
            else:
                if p == 2 :
                    krho=4;
                else:
                    krho=p+1;
        else:    
            if p == 3 :
                krho=case3(D,j,c4,c6);              
            elif p == 2 :
                krho=case2(D);              
            else:
                krho=casep(p,valuation(j,p),valuation(D,p),valuation(c4,p),valuation(c6,p));

        Nrho=conductor(p, E);
    return krho, Nrho;
