/* This code is an implementation of the main theorems in "Détermination du poids et du conducteur associés aux représentations des points de p-torsion d'une courbe elliptique" by Alain Kraus in Rozprawy Matematycznetomnrw serii: 364 wydano: 1997. 

How to use it:

load "kraus.m";
a,b:=serre_invariants(EllipticCurve("15a1"),5);

It will return a=6, b=3, where a is the Serre's weight and b is the Serre's conductor of the representation on the 5 torsion.

Author: Samuele Anni
 */



case2:=function(D)
    /* D is the discriminant of the minimal model of the elliptic curve over QQ
       This is an auxiliary function that returns the weight for mod 2 representations associated to the elliptic curve over QQ which has additive reduction at 2
       Reference: Théorème 3 of Kraus' paper. */
    vd:=Valuation(D,2);
    if vd mod 2 eq 0 then
        krho:=2;
    else
        krho:=4;
    end if;
    return krho;    
end function;


case3:=function(D,j,c4,c6)
    /* D is the discriminant, j the j-invariant, c4 and c6 the C4 and C6 invariants of the minimal model of the elliptic curve over QQ
       This is an auxiliary function that returns the weight for mod 3 representations associated to the elliptic curve over QQ which has additive reduction at 3
       Reference: Théorème 2 of Kraus' paper. */
    vd:=Valuation(D,3);
    D1:=Integers()!(D*3^(-vd));
    if Valuation(j,3) lt 0 then
        if vd mod 3 eq 0 then
            if (D1 mod 9 eq 1) or (D1 mod 9 eq 8) then
                krho:=2;
            else
                krho:=6;
            end if; 
        else
            krho:=8;
        end if;
    else
        if vd eq 3 then
            if Valuation(c6,3) eq 3 then
                if (D1 mod 9 eq 1) or (D1 mod 9 eq 8) then
                    krho:=2;
                else
                    krho:=6;
                end if;
            else
                krho:=6;   
            end if;
        elif vd eq 4 then 
            krho:=8;
        elif vd eq 5 then
            if Valuation(c6,3) eq 3 then
                krho:=8;
            else
                krho:=4;
            end if;
        elif vd eq 6 then
            if Valuation(c6,3) eq 3 then
                if (D1 mod 9 eq 1) or (D1 mod 9 eq 8) then
                    krho:=2;
                else
                    krho:=6;
                end if;
            else
                krho:=6;   
            end if;
        elif vd eq 7 then
            krho:=8;   
        elif vd eq 9 then
            krho:=2; 
        elif vd eq 10 then 
            krho:=4;
        elif vd eq 11 then 
            if Valuation(c6,3) eq 6 then
                krho:=4;
            else
                krho:=8;
            end if;
        elif vd eq 12 then
            krho:=2;
        else
            krho:=4;
        end if;
    end if;
    return krho;    
end function;


casep:=function(p,vj,vd,vc4,vc6)
    /* p is a prime p>3, vj, vd, vc4, vc6 are the p-adic valuations of the j-invariant, the  discriminant, the C4 and C6 invariants of the minimal model of the elliptic curve over QQ
       This is an auxiliary function that returns the weight for mod p representations associated to the elliptic curve over QQ which has additive reduction at p
       Reference: Théorème 1 of Kraus' paper. */
    if vj lt 0 then
        if vj mod p eq 0 then
            krho:=(p^2+3)/2;
        else   
            krho:=(p+1)^2/2; 
        end if;
    else
        if vd eq 2 then
            if p eq 7 then
                if vc4 eq 1 then
                    krho:=14;
                else
                    krho:=2;
                end if;
            elif p mod 3 eq 1 and not p eq 7 then
                krho:=(p^2+4*p+7)/6;
            else
                krho:=(p^2+6*p+5)/6;
            end if;
        elif vd eq 3 then
            if p eq 5 then
                if vc6 eq 2 then
                    krho:=10;
                else
                    krho:=2;
                end if;
            else
                if p mod 4 eq 1 then
                    krho:=(p^2+2*p+5)/4;
                else
                    krho:=(p^2+4*p+3)/4;
                end if;
            end if;
        elif vd eq 4 then
            if p mod 3 eq 1 then
                krho:=(p^2+p+4)/3;
            else
                krho:=(p^2+3*p+2)/3;
            end if;
        elif vd eq 6 then
            krho:=(p^2+3)/2;
        elif vd eq 8 then
            if p mod 3 eq 1 then
                krho:=(p^2+4*p+1)/3;
            else
                if vc4 eq 3 then
                    krho:=(p^2+3*p+2)/3;
                else
                    krho:=(p^2+5)/3;
                end if;
            end if;
        elif vd eq 9 then
           if p mod 4 eq 1 then
                krho:=(p^2+6*p+1)/4;
            else
                if vc6 eq 5 then
                    krho:=(p^2+4*p+3)/4;
                else
                    krho:=(p^2+7)/4;
                end if;
            end if;
        else
            if p eq 5 then
                krho:=2;
            else
                if p mod 3 eq 2 and not p eq 5 then
                    if vc4 eq 4 then
                        krho:=(p^2+6*p+5)/6;
                    else
                        krho:=(p^2+11)/6;
                    end if;
                else
                    krho:=(p^2+10*p+1)/6;
                end if;
            end if;
        end if;
    end if;
    return krho;
end function;


conductor:=function(p, EC)
    /* p is a prime, EC an elliptic curve over QQ
       This is an auxiliary function that returns the conductor for mod p representations associated to EC
       Reference: Proposition part II.A of Kraus' paper. */
    E:=MinimalModel(EC);
    D:=Integers()!Discriminant(E);
    N:=Integers()!Conductor(E);
    FN:=Factorization(N);
    cd:=1;
    for l in [i[1] : i in FN] do
        if not l eq p then
            fl:=Valuation(N,l);
            flp:=0;
            if ReductionType(E,l)[1] in ["N","S"] then
                if Valuation(D,l) mod p eq 0 then
                    flp:=1;
                end if;
            else
                if p gt 3 then
                    flp:=0;
                elif p eq 3 then
                    if KodairaSymbol(E, 3) in [KodairaSymbol(i) : i in ["IV", "IV*"]] then
                        flp:=1;
                    else
                        flp:=0;
                    end if;
                else
                    if KodairaSymbol(E, l) eq KodairaSymbol("In") then
                        if Valuation(D, l) mod 2 eq 0 and Valuation(D, l) gt 0 then
                            flp:=2;
                        else
                            flp:=1;
                        end if;
                    elif KodairaSymbol(E, l) in [KodairaSymbol(i) : i in ["III", "III*"]] then
                        flp:=1;
                    elif KodairaSymbol(E, l) in [KodairaSymbol(i) : i in ["II","II*","IV", "IV*"]] then
                        flp:=0;
                    end if;
                end if;
            end if;
            cd:=cd*l^(fl-flp);
            //print fl;
            //print flp;
            //print cd;
        end if;
    end for;
    return cd;
end function;  



serre_invariants:=function(EC,p)
    /* EC an elliptic curve over QQ, p is a prime
       This function returns the Serre's weight and Serre's conductor for the mod p representations associated to EC
    */

    if not Degree(BaseRing(EC)) eq 1 then
        print "only implemented for elliptic curves over the rationals";
    else
        E:=MinimalModel(EC);
        D:=Integers()!Discriminant(E);
        j:=jInvariant(E);
        c4:=cInvariants(E)[1];
        c6:=cInvariants(E)[2]; 
	/* Remark before part I.A*/
        if ReductionType(E,p)[1] eq "G" then
            krho:=2;
        elif ReductionType(E,p)[1] in ["N","S"] then
	    if Valuation(j,p) gt 0 then 
	        krho:=2;
	    else
		if p eq 2 then 
		    krho:=4;
		else
		    krho:=p+1;	
		end if;
	    end if;	
        else    
            if p eq 3 then
                krho:=case3(D,j,c4,c6);              
            elif p eq 2 then
                krho:=case2(D);              
            else
                krho:=casep(p,Valuation(j,p),Valuation(D,p),Valuation(c4,p),Valuation(c6,p));
            end if;
        end if;
        Nrho:=conductor(p, E);
    end if;
    return Integers()!krho, Integers()!Nrho;
end function;




