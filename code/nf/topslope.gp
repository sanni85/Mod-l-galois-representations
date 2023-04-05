/* gp program to find the top slope for a polynomial f and a prime p

   Bad news, it uses function and data stored on JJ's hard drive, so
   it won't work for most people */
topslope(f,p)=
{
  my(fp, all, ts=0);
  fp = lift(factorpadic(f,p,500)[,1]~ * Mod(1,p^500));
  for(j=1,#fp,
    all=load(p,poldegree(fp[j]));
    ppol = findinlist(fp[j], all);
    ppol = getpolyent(ppol[1],all);
    if(ppol[11][1]>1 && ts == 0, ts = 1);
    if(#ppol[10]>0 && ppol[10][#ppol[10]] > ts, ts=ppol[10][#ppol[10]]);
  );
  return(ts);
}

