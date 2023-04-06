/* Prepare data for entry to the lmfdb 
   Each entry is a list of at the bottom 

Good/bad primes computed up to 100

We don't do full label or num here.  Leave that for sage step

*/

\r ~/gp_progs/generate.gp

assoc1(li, val)=for(j=1,#li,if(li[j][1]==val, return(li[j][2])));return(-1);
pos(li, val)=for(j=1,#li,if(li[j][1]==val, return(j)));return(-1);

/* Hardwiring finite field 
   and s3dict

*/
resp = 2;
projlabel="GSp";

{
s6mats = [[[1,1,1,1,1,1], [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]],
          [[2,1,1,1,1], [1,1,0,0,0,1,0,0,0,0,1,1,0,0,0,1]],
		  [[2,2,1,1], [1,0,1,1,0,1,0,1,0,0,1,0,0,0,0,1]],
		  [[2,2,2], [1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1]],
		  [[3,1,1,1], [1,1,1,1,1,0,0,1,0,0,1,1,0,0,1,0]],
		  [[3,2,1], [0,0,1,0,1,0,1,0,1,0,1,1,0,1,1,1]],
		  [[3,3], [0,1,1,1,1,0,1,1,1,1,0,1,1,0,0,1]],
		  [[4,1,1], [0,1,0,1,1,0,1,1,0,0,0,1,0,0,1,0]],
		  [[4,2], [0,1,0,0,1,0,0,1,0,0,0,1,0,0,1,0]],
		  [[5,1], [1,1,1,0,1,0,1,1,0,1,1,0,1,0,0,1]],
		  [[6], [1,1,0,1,1,0,0,1,0,1,1,1,1,0,1,1]]];
}
gen1=[2,4,6,8,11];
gen2=Set([[6,7], [6,9]]);

getgenprimes(ggps)=
{
  my(hit=[],revmats, place);
  revmats = vector(#s6mats,h,[s6mats[h][2], s6mats[h][1]]);
  for(j=1,#ggps,
    place = pos(revmats, ggps[j][2]);
	if(setsearch(gen1, place), return([ggps[j][1]]));
	if(pos(hit, place)<0,
	  hit = concat(hit, [[place, ggps[j][1]]]);
	);
	for(k=1,#hit-1, 
	  if(setsearch(gen2, vecsort([place, hit[k][1]]))>0,
	    return(vecsort([ggps[j][1], hit[k][2]]))));
  );
  error("Cannot find generators");
}

/* topslope for ell=2 for a cubic
   Makes use of special facts about the slope of a ramified quadratic over Q_2
*/
topslope(f)=
{
  my(fp, all, ts=0);
  fp = lift(factorpadic(f,2,500)[,1]~ * Mod(1,2^500));
  for(j=1,#fp,
    all=load(2,poldegree(fp[j]));
	ppol = findinlist(fp[j], all);
	ppol = getpolyent(ppol[1],all);
	if(ppol[11][1]>1 && ts == 0, ts = 1);
	if(#ppol[10]>0 && ppol[10][#ppol[10]] > ts, ts=ppol[10][#ppol[10]]);
  );
  return(ts);
}

cycletype(nf,p)=
{
  my(ipd=idealprimedec(nf,p),fs);
  return(vecsort(apply(z->z[4], ipd),,4)); /* sort reverse order */
}

ef(pol,p)= return(vecsort(apply(z->[z.e,z.f], idealprimedec(nfinit([pol,[p]]),p)),[1,2]));


/* ent is [polynomial, conductor, ?) */
doit(ent)=
{
  my(c, minpol,v,pol,gps, bps,nf,ggps,bps1,galinfo,reptype, fc, fcps, fcsf,ts,cfact);
  my(genlist);
  c=ent[2]; pol=ent[1]; minpol=ent[3];
  fc = factor(c);
  fcps = vecsort(fc[,1]~);
/*  cfact = vector(#fcps, h, [fc[h,1],fc[h,2]]); */
  fcsf = prod(j=1,#fcps, fc[j,2]) == 1;
  v=vector(#pols);
  bps1 = factor(c)[,1]~;
  gps = setminus(Set(primes(25)),Set(bps1));
  gps = vecsort(setminus(gps, [resp]));
  bps = apply(z->[z,"","","",""], bps1);
  nf=nfinit(pol);
  baselabel=Str(resp "." 4 "." c);
  /* ******************** */
  ggps = vector(#gps);
  for(j=1,#ggps, ggps[j]=[gps[j], assoc1(s6mats, cycletype(nf, gps[j]))]);
  genlist = getgenprimes(ggps);
  ts = topslope(pol);
  imindex = 1;
  imagelabel = "1.1.0.1";
  imorder = 720;
  absirr = 1;
  projtype = "big";
  reptype = "big";

	/* skipping num and related objects */
  v = ["GSp", bps, resp, 1, 2, c, fcps, fcsf, #fcps, 0, "2.1.1", 4, ggps,
	        imindex, imagelabel, imorder, reptype, absirr, 1, 0, 1, Vecrev(minpol),
			baselabel, 1, Vecrev(minpol), projtype, Str(ts), ts*1., genlist,
			apply(z->z[2], ggps)];

  return(v);
}

/*
   ['algebraic_group', 'bad_prime_list', 'base_ring_characteristic',
    'base_ring_is_field', 'base_ring_order', 'conductor',
    'conductor_primes', 'conductor_num_primes', 'conductor_squarefree',
 'cyclotomic_exponent', 'determinant_label', 'dimension', 'good_prime_list',
 'image_index', 'image_label', 'image_order', 'image_type', 'is_absolutely_irreducible',
 'is_irreducible', 'is_solvable', 'is_surjective', 'kernel_polynomial',
 'label', 'numberfield_label', 'projective_is_surjective',
 'projective_kernel_polynomial', 'projective_type', 'related_objects']


   */
