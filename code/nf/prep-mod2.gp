/* Prepare data for entry to the lmfdb 
   Each entry is a list of at the bottom 

Good/bad primes computed up to 100

We don't do full label or num here.  Leave that for sage step

*/

assoc1(li, val)=for(j=1,#li,if(li[j][1]==val, return(li[j][2])));return(-1);

/* Hardwiring finite field 
   and s3dict

*/
resp = 2;
projlabel="GL(2,2)";

{
s3dict = [[1, [[1,0,1], "(x+1)^2", 1, 1, [[1,0],[0,1]]]], 
          [2, [[1,0,1], "(x+1)^2", 2, 2, [[0,1],[1,0]]]],
          [3, [[1,1,1], "x^2+x+1", 3, 3, [[0,1],[1,1]]]]];
}
/* for new scheme */
s3mats = [[1, [1,0,0,1]], [2, [0,1,1,0]],  [3, [0,1,1,1]]];
c3mats = [[1, [1,0,0,1]], [0, [0,1,1,1]], [0, [1,1,1,0]]];

getc3mats(nf, gps)=
{
  my(genlist=[],sgp=[1], modu, expo,ans=vector(#gps),pfrob);
  modu = factor(nf.disc); /* Will assume this is a prime power */
  c3m = c3mats;
  modu = modu[1,1];
  expo = (modu-1)/3;
  for(j=1,#gps,
    pfrob = lift(Mod(gps[j],modu)^expo);
	if(setsearch(sgp, pfrob,1), sgp=Set(concat(sgp, [pfrob]));
	  if(#genlist == 0, 
	    genlist=[gps[j]];
	    c3m[2][1] = pfrob,
		c3m[3][1] = pfrob);
	);
	ans[j] = assoc1(c3m, pfrob);
  );
  return([ans, genlist]);
}

gets3gens(nf, gps)=
{
  my(ans=[],need=[2,3],fdat);
  for(j=1,#gps,
    fdat = apply(z->z.f, idealprimedec(nf, gps[j]));
	fdat = prod(j=1,#fdat, fdat[j]);
	if(setsearch(need, fdat),
	  need = setminus(need, [fdat]);
	  ans = concat(ans, [gps[j]]);
	  if(need==[], return(ans))
	  );
  );
  return(ans);
}

/* topslope for ell=2 for a cubic
   Makes use of special facts about the slope of a ramified quadratic over Q_2
*/
topslope(f)=
{
  my(efdat,tote);
  efdat = ef(f,2);
  tote = prod(j=1,#efdat,efdat[j][1]);
  if(tote==1, return(0));
  if(tote==3, return(1));
  return(pdiscf(f,2));
}

goodord(nf,p)=
{
  my(ipd=idealprimedec(nf,p));
  return(lcm(apply(z->z[4], ipd)));
}

ef(pol,p)= return(vecsort(apply(z->[z.e,z.f], idealprimedec(nfinit([pol,[p]]),p)),[1,2]));

weight(f,p)=
{
  my(efdat,lcme=1,pdf);
  efdat = ef(f,p);
  lcme = lcm(apply(z->z[1],efdat));
  if(lcme%resp != 0, return(0)); /* Could be tame */
  /* Only other choice here is e's: 2, 1, so disc comes from the one wild piece */
  return(pdiscf(f,p));
}

/* ent is [conductor, [list of polynomials]) */
doit(ent)=
{
  my(c, pols,v,pol,gps, bps,nf,ggps,bps1,galinfo,reptype, fc, fcps, fcsf,ts,cfact);
  c=ent[1]; pols=ent[2];
  fc = factor(c);
  fcps = vecsort(fc[,1]~);
  /* cfact = vector(#fcps, h, [fc[h,1],fc[h,2]]);  */
  fcsf = prod(j=1,#fcps, fc[j,2]) == 1;
  v=vector(#pols);
  bps1 = factor(c)[,1]~;
  gps = setminus(Set(primes(25)),Set(bps1));
  gps = vecsort(setminus(gps, [resp]));
  bps = apply(z->[z,"","","",""], bps1);
  for(j=1,#pols,
	pol = pols[j];
	galinfo = polgalois(pol);
	t = galinfo[3];
	if(t==1, 
	  projlabel="C_3";
	  imindex = 2;
	  imagelabel="2Cn";
	  imorder=3;
	  projtype="cyclic";
	  absirr=0;
	  reptype = "cyclic"
	  , /* else */
	  projlabel="GL(2,2)";
	  projtype="big";
	  imindex = 1;
	  imorder=6;
	  absirr=1;
	  imagelabel="2G";
	  reptype = "big");
	nf=nfinit(pol);
	baselabel=Str(resp "." 2 "." c);
	if(t==2,
	  ggps = vector(#gps,h,[gps[h], assoc1(s3mats, goodord(nf,gps[h]))]);
	  genlist = gets3gens(nf, gps);
	  ,
	  ggps = getc3mats(nf, gps);
	  genlist = ggps[2];
	  ggps = ggps[1];
	  );
	ts = topslope(pol);
	frobmats = ggps;
	ggps = vector(#ggps, h, [gps[h], ggps[h]]);
	/* skipping num and related objects */
    v[j] = ["GL", bps, resp, 1, 2, c, fcps, fcsf, #fcps, 0, "2.1.1", 2, ggps,
	        imindex, imagelabel, imorder, reptype, absirr, 1, 1, t==2, Vecrev(pol),
			baselabel, t==2, Vecrev(pol), projtype, Str(ts), ts*1., genlist,
			frobmats];

/*
	galinfo[1], reptype, weight(pol,resp), bps1, resp, 1, 2, "big", ggps,
	      Str(pol), [resp, 1], projtype, bps, c, 1, Str(pol), absirr, "solvable", 2,
		  "1.1.1.1", baselabel, projlabel, imagelabel];
		  */
  );
  return(v);
}

/*
   ['algebraic_group', 'bad_prime_list', 'base_ring_characteristic',
    'base_ring_is_field', 'base_ring_order', 'conductor',
    'conductor_num_primes', 'conductor_primes', 'conductor_squarefree',
 'cyclotomic_exponent', 'determinant_label', 'dimension', 'good_prime_list',
 'image_index', 'image_label', 'image_order', 'image_type', 'is_absolutely_irreducible',
 'is_irreducible', 'is_solvable', 'is_surjective', 'kernel_polynomial',
 'label', 'numberfield_label', 'projective_is_surjective',
 'projective_kernel_polynomial', 'projective_type', 'frobenius_matrices']


   ['image_order', 'rep_type', 'weight', 'primes_conductor', 'field_char', 
   'degree_proj_field', 'field_order', 'image_type', 'good_prime_list', 
   'poly_ker', 'field', 'projective_type', 'bad_prime_list', 
   'conductor', 'field_deg', 'poly_proj_ker', 'abs_irr', 'image_at', 'dim', 
   'base_field', 'base_label', 'projective_label', 'image_label']
   */
