resp:=3;
G:=GL(2,3);
sublist := Subgroups(G);
SetColumns(0);

load "mintwin.m";

function replacestring(s,fs,ts)
  if Type(fs) eq SeqEnum then
    for i:=1 to #fs do
      s:=replacestring(s,fs[i],ts[i]);
    end for;
    return s;
  end if;
  s:=CodeToString(255) cat Sprint(s) cat CodeToString(255);
  while Position(s,fs) ne 0 do
    p:=Position(s,fs);
    p:=[p-1,#fs,#s-p-#fs+1];
    s1,s2,s3:=Explode(Partition(Eltseq(s),p));
    s:=&cat(s1) cat ts cat &cat(s3);
  end while;
  return s[[2..#s-1]];
end function;


function tostring(e)
  newent:= Sprintf("[\"%o\", %o, %o, %o, %o, %o, %o,%o,%o,%o,\"%o\", %o, %o, %o,\"%o\",%o,\"%o\",%o,%o,%o,%o,%o,\"%o\",%o, %o,\"%o\",\"%o\",%o,%o,%o,\"%o\",\"%o\"]",
    e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9], e[10], e[11],e[12],e[13], e[14], e[15], e[16],e[17], e[18], e[19], e[20], e[21], e[22], e[23], e[24], e[25],e[26],e[27],e[28],e[29],e[30],e[31], e[32]);
  newent:=replacestring(newent, "<", "[");
  newent:=replacestring(newent, ">", "]");
  newent:=DelSpaces(newent);
  return newent;
	/*
  v := <"GL", bps, resp, 1, 3, c, fcps, fcsf, #fcps, cyclopow, detlabel, 2, ggps,
                imindex, imlabel, imorder, reptype, absirr, 1, 1, issurj, 
				Coefficients(minpol),
                baselabel, issurj, Coefficients(ppol), projtype, ts, 1.*ts, 
				genlist, [z[2] : z in ggps]>;
				*/
end function;

function getgens(g, flist)
  og:=Order(g);
  /* Try for single element */
  for j:=1 to #flist do
    if og eq Order(ncl<g|[flist[j][2]]>) then
	  return [flist[j][1]];
	end if;
  end for;
  /* Try 2 elements */
  for j:=1 to (#flist -1) do
    for k:=j+1 to #flist  do
      if og eq Order(ncl<g|[flist[j][2],flist[k][2]]>) then
	    return [flist[j][1],flist[k][1]];
	  end if;
    end for;
  end for;
  /* Throw an error since we should not get here */
  assert false;
end function;

function makevec(ent)
  pol:=ent[1];
  ts := ent[3];
  ts:=Coefficient(ts,0);
  c:=ent[2];
  c:=Coefficient(c,0);
  K:=NumberField(pol);
  gal:=GaloisGroup(K);
  isisom:=false; j:=0;
  minpol:=mintwin(pol);
  while not isisom do
	j:=j+1;
    isisom,mapp := IsIsomorphic(gal, sublist[j]`subgroup);
  end while;
  t := Degree(pol) lt 8 select -1 /* C_4 */ else
    TransitiveGroupIdentification(gal);
  if t eq 1 then /* C_8 */
      imorder := 8; absirr := 0; isirr := 1; imindex:=6; imlabel:="3Cn";
	  cyclopow:=1;reptype:="cyclic";projtype:="Cn";
	  gal1 := "8.1"; pgal1 := "4.1";
	  ppol:=Polredabs(res(pol,4)[1]);
  elif t eq 4 then /* D_4 */
	  imorder := 8; absirr := 1; isirr := 1; imindex:=6; imlabel:="3Ns";
	  cyclopow:=1;reptype:="dihedral";projtype:="Dn";
	  gal1 := "8.3"; pgal1 := "4.2";
	  ppol:=Polredabs(res(pol,4)[1]);
	  ppol:=Polredabs(res2(pol,4,2)[1]);
  elif t eq 5 then /* Q_8 */
	  imorder := 8; absirr := 1; isirr := 1; imindex:=6; imlabel:="3Nn[2]";
	  cyclopow:=0;reptype:="other";projtype:="Dn";
	  gal1 := "8.4"; pgal1 := "4.2";
	  ppol:=Polredabs(res2(pol,4,2)[1]);
  elif t eq 8 then /* QD_16 */
	  imorder := 16; absirr := 1; isirr := 1; imindex:=3; imlabel:="3Nn";
	  cyclopow:=1;reptype:="other";projtype:="Dn";
	  gal1 := "16.8"; pgal1 := "4.3";
	  ppol:=Polredabs(res2(pol,4,3)[1]);
	  ppol:=mintwin(ppol);
  elif t eq 12 then /* SL(2,3) */
	  imorder := 24; absirr := 1; isirr := 1; imindex:=2; imlabel:="3G[]";
	  cyclopow:=0; reptype:="big"; projtype:="A4";
	  gal1 := "24.3"; pgal1 := "12.3";
	  ppol:=Polredabs(res2(pol,4,4)[1]);
  elif t eq 23 then /* GL(2,3) */
	  imorder := 48; absirr := 1; isirr := 1; imindex:=1; imlabel:="3G";
	  cyclopow:=1; reptype:="big"; projtype:="S4";
	  gal1 := "48.29"; pgal1 := "24.12";
	  ppol:=Polredabs(res2(pol,4,5)[1]);
  else   /* C_4 */
	  imorder := 4; absirr := 0; isirr := 1; imindex:=12; imlabel:="3Cn[2]";
	  cyclopow:=1;reptype:="cyclic";projtype:="Cn";
	  gal1 := "4.1"; pgal1 := "2.1";
	  ppol:=Polredabs(res(pol,2)[1]);
  end if;

  detlabel := cyclopow eq 0 select
    Sprintf("%o.%o.1", resp, c) else Sprintf("%o.%o.%o", resp, c,c-1);
  issurj := t eq 24 select 1 else 0;
  fc := Factorization(Integers()!c);
  fcps := [z[1] : z in fc];
  Sort(~fcps);
  fcsf := &*[z[2] : z in fc] eq 1 select 1 else 0;

  bps1 := Set(fcps) diff {resp};
  gps := Set(PrimesUpTo(100)) diff bps1;
  bps1 := Sort([z : z in bps1]);
  gps := Sort([z : z in gps diff {resp}]);
  bps := <<z,"","","",""> : z in bps1>;
  baselabel:=Sprintf("%o.2.%o", resp, c);

  ggps1:= <<z, mapp(FrobeniusElement(K,z))> : z in gps>;
  genlist:= getgens(sublist[j]`subgroup, ggps1);
  ggps:= <<z[1], Eltseq(z[2])> : z in ggps1>;
  v := <"GL", bps, resp, 1, 3, c, fcps, fcsf, #fcps, cyclopow, detlabel, 2, ggps,
                imindex, imlabel, imorder, reptype, absirr, 1, 1, issurj, 
				Coefficients(minpol),
                baselabel, issurj, Coefficients(ppol), projtype, ts, 1.*ts, 
				genlist, [z[2] : z in ggps], gal1, pgal1>;
  if t eq 23 then
    elts := {z[2] : z in ggps1 | Order(z[2]) eq 8};
	elts := [z : z in elts];
	if #elts ne 2 then
	  elts := [FrobeniusElement(K, p) : p in PrimesUpTo(500)];
      elts := {mapp(z) : z in elts | Order(z) eq 8};
	  elts := [z : z in elts];
	  assert #elts eq 2;
	end if;
	ggps2 := [z : z in ggps1];
	for j:= 1 to #ggps2 do
	  if Order(ggps2[j][2]) eq 8 then
	    ggps2[j][2] := elts[1] eq ggps1[j][2] select elts[2] else elts[1];
	  end if;
	end for;
    v2 := v;
	v2[13] := <<z[1], Eltseq(z[2])> : z in ggps2>;
	v2[30] := [z[2] : z in v2[13]];
	return <v, v2>;
  end if;
  if t eq 1 then /* C_8 has 2 reps */
	cc:=Classes(G);
	cm:=ClassMap(G);
	ggps2 := <<z[1], cc[cm(z[2]^(-1))][3]> : z in ggps1>;
    v2 := v;
	v2[13] := <<z[1], Eltseq(z[2])> : z in ggps2>;
	v2[30] := [z[2] : z in v2[13]];
	return <v, v2>;
  end if;
  return <v>;
end function;

