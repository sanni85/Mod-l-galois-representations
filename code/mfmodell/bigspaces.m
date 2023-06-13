load "mf.m";

function MakeSpaceLabel(label)
  t := Split(label,".");
  if #t eq 4 then
    N,k,c,i := Explode(t);
    return &cat [N, ".", k, ".", c];
  else
    return label;
  end if;
end function;

DATA_DIR := "/scratch/home/jcremona/Mod-l-galois-representations/data/mfmodell/0/"; 


procedure OneSpace(label: Coeffs:=1000, DegreeBound:=0, Detail:=0)
// e.g. label := "89.4.b.a" or "89.4.b"
  t := Split(label,".");
  if #t eq 4 then
    N,k,c,i := Explode(t);
    space_label := &cat [N, ".", k, ".", c];
  else
    N,k,c := Explode(t);
    space_label := label;
  end if;
  chi := DirichletCharacter(&cat [N, ".", c]);
  weight := StringToInteger(k);
  print "Computing space at level ",N,"; character ",chi,"; weight ",weight;
  s := NewspaceData(chi, weight, 0:
		    ComputeEigenvalues:=true,
		    NumberOfCoefficients:=Coeffs,
		    ReturnDecomposition:=true,
		    DegreeBound:= DegreeBound,
		    Detail:=Detail);
  SetColumns(0);
  outfilename := DATA_DIR cat space_label;
  print "--------------------------------------------";
  print "...done, outputting to ", outfilename;
  PrintFile(outfilename,s: Overwrite:=true);
  print "...done";
end procedure;

procedure DoAll(fname :  spaces_done := [], DegreeBound:=0, Detail:=0)
  s:=Read(fname);
  for dat in Split(s) do
      lab,dim := Explode(Split(dat," "));
      print "############################################";
      print lab, " has dimension ", dim; 
      N,k,c,i := Explode(Split(lab,"."));
      space_label := &cat [N, ".", k, ".", c];
      filename := DATA_DIR cat space_label;
      if ReadTest(filename) then
        print "data file ",filename," already exists, skipping";
        continue;
      end if;
      chi_label := &cat [N, ".", c];
      chi := DirichletCharacter(chi_label);
      space := [chi_label,k];
      if space in spaces_done then
        print "Skipping ",lab," as space ",space, " done already";
        continue;
      end if;
      print "Calling OneSpace(", lab, ")";
      OneSpace(lab : Coeffs:=1000,  DegreeBound:= DegreeBound, Detail:=Detail);
      Append(~spaces_done,space);
  end for;
end procedure;
