procedure OneSpace(label: Coeffs:=100)
// e.g. label := "89.4.b.a";
  t := Split(label,".");
  if #t eq 4 then
    N,k,c,i := Explode(t);
  else
    N,k,c := Explode(t);
  end if;
  chi := DirichletCharacter(&cat [N, ".", c]);
  weight := StringToInteger(k);
  print "Computing space at level ",N,"; character ",chi,"; weight ",weight;
  s := NewspaceData(chi, weight, 0:
		    ComputeEigenvalues:=true,
		    NumberOfCoefficients:=Coeffs,
		    ReturnDecomposition:=true);
  SetColumns(0);
  print "...done, outputting to ", label;
  PrintFile(label,s: Overwrite:=true);
  print "...done";
end procedure;

procedure DoAll(fname)
  s:=Read(fname);
  for dat in Split(s) do
      lab,dim := Explode(Split(dat," "));
      print lab, " has dimension ", dim; 
  OneSpace(lab : Coeffs:=1000);
  end for;
end procedure;
  
