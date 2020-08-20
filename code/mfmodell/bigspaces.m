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

procedure DoAll(fname :  spaces_done := [])
  s:=Read(fname);
  for dat in Split(s) do
      lab,dim := Explode(Split(dat," "));
      print lab, " has dimension ", dim; 
      N,k,c,i := Explode(Split(lab,"."));
      chi_label := &cat [N, ".", c];
      chi := DirichletCharacter(chi_label);
      space := [chi_label,k];
      if space in spaces_done then
        print "Skipping ",lab," as space ",space, " done already";
      else
	OneSpace(lab : Coeffs:=1000);
      end if;
  end for;
end procedure;

sp_done := [["37.d", "5"],  ["51.b", "5"],  ["53.b", "6"], ["56.g", "5"], ["89.b", "4"], ["89.a", "6"], ["91.c", "4"], ["97.b", "4"],
			 ["97.a", "6"] , ["91.e", "4"], ["92.c", "3"] ];
