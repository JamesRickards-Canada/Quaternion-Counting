
\\Given a list of definite and indefinite discriminants, returns the proportion that is ramified at p. Returns [def, indef, total].
prop_ram(~def, ~indef, p)={
  my(cdef=0, cindef=0);
  for(i=1,#def,if(def[i]%p==0,cdef++));
  for(i=1,#indef,if(indef[i]%p==0,cindef++));
  return(1.*[cdef/#def, cindef/#indef, (cdef+cindef)/(#def+#indef)]);
}