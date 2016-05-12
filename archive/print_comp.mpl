print_comp:=proc(A,coord::list,dim::integer,{diag:=false,tp:=dn})

  local i, j;
  local s;
  local cn,cm;

  if dim = 0 then
  s:= "";
  s:= cat(s,convert(A,string));
  s:= cat(s,"_expr");
  s:= cat(s,sprintf(":="));
  s:= cat(s,convert(grcomponent(A),string));
  s:= cat(s,sprintf(";\n"));
  return s;
  end if;
  

  if dim = 1 then
  s:="";
   for i from 1 to nops(coord)-1 do
    cn:=coord[i];
    s:= cat(s,convert(A||cn,string));
    s:= cat(s,"_expr");
    s:= cat(s,sprintf(":="));
    s:=cat(s,convert(grcomponent(A(tp),[coord[i]]),string));
    s:= cat(s,sprintf(";\n"));
   end do;
    return s;
  end if;

  if dim = 2 then
  s:="";
  if not(diag) then 
   for i from 1 to nops(coord)-1 do
   for j from 1 to nops(coord)-1 do
    cn:=coord[i];
    cm := coord[j];
    s:=cat(s,convert(A||cn||cm,string));
    s:= cat(s,"_expr");
    s:=cat(s,sprintf(":="));
    s:=cat(s,convert(grcomponent(A(tp,tp),[coord[i],coord[j]]),string));
    s:=cat(s,sprintf(";\n"));
   end do;
   end do;

  else
   for i from 1 to nops(coord)-1 do
    cn:=coord[i];
    s:=cat(s,convert(A||cn||cn,string));
    s:= cat(s,"_expr");
    s:=cat(s,sprintf(":="));
    s:=cat(s,convert(grcomponent(A(tp,tp),[coord[i],coord[i]]),string));
    s:=cat(s,sprintf(";\n"));
   end do;
   end if;

   return s;
  end if;

end proc;


