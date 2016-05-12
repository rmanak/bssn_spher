RG := proc(A::algebraic,var::symbol,n::integer)
   local B1,B2,B3,i;
   local B;
   B:=A;
   for i from 1 to n do
   B1:=simplify(B*var^(n-i+1));
   B2:=simplify(eval(B1,var=0));
   if (B2 <> 0) then
   B3:=eval(B2,0=var); 
   lprint(B3/var^(n-i+1));
   B:=simplify(B-B3/var^(n-i+1));
   end if;
   end do;
   return B;
end proc;
