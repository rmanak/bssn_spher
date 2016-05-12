read "print_comp.mpl":
all_var_set:=[ 
   [DPH,0,dn] , [dg,2,dn] , 
         [met,2,dn], [invmet,2,up], [DDL,2,dn],
   [DTK,0,dn],
         [RR,2,dn], [U,2,dn] , [P,2,dn],
   [DAA,2,dn],
   [DLAM,1,up],
   [C1,0,dn],
   [C5,0,dn],
   [HS,0,dn]
             ];
s:="";
for i from 1 to nops(all_var_set) do
  s:=cat(s,print_comp(all_var_set[i][1],[x,theta,phi],all_var_set[i][2],diag=true,tp=all_var_set[i][3])):
end do:
FP:=fopen("allequations.mpl",WRITE);
writeline(FP,s);
fclose(FP);
