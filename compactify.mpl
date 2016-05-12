######################################################
# Performing the compactificatoin and non-uniform
# coordinate transformation on the tensor expressions
######################################################

all_pdes:=[ 
   [DPH,0,dn] , [dg,2,dn] , 
         [DDL,2,dn],
   [DTK,0,dn],
         [RR,2,dn],  [P,2,dn],
   [DLAM,1,up],
       [C1,0,dn],
   [HS,0,dn], 
   [DB,0,dn],
   [HS3,0,dn], 
   [MK,1,up],
   [EQRPIb2,0,dn], 
   [EQIPIb2,0,dn],
   [DRPSI,0,dn],
   [DIPSI,0,dn]
             ];

# Going from x -> ctfm(x) using DEtools
with(DEtools);
addcoords(compact,[t,x],[t,ctfm(x)]);


# See built_eqns.mpl for how 'tt' table is defined

for ii from 1 to nops(all_pdes) do
   nm := all_pdes[ii][1]:
   rank := all_pdes[ii][2]:
   tp := all_pdes[ii][3]:
   printf("Transforming %a\n",nm);
   if rank <> 0 then

    for jj from 1 to nops(coord) do
     if (tt[nm][all][coord[jj]] <> 0 ) then
      tt[nm][allcomp][coord[jj]] := PDEchangecoords(tt[nm][all][coord[jj]],[t,x],compact,[t,x]);
      tt[nm][regcomp][coord[jj]] := PDEchangecoords(tt[nm][reg][coord[jj]],[t,x],compact,[t,x]);

	   # Symbolically replacing derivatives of compactification function 
      # since we do not want to finite difference them)

      tt[nm][allcomp][coord[jj]] := subs({diff(ctfm(x),x)=ctfmp(x),diff(ctfm(x),x,x)=ctfmpp(x)},tt[nm][allcomp][coord[jj]]);
      tt[nm][regcomp][coord[jj]] := subs({diff(ctfm(x),x)=ctfmp(x),diff(ctfm(x),x,x)=ctfmpp(x)},tt[nm][regcomp][coord[jj]]);

      # DEtools will replace x with ctfm(x) symbolically in our expressions, but ctfm(x) is our new x, fixing that
      for kk from 1 to nops(grid_functions) do
        fn := grid_functions[kk];
        tt[nm][allcomp][coord[jj]] := subs(fn(t,ctfm(x))=fn(t,x),tt[nm][allcomp][coord[jj]]);
        tt[nm][regcomp][coord[jj]] := subs(fn(t,ctfm(x))=fn(t,x),tt[nm][regcomp][coord[jj]]);
      end do;

     end if;

    end do;
   else

     # Same as above, symbolically replacing derivatives of compactification function
     if ( tt[nm][all] <> 0 ) then
      tt[nm][allcomp] := PDEchangecoords(tt[nm][all],[t,x],compact,[t,x]);
      tt[nm][regcomp] := PDEchangecoords(tt[nm][reg],[t,x],compact,[t,x]);
      tt[nm][allcomp] := subs({diff(ctfm(x),x)=ctfmp(x),diff(ctfm(x),x,x)=ctfmpp(x)},tt[nm][allcomp]);
      tt[nm][regcomp] := subs({diff(ctfm(x),x)=ctfmp(x),diff(ctfm(x),x,x)=ctfmpp(x)},tt[nm][regcomp]);

      # Same as above, fixing DEtools mess!
      for kk from 1 to nops(grid_functions) do 
       fn := grid_functions[kk];
       tt[nm][allcomp] := subs(fn(t,ctfm(x))=fn(t,x),tt[nm][allcomp]);
       tt[nm][regcomp] := subs(fn(t,ctfm(x))=fn(t,x),tt[nm][regcomp]);
      end do;

     end if;

   end if;
end do:

# Checking if everything is transformed properly at the limit ctfm(x) = x

for ii from 1 to nops(all_pdes) do
   nm := all_pdes[ii][1]:
   rank := all_pdes[ii][2]:
   tp := all_pdes[ii][3]:
   if rank <> 0 then

			 for jj from 1 to nops(coord) do
            if ( tt[nm][all][coord[jj]] <> 0 ) then
				res := simplify(tt[nm][all][coord[jj]] - eval(tt[nm][allcomp][coord[jj]],{ctfm(x) = x, ctfmp(x) = 1, ctfmpp(x) = 0}));
				printf("res = %a\n",res);
				res := simplify(tt[nm][reg][coord[jj]] - eval(tt[nm][regcomp][coord[jj]],{ctfm(x) = x, ctfmp(x) = 1, ctfmpp(x) = 0}));
				printf("res = %a\n",res);
            end if;
			 end do;

   else
      res := simplify(tt[nm][all] - eval(tt[nm][allcomp],{ctfm(x) = x, ctfmp(x) = 1, ctfmpp(x) = 0}));
      printf("res = %a\n",res);
      res := simplify(tt[nm][reg] - eval(tt[nm][regcomp],{ctfm(x) = x, ctfmp(x) = 1, ctfmpp(x) = 0}));
      printf("res = %a\n",res);
   end if;
end do:

##################################################
# Performing a collection operation, i.e. 
# powers of ctfm(x) function, this will help the
# intel c and fortran compiler a lot! 
#################################################

for ii from 1 to nops(all_pdes) do
   nm := all_pdes[ii][1]:
   rank := all_pdes[ii][2]:
   tp := all_pdes[ii][3]:
   if rank <> 0 then

			 for jj from 1 to nops(coord) do
            if (tt[nm][all][coord[jj]] <> 0 ) then
  			   tt[nm][all][coord[jj]] := collect(collect(eval(expand(tt[nm][allcomp][coord[jj]]),{ctfmpp(x)=octfmp(x)*ctfmp(x)^2}),1/ctfmp(x)),1/ctfm(x));
			   tt[nm][reg][coord[jj]] := collect(collect(eval(expand(tt[nm][regcomp][coord[jj]]),{ctfmpp(x)=octfmp(x)*ctfmp(x)^2}),1/ctfmp(x)),1/ctfm(x));
				print(tt[nm][all][coord[jj]]);
            print(tt[nm][reg][coord[jj]]);
            end if;
			 end do;
	 
   else
      tt[nm][all] :=  collect(collect(eval(expand(tt[nm][allcomp]),{ctfmpp(x)=octfmp(x)*ctfmp(x)^2}),1/ctfmp(x)),1/ctfm(x));
      tt[nm][reg] :=  collect(collect(eval(expand(tt[nm][regcomp]),{ctfmpp(x)=octfmp(x)*ctfmp(x)^2}),1/ctfmp(x)),1/ctfm(x));
	   print(tt[nm][all]);
	   print(tt[nm][reg]);
   end if;
end do:

# Adding new members of grid functions to our set:
grid_functions := grid_functions union {ctfm,ctfmp,octfmp};

