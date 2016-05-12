################################################################
# Utility function for irregular term extraction
# Usage:
#
#    gencoeffs(a/x+c+b*x+a*x^3,x);
# 
#     Returns ->>  [-1, a], [0, c], [1, b], [2, 0], [3, a]
#
#  (i.e. it returns a list of powers of x and their coefficients
################################################################

gencoeffs := proc(expr::algebraic, var::name)
   local i;
   if (expr <> 0 ) then
   return seq([i,frontend(coeff,[expr,var,i])],
       i=frontend(ldegree,[expr,var]) .. frontend(degree,[expr,var]));
   else 
    return [0,0];
   end if;
end proc;


#######################################################
# List of grid functions symbols defined in the problem
# See BSSN_Spher.mw or BSSN_Spher.mpl
#######################################################

grid_functions:={alpha,K,beta,phi,Axx,Athth,A,B,metxx,metthetatheta,DDLxx,DDLthetatheta,rho,TS,Lamx,em4phi,RRxx,RRthetatheta,JSxx,JSthth,Uxx,Pxx,Uthetatheta,Pthetatheta,C5s,C1s,psi,Sx,PI2,PHI2,U,rPI,iPI,rPHI,iPHI,a1,rpsi,ipsi,rPIb2,iPIb2,b1,a2,rpsidot,ipsidot,C5psi,divbeta};


################################################################
# List of tensors, will be used for getting the RHS of equations
# and looping over all of the equations
################################################################

all_tensors:=[ 
      [EMPH,0,dn],
   [DPH,0,dn] , [dg,2,dn] , 
         [met,2,dn], [DDL,2,dn],
         [DB,0,dn],
   [DTK,0,dn],
         [RR,2,dn], [U,2,dn] , [P,2,dn],
   [DAA,2,dn],
   [DLAM,1,up],
       [C1,0,dn],
       [C5,0,dn],
   [HS,0,dn],
   [HS3,0,dn],
   [MK,1,up]
             ];

# List of coordinates used: (x is r in spherical coordinate here)
coord:=[x,theta,phi];

######################################################################
# The coordiante coefficients that are defined 
# in the definition of the generic 3-metric in curvilinear coordinate.
######################################################################

coordcoef:=[1,x^2,x^2*sin(theta)^2];


####################################################################
# Initiating a table for the diff expression of tensors
# tt[A][all][x] will store the x compononent of tensor A
# tt[A][sep][x] will store the x component of tensor A in seperated
# powers of 'r' 
# tt[A][reg][x] will store the regular part of above
# tt[A][dis][x] will become the finite difference approximation of
# the above. All in all, 'tt' is our working table of everything!
####################################################################

tt:=table([]):



#########################################
# All of the following atomic variables
# will be used for indexing in table 'tt' 
# and shall not be assigned anythin
#########################################

unassign('all');
unassign('sep');
unassign('tres');
unassign('reg');
unassign('dis');
unassign('disreg');
unassign('disinf');
unassign('inf');

#######  IMPORTANT NOTE: ###############################
# Never use i,j,k as indexing as when working with FD as
# they are protected by default to be used as indexing 
# of coordinates x, y , z
########################################################


for jj from 1 to nops(all_tensors) do


  # Tensor name  
  nm:= all_tensors[jj][1]:

  # Tensor rank
  rank := all_tensors[jj][2]:

  # Tensor type
  tp  := all_tensors[jj][3]:

  # rank-1 tensors
  if rank = 1 then
  tt[nm][all] := table( [ seq( coord[ii] = expand(grcomponent(nm(tp),[coord[ii]])/coordcoef[ii])  ,ii=1..nops(coord)) ] ):
  end if;

  #Assuming rank-2 tensors are all diagnoal (it is the case in spherical symmetry)
  if rank = 2 then
  tt[nm][all] := table( [ seq( coord[ii] = expand(grcomponent(nm(tp,tp),[coord[ii],coord[ii]])/coordcoef[ii])  ,ii=1..nops(coord))  ] ):
  end if;

  if rank = 0 then
  tt[nm][all] := grcomponent(nm):
  end if;

  #Breaking down the expressions to regular and irregular parts

  if rank <> 0 then
  tt[nm][sep] := table (  [  seq( coord[ii] =  [ gencoeffs(tt[nm][all][coord[ii]],'x')   ]  , ii=1..nops(coord) )   ]  ) ;
  
  end if;

  if rank = 0 then
    tt[nm][sep] := [gencoeffs(tt[nm][all],'x')]; 
  end if;
end do:

# Reading EOM for matter and adding it to tt table
read("eq_evol_matter.mpl");

evol_matters := [ [EQRPIb2,0,dn] , [EQIPIb2,0,dn] , [EQRPSI,0,dn]  ,[EQIPSI,0,dn], [DRPSI,0,dn], [DIPSI,0,dn] ]:

tt[EQRPIb2][all] := rhs(  eq_evol_rPIb2  ):
tt[EQIPIb2][all] := rhs( eq_evol_iPIb2   ):
tt[EQRPSI][all] := rhs(eq_evol_rpsi):
tt[EQIPSI][all] := rhs(eq_evol_ipsi):
tt[DRPSI][all] := diff(rpsi(t,x),x)+myzero*x:
tt[DIPSI][all] := diff(ipsi(t,x),x)+myzero*x:


for jj from 1 to nops(evol_matters) do
  nm:= evol_matters[jj][1]:
  #Breaking down the expressions to regular and irregular parts
  tt[nm][sep] := [gencoeffs(tt[nm][all],'x')]; 
end do:



# Updating all_tensors with newly added tensors:

all_tensors:=[ 
      [EMPH,0,dn],
   [DPH,0,dn] , [dg,2,dn] , 
         [met,2,dn], [DDL,2,dn],
         [DB,0,dn],
   [DTK,0,dn],
         [RR,2,dn], [U,2,dn] , [P,2,dn],
   [DAA,2,dn],
   [DLAM,1,up],
       [C1,0,dn],
       [C5,0,dn],
   [HS,0,dn], 
   [HS3,0,dn], 
   [MK,1,up],
   [EQRPIb2,0,dn], 
   [EQIPIb2,0,dn],
   [EQIPSI,0,dn],
   [EQRPSI,0,dn],
   [DRPSI,0,dn],
   [DIPSI,0,dn]
             ];


#Checking that gencoeffs generated coefficients correctly:
check_res := true;
if (check_res) then

		for jj from 1 to nops(all_tensors) do
		  nm:= all_tensors[jj][1]:
		  rank := all_tensors[jj][2]:
		  tp  := all_tensors[jj][3]:
		  if rank <> 0 then  
			tt[nm][tres] := table ( [ seq (  coord[ii] = simplify(expand( tt[nm][all][coord[ii]] - sum(x^tt[nm][sep][coord[ii]][k][1]*tt[nm][sep][coord[ii]][k][2],k=1..nops(tt[nm][sep][coord[ii]])   ) ))      ,ii=1..nops(coord) ) ] );
		  else 
			tt[nm][tres] := simplify(expand( tt[nm][all] - sum(x^tt[nm][sep][k][1]*tt[nm][sep][k][2],k=1..nops(tt[nm][sep])   ) ))
		  end if;
		end do:



		for jj from 1 to nops(all_tensors) do
		  nm:= all_tensors[jj][1]:
		  rank := all_tensors[jj][2]:
		  tp  := all_tensors[jj][3]:
		  if rank <> 0 then
			 for ii from 1 to nops(coord) do
			  print(tt[nm][tres][coord[ii]]);
			 end do
		  else 
			  print(tt[nm][tres]);
		  end if;
		  
		end do:
end if;


for jj from 1 to nops(all_tensors) do
  nm:= all_tensors[jj][1]:
  rank := all_tensors[jj][2]:
  tp  := all_tensors[jj][3]:
  if rank <> 0 then
    for ii from 1 to nops(coord) do
        if tt[nm][sep][coord[ii]][1][1] >= 0 then
               print(nm + comp + coord[ii] + reg );
        else
                print( nm + comp + coord[ii] + ireg);
                kk := 1;
                while ( tt[nm][sep][coord[ii]][kk][1] < 0 ) do
                   print( tt[nm][sep][coord[ii]][kk][2]*x^tt[nm][sep][coord[ii]][kk][1] );
                    kk := kk + 1;
                end do
        end if;
    end do;
  else 
       if tt[nm][sep][1][1] >= 0 then 
          print(nm + comp + reg);
       else
          print(nm + comp + ireg);
          kk := 1;
          while ( tt[nm][sep][kk][1] < 0 ) do
                   print( tt[nm][sep][kk][2]*x^tt[nm][sep][kk][1] );
                    kk := kk + 1;
          end do
        end if;
  end if;
  
end do:

# End of residual check


#### NOTE: ############################################
# Below was effectively disabled by setting all to 0
# i.e. we decided to leave the expressions as they are
#######################################################

# Table to regularize the negative powers of see the following
reg_tbl:= table( [  
 [dg,theta,-1] = 0 , 
 [dg,phi,-1]   = 0 ,
 [DDL,theta,-1] = 0 ,
 [DDL,phi,-1]  = 0 , 
 [RR,x,-2] = 0 ,
 [RR,x,-1] = 0,
 [RR,theta,-2] = 0,
 [RR,theta,-1] = 0,
 [RR,phi,-2] = 0,
 [RR,phi,-1] = 0,
 [P,theta,-1] = 0,
 [P,phi,-1]   = 0,
 [DLAM,x,-2] = 0,
 [DLAM,x,-1] = 0,
 [C1,-1] = 0, 
 [HS,-1] = 0,
 [HS3,-1] = 0,
 [EQRPIb2,-1] = 0 ,
 [EQIPIb2,-1] = 0,
 [MK,x,-1] = 0,
 [DB,-1] = 0
] );


#################################################
# Building list of regular terms in expressions:
# Irregular term replacement was disabled
#################################################
for ii from 1 to nops(all_tensors) do
  nm:= all_tensors[ii][1]:
  rank := all_tensors[ii][2]:
  tp  := all_tensors[ii][3]:
  if rank <> 0 then 
    for jj from 1 to nops(coord) do
        if tt[nm][sep][coord[jj]][1][1] >= 0 then
           tt[nm][reg][coord[jj]] := tt[nm][all][coord[jj]];
        else
           acm := 0;
           for kk from 1 to nops( tt[nm][sep][coord[jj]] ) do
             ord :=  tt[nm][sep][coord[jj]][kk][1];
             trm :=  tt[nm][sep][coord[jj]][kk][2]; 
             if ord < 0 then
               if not(assigned(reg_tbl[ [nm,coord[jj],ord] ] ) ) then
                  print("uncomplete reg table for:");
                  print(nm + coord[jj] + ord);
                  error("quiting...");
               else
                  acm := acm + reg_tbl[ [nm,coord[jj],ord] ];
               end if;
             else
               acm := acm + trm*x^ord;
             end if;
           end do;
           tt[nm][reg][coord[jj]] := acm;
        end if
     end do;
  else

        if tt[nm][sep][1][1] >= 0 then
           tt[nm][reg] := tt[nm][all];
        else
           acm := 0;
           for kk from 1 to nops( tt[nm][sep] ) do
             ord :=  tt[nm][sep][kk][1];
             trm :=  tt[nm][sep][kk][2]; 
             if ord < 0 then
               if not(assigned(reg_tbl[ [nm,ord] ] ) ) then
                  print("uncomplete reg table for:");
                  print(nm + ord);
                  error("quitting...");
               else
                  acm := acm + reg_tbl[ [nm,ord] ];
               end if;
             else
               acm := acm + trm*x^ord;
             end if;
           end do;
           tt[nm][reg]:= acm;
        end if

  end if;
end do:


# Treating 3-Reimmenn tensor seperately

not_all_tensors:=[
         [RR,2,dn]
             ];

R_rest := table([]);
R_ireg := table([]);

for ii from 1 to nops(not_all_tensors) do
  nm:= not_all_tensors[ii][1]:
  rank := not_all_tensors[ii][2]:
  tp  := not_all_tensors[ii][3]:
  if rank <> 0 then
    for jj from 1 to nops(coord) do
        rr_ireg := 0;
        rr_rest := 0;
        for kk from 1 to nops( tt[nm][sep][coord[jj]] ) do
             ord :=  tt[nm][sep][coord[jj]][kk][1];
             trm :=  tt[nm][sep][coord[jj]][kk][2];
             if ord < 0 then
                  rr_ireg:= rr_ireg + trm*x^ord;
             else
                  rr_rest := rr_rest + trm*x^ord;
             end if;
           end do;
        R_ireg[coord[jj]] := rr_ireg;
        R_rest[coord[jj]] := rr_rest;
     end do;
  end if;
end do:

check := true;

#############################################
# Following will print the difference 
# between all and reg aka the irregular terms
# in the tensor components
#############################################

if check then
 for ii from 1 to nops(all_tensors) do
  nm:= all_tensors[ii][1]:
  rank := all_tensors[ii][2]:
  tp  := all_tensors[ii][3]:
  if rank <> 0 then
     for jj from 1 to nops(coord) do 
       print(nm + coord[jj] );
       ress:=simplify(expand(tt[nm][all][coord[jj]] - tt[nm][reg][coord[jj]]));
       print(expand(ress));
     end do:
  else
      print(nm);
      ress:=simplify(expand(tt[nm][all] - tt[nm][reg]));
      print(expand(ress));
  end if;

 end do; 

end if;
