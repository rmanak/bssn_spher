############################################
#  Performing finite differencing  using
#  FD toolkit. See previous script 
#  'compactification.mpl' for how the tt table
#  is generated
############################################

read "/home/arman/FD/FD.mpl":
CFD();
MFD();

interface(warnlevel=0);
with(StringTools):

##############################################################
# By default FD uses second order centered in Gen_Sten for RHS
# which is exactly what we want 
##############################################################

for ii from 1 to nops(all_tensors) do
  nm:= all_tensors[ii][1]:
  rank := all_tensors[ii][2]:
  tp  := all_tensors[ii][3]:

  printf("Discretizing %a\n",nm);
  if rank <> 0 then
      # Note the collect operator on powers of 1/hx, it will help the compiler's O3 flag a lot!
      tt[nm][dis] := table( [ seq(  coord[jj] = collect(Gen_Sten(tt[nm][all][coord[jj]]),1/hx) ,jj=1..2 ) ] );
      tt[nm][disreg] := table( [seq( coord[jj] = collect(Gen_Sten(tt[nm][reg][coord[jj]]),1/hx) ,jj=1..2 ) ] );
  else
      tt[nm][dis] := collect(Gen_Sten(tt[nm][all]),1/hx);
      tt[nm][disreg] := collect(Gen_Sten(tt[nm][reg]),1/hx);
  end if;

end do:

###################################################
# Bounday condition at infinity for the tensors:
# Note that we use 'myzero*x' to indicate '0' grid
# function at infinity, See documentation of FD
# for why this is neccessary
# The overhead of this is almost zero!
###################################################

inf_b_cond:= table( [
EMPH = myzero*x+1,
DPH = myzero*x,
DB  = myzero*x,
dg  = [myzero*x,myzero*x],
met = [myzero*x+1,myzero*x+1],
DDL = [myzero*x,myzero*x],
DTK = myzero*x,
RR  = [myzero*x,myzero*x],
U   = [myzero*x,myzero*x],
P   = [myzero*x,myzero*x],
DAA = [myzero*x,myzero*x],
DLAM = [myzero*x,myzero*x],
C1 = myzero*x,
C5 = myzero*x,
HS = psi(t,x)-1 + myzero*x,
HS3 = psi(t,x)-1 + myzero*x,
EQRPIb2 = myzero*x,
EQIPIb2 = myzero*x,
EQIPSI = myzero*x,
EQRPSI = myzero*x,
DRPSI = myzero*x,
DIPSI = myzero*x
]);

# Discretizing the infinity boundary conditions

for ii from 1 to nops(all_tensors) do
  nm:= all_tensors[ii][1]:
  rank := all_tensors[ii][2]:
  tp  := all_tensors[ii][3]:

  printf("Discretizing %a\n",nm);
  if rank <> 0 then
      tt[nm][disinf] := table( [ seq(  coord[jj] = Gen_Sten(inf_b_cond[nm][jj]) ,jj=1..2 ) ] );
  else
      tt[nm][disinf] := Gen_Sten(inf_b_cond[nm]);
  end if;

end do:


# Updating FD's table to default centered in all coordinates

pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl):



############################################
# A table to hold the LHS of BSSN equations
# in their continuous form
############################################

evol_eq_LHS_tbl := table ( [ 
DPH = diff(phi(t,x),t),
dg = [diff(A(t,x),t),diff(B(t,x),t)],
DTK = diff(K(t,x),t),
DAA = [diff(Axx(t,x),t),diff(Athth(t,x),t)],
DLAM = [diff(Lamx(t,x),t),0],
EQRPIb2 = diff(rPIb2(t,x),t),
EQIPIb2 = diff(iPIb2(t,x),t),
EQRPSI  = diff(rpsi(t,x),t),
EQIPSI  = diff(ipsi(t,x),t)
] );


###############################################
# A table to hold what each tensorial equation
# is actually will be solved for, this will be
# needed by FD to decide what to solve the FDA
# for!
###############################################

solve_for_tbl := table ( [ 
   DPH = phi(n+1,i),
   dg  = table( [ x = A(n+1,i) , theta=B(n+1,i)  ] ),
   DTK = K(n+1,i), 
   DAA = table( [ x = Axx(n+1,i), theta=Athth(n+1,i) ] ),
   DLAM = table( [ x = Lamx(n+1,i), theta = 0 ] ),
   EQRPIb2 = rPIb2(n+1,i),
   EQIPIb2 = iPIb2(n+1,i),
   EQRPSI = rpsi(n+1,i),
   EQIPSI = ipsi(n+1,i)
] );

#################################################
# A table to hold what are the actual discrete
# form of the functions are, this will be needed
# by FD
#################################################

dis_func_tbl := table ( [ 
   DPH = phi(n,i),
   dg  = table( [ x = A(n,i) , theta=B(n,i)  ] ),
   DTK = K(n,i), 
   DAA = table( [ x = Axx(n,i), theta=Athth(n,i) ] ),
   DLAM = table( [ x = Lamx(n,i), theta = 0 ] ),
   EQRPIb2 = rPIb2(n,i),
   EQIPIb2 = iPIb2(n,i),
   EQRPSI = rpsi(n,i),
   EQIPSI = ipsi(n,i)
] );


##########################################################
# The implementaion of the famous Kreiss/Oliger disspation 
# using FD's 'FD' operator, will be used in the RHS of the
# evolution equations to suppress numerical noise
# See FD's documentation on how 'FD' operator works, but
# if you are familiar with RNPL, it should be pretty clear
# what is happening here!
##########################################################

DISS := f-> -zepsdis / (16*ht) * ( 6*FD(f,[[0],[0,0,0]]) 
                                   + FD(f,[[0],[2,0,0]]) 
                                   + FD(f,[[0],[-2,0,0]])
                                 - 4*FD(f,[[0],[1,0,0]]) 
                                 - 4*FD(f,[[0],[-1,0,0]]) 
                                );

# Defining average operator for Crank-Nickelson second order update
AVGT := f -> (FD(f,[[1],[0,0,0]]) + FD(f,[[0],[0,0,0]]))/2 ;

# Advanced time operator
NT := f -> FD(f,[[1],[0,0,0]]);

# 2nd order interpolation at the boundary (used at origin)
QFIT := f-> FD(f,[[0],[0,0,0]]) - 4*FD(f,[[0],[1,0,0]])/3 + FD(f,[[0],[2,0,0]])/3 ;

######################################################################
# Shift operator, will be used for the point next to origin condition
# See Eq. (58) in the paper: http://arxiv.org/pdf/1508.01614v2.pdf
#######################################################################

SB := f -> FD(f,[[0],[-1,0,0]]);

# Updating FD's table to default in all coordinates
pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl):

# The evolutionary type equations:
evol_trm := [[DPH,0],[dg,2],[DTK,0],[DAA,2],[DLAM,1], [EQRPIb2,0] , [EQIPIb2,0], [EQRPSI,0], [EQIPSI,0] ];

# Initiating a table to hold the actual equations (residuals) in discrete form
evol_eq_tbl := table ( [] );

# Updating to forward in time differencing scheme, i.e. Crank-Niecklson when combined with AVGT
FD_table[t] := [ [0] , [0,1] ];

# List of all even/odd functions, will be used by FD for imposing inner-boundary conditions
set_of_even_funcs:={alpha,K,phi,Axx,Athth,A,B,psi,a1,a2,b1,ipsi,rpsi,iPIb2,rPIb2,rPI,iPI,divbeta}:
set_of_odd_funcs:={beta,Lamx,rPHI,iPHI}:

for ii from 1 to nops(evol_trm) do
  nm := evol_trm[ii][1];
  rank := evol_trm[ii][2];
  if rank <> 0 then
    for jj from 1 to 2 do
       cn := coord[jj];
       #######  THIS IS THE EQUATION IN FDA FORM ############
       evol_eq_tbl[nm][dis][cn] := Gen_Sten(evol_eq_LHS_tbl[nm][jj]) - AVGT(tt[nm][dis][cn]);
       evol_eq_tbl[nm][disreg][cn] := Gen_Sten(evol_eq_LHS_tbl[nm][jj]) - AVGT(tt[nm][disreg][cn]);
       # Applying boundary conditions using FD's facilities at the origin
       t_mp_1:= A_FD_Even(tt[nm][disreg][cn],x,set_of_even_funcs,0,"forward");
       t_mp_2:= A_FD_Odd(t_mp_1,x,set_of_odd_funcs,0,"forward");
       # Storing  (in 'disorg') a version of FDA that can be used next to the origin 
       evol_eq_tbl[nm][disorg][cn] := Gen_Sten(evol_eq_LHS_tbl[nm][jj]) - AVGT(t_mp_2);
       # Storing the form of the equation at infinity:
       evol_eq_tbl[nm][disinf][cn] := Gen_Sten(evol_eq_LHS_tbl[nm][jj]) - AVGT(tt[nm][disinf][cn]);
    end do
  else
      # Same as above
      evol_eq_tbl[nm][dis] := Gen_Sten(evol_eq_LHS_tbl[nm]) - AVGT(tt[nm][dis]);
      evol_eq_tbl[nm][disreg] := Gen_Sten(evol_eq_LHS_tbl[nm]) - AVGT(tt[nm][disreg]);
      t_mp_1:= A_FD_Even(tt[nm][disreg],x,set_of_even_funcs,0,"forward");
      t_mp_2:= A_FD_Odd(t_mp_1,x,set_of_odd_funcs,0,"forward");
      evol_eq_tbl[nm][disorg] := Gen_Sten(evol_eq_LHS_tbl[nm]) - AVGT(t_mp_2);
      evol_eq_tbl[nm][disinf] := Gen_Sten(evol_eq_LHS_tbl[nm]) - AVGT(tt[nm][disinf]);
  end if;
end do;

# The evaluation type equations, i.e. 'work variables'
eval_term := [[EMPH,0] , [DB,0] ,[met,2],[DDL,2] ,[RR,2],[U,2], [P,2] ,[C1,0], [C5,0], [DIPSI,0], [DRPSI,0]];

eval_expr_tbl := table ( [] );

for ii from 1 to nops(eval_term) do 
  nm := eval_term[ii][1];
  rank := eval_term[ii][2];
  if rank <> 0 then
    for jj from 1 to 2 do
      cn := coord[jj];
      eval_expr_tbl[nm][dis][cn] := (tt[nm][dis][cn]);
      eval_expr_tbl[nm][disreg][cn] := (tt[nm][disreg][cn]);
      eval_expr_tbl[nm][disinf][cn] := (tt[nm][disinf][cn]);
      # Same as evolution case
      t_mp_1:= A_FD_Even(tt[nm][disreg][cn],x,set_of_even_funcs,0,"forward");
      t_mp_2:= A_FD_Odd(t_mp_1,x,set_of_odd_funcs,0,"forward");
      eval_expr_tbl[nm][disorg][cn] := t_mp_2;
    end do;
  else
      eval_expr_tbl[nm][dis]:= (tt[nm][dis]);
      eval_expr_tbl[nm][disreg]:= (tt[nm][disreg]);
      eval_expr_tbl[nm][disinf]:= (tt[nm][disinf]);
      t_mp_1:= A_FD_Even(tt[nm][disreg],x,set_of_even_funcs,0,"forward");
      t_mp_2:= A_FD_Odd(t_mp_1,x,set_of_odd_funcs,0,"forward");
      eval_expr_tbl[nm][disorg] := t_mp_2;    
  end if;
end do;



# Switching to forward in 'x' (ie 'r') for doing FDA at the origin
pl:=table([ t=[-1,-1],x=[0,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl):

# Again in time, we use CN scheme
FD_table[t] := [ [0] , [0,1] ];

#######################################################
# Table for boundary conditions at the origin
# See the paper: http://arxiv.org/pdf/1508.01614v2.pdf 
#######################################################

org_bdc_tbl := table ( [
   # d phi / dx must be zero at origin
   DPH = NT(Gen_Sten(diff(phi(t,x),x))),
   # conformal metric components are 1 at the origin
   dg  = table( [ x = A(n+1,i) - 1 +myzero*x(i) , theta= B(n+1,i) -1 + myzero*x(i)  ] ),
   # exterinsic curvature has zero spatial derivative at origin:
   DTK = NT(Gen_Sten(diff(K(t,x),x))),
   # (from dg eq above) conformal extrinsic curvature are fixed 0 at the origin
   DAA = table( [ x = Axx(n+1,i) +myzero*x(i), theta= Athth(n+1,i)+myzero*x(i) ] ),
   # Conformal connection is vector thus zero at the origin
   DLAM = table( [ x = Lamx(n+1,i) +myzero*x(i), theta= 0 ] ),
   # We use interpolation for matter fields
   EQRPIb2 = NT(QFIT(rPIb2(n,i))),
   EQIPIb2 = NT(QFIT(iPIb2(n,i))),
   EQRPSI  = NT(QFIT(rpsi(n,i))),
   EQIPSI  = NT(QFIT(ipsi(n,i)))
] );

##################################################################
# Extra boundary conditions for the point next to the origin
# as explained in the paper: http://arxiv.org/pdf/1508.01614v2.pdf
# See eq. 58 in the paper
###################################################################

orgp1_bdc_tbl := table ( [
   # Just using the original version, (no extra condition)
   DPH = evol_eq_tbl[DPH][dis],
   # Using higher order derivative condition (see the paper)
   dg  = table( [ x = NT(SB(Gen_Sten(diff(A(t,x),x)))) ,theta= NT(SB(Gen_Sten(diff(B(t,x),x)))) ] ),
   # No BC here
   DTK = evol_eq_tbl[DTK][dis],
   # Same as dg
   DAA  = table( [ x = NT(SB(Gen_Sten(diff(Axx(t,x),x)))) +hx*myzero ,theta= hx*myzero+NT(SB(Gen_Sten(diff(Athth(t,x),x))))  ] ),
   # Same as dg
   DLAM = table( [ x = NT(SB(Gen_Sten(diff(Lamx(t,x),x,x))) ), theta= 0 ] ),
   # No BC here
   EQRPIb2 = evol_eq_tbl[EQRPIb2][dis],
   EQIPIb2 = evol_eq_tbl[EQIPIb2][dis],
   EQRPSI  = evol_eq_tbl[EQRPSI][dis],
   EQIPSI  = evol_eq_tbl[EQIPSI][dis]

] );


# Adding the newly defined FDA's to the original table:
for ii from 1 to nops(evol_trm) do
  nm := evol_trm[ii][1];
  rank := evol_trm[ii][2];
  if rank <> 0 then
    for jj from 1 to 2 do
       cn := coord[jj];
        evol_eq_tbl[nm][disorg][cn] := org_bdc_tbl[nm][cn];
        evol_eq_tbl[nm][disorgp1][cn] := orgp1_bdc_tbl[nm][cn];
    end do
  else
       evol_eq_tbl[nm][disorg] := org_bdc_tbl[nm];
       evol_eq_tbl[nm][disorgp1] := orgp1_bdc_tbl[nm];
  end if;
end do:

# BC for work variables, in reality the values will never be used
e_org_bdc_tbl := table ( [
EMPH = myzero*x+1,
met = [myzero*x+1,myzero*x+1],
DDL = [myzero*x,myzero*x],
RR  = [myzero*x,myzero*x],
U   = [myzero*x,myzero*x],
P   = [myzero*x,myzero*x],
C1 = myzero*x,
C5 = myzero*x,
DRPSI = myzero*x,
DIPSI = myzero*x,
DB  = myzero*x
 ] );

pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl):

# Adding them to the equation table:
for ii from 1 to nops(eval_term) do
  nm := eval_term[ii][1];
  rank := eval_term[ii][2];
  if rank <> 0 then
    for jj from 1 to 2 do
       cn := coord[jj];
        eval_expr_tbl[nm][disorg][cn] := Gen_Sten(e_org_bdc_tbl[nm][jj]);
    end do
  else
        eval_expr_tbl[nm][disorg] := Gen_Sten(e_org_bdc_tbl[nm]);
  end if;
end do:


# Treating Hamiltonian constraint seperately
pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl):


nm := HS;

const_eq_tbl := table ( [] );

const_eq_tbl[nm][dis] := tt[nm][dis];
const_eq_tbl[nm][disreg] := tt[nm][disreg];
t_mp_1 := A_FD_Even(tt[nm][disreg],x,set_of_even_funcs,0,"forward");
t_mp_2 := A_FD_Odd(t_mp_1,x,set_of_odd_funcs,0,"forward");
const_eq_tbl[nm][disorg] := t_mp_2;

pl:=table([ t=[-1,-1],x=[-1,0],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl):

# BC for Hamiltonian constraint at infinity:
const_eq_tbl[nm][inf] := Gen_Sten(ctfm(x)*diff(psi(t,x),x)/ctfmp(x) + psi(t,x) - 1);

pl:=table([ t=[-1,-1],x=[0,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl):

const_eq_tbl[nm][disorg] := Gen_Sten(diff(psi(t,x),x));

pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl):


# Will store header files, to pre-generate the C routine
Other_Headers := NULL;

# DDS for HC is easy:

spec := [   
  { i = [1,1,1],  b=xmin } =  const_eq_tbl[nm][disorg] ,
  { i = [2,Nx-1,1] } = const_eq_tbl[nm][dis] , 
  { i = [Nx,Nx,1], b=xmax  } = const_eq_tbl[nm][inf]
 ];

# Residual evaluator for HC
A_Gen_Res_Code(spec,input="d",proc_name="resid_hs");
Other_Headers := Other_Headers , "resid_hs";

# Solver for HC
A_Gen_Solve_Code(spec,{psi(n,i)},proc_name="solve_hs");
Other_Headers := Other_Headers , "solve_hs";

# Now working with momentum constraint
pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl):

nm := MK;

const_eq_tbl := table ( [] );

const_eq_tbl[nm][dis] := tt[nm][dis][x];
const_eq_tbl[nm][disreg] := tt[nm][disreg][x];
t_mp_1 := A_FD_Even(tt[nm][disreg][x],x,set_of_even_funcs,0,"forward");
t_mp_2 := A_FD_Odd(t_mp_1,x,set_of_odd_funcs,0,"forward");
const_eq_tbl[nm][disorg] := t_mp_2;

spec := [   
  { i = [1,1,1],  b=xmin } =  const_eq_tbl[nm][disorg] ,
  { i = [2,Nx-1,1] } = const_eq_tbl[nm][dis] , 
  { i = [Nx,Nx,1], b=xmax  } = myzero*x(i)
 ];

# Evaluator routine for momentum constraint, for testing purposes
A_Gen_Eval_Code(spec,input="d",proc_name="ire_mk");


##########################################################
# See the bssn_spher.mw for what HS3 is, another form of
# Hamiltonian constraint
##########################################################
pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl):


nm := HS3;

const_eq_tbl := table ( [] );

const_eq_tbl[nm][dis] := tt[nm][dis];

pl:=table([ t=[-1,-1],x=[-1,0],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl);

const_eq_tbl[nm][inf] := Gen_Sten(psi(t,x) - 1 +myzero*ctfm(x));

pl:=table([ t=[-1,-1],x=[0,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl):

const_eq_tbl[nm][disorg] := Gen_Sten(diff(psi(t,x),x)):


pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl):


Other_Headers := NULL;

spec := [   
  { i = [1,1,1],  b=xmin } =  const_eq_tbl[nm][disorg] ,
  { i = [2,Nx-1,1] } = const_eq_tbl[nm][dis] , 
  { i = [Nx,Nx,1], b=xmax  } = const_eq_tbl[nm][inf]
 ]:

A_Gen_Res_Code(spec,input="d",proc_name="resid_hs3");
Other_Headers := Other_Headers , "resid_HS3";

A_Gen_Solve_Code(spec,{psi(n,i)},proc_name="solve_hs3");
Other_Headers := Other_Headers , "solve_HS3";


# End of constraint equation solvers / evaluators


# Initializer expressions for grid functions:

init_alpha_expr:=myzero*x+1;
init_beta_expr:=myzero*x;
init_A_expr:=myzero*x+1;
init_B_expr:=myzero*x+1;
init_Axx_expr:=myzero*x;
init_Athth_expr:=myzero*x;
init_TK_expr:=myzero*x;
init_Lamx_expr:=myzero*x;

# Some redefintion of variables:

phi_expr:=log(psi(t,x));
psi_expr := exp(phi(t,x));
a2_expr:= psi(t,x)^4*A(t,x);
a1_expr := sqrt(a2(t,x));
b1_expr := psi(t,x)^2*sqrt(B(t,x));

lv := 'lv'; 
rv := 'rv';

# Initializer for matter field (Gaussian pulse on real coordinate)
rpsi_expr :=  amp*exp( - ( ctfm(x) - x0 )^2 / deltx^2 );
ipsi_expr := myzero*x*x;


# Note the correction due to compactification
init_rPHI_expr := diff(rpsi_expr, x)/ctfmp(x);
init_iPHI_expr := diff(ipsi_expr, x)/ctfmp(x);

rpsidot_expr := myzero*x;
ipsidot_expr := omega*rpsi_expr;

# PI = a/alpha(psidot - beta*psiprime)

init_rPI_expr := a1(t,x)/alpha(t,x)*(rpsidot(t,x) - beta(t,x)*rPHI(t,x));
init_iPI_expr := a1(t,x)/alpha(t,x)*(ipsidot(t,x) - beta(t,x)*iPHI(t,x));

init_rPIb2_expr := rPI(t,x)*b1(t,x)^2;
init_iPIb2_expr := iPI(t,x)*b1(t,x)^2;

rPI_comp_expr :=  rPIb2(t,x)/b1(t,x)^2;
iPI_comp_expr :=  iPIb2(t,x)/b1(t,x)^2;

PI2_expr := (rPI(t,x)^2 + iPI(t,x)^2);
PHI2_expr := (rPHI(t,x)^2 + iPHI(t,x)^2);

# Energy-momentum tensor expressions:

U_expr := mass*(rpsi(t,x)^2 + ipsi(t,x)^2);

rho_expr:= ( PI2(t,x) + PHI2(t,x) ) / ( 2*a2(t,x)) + U(t,x)/2;

JSxx_expr:=  ( PI2(t,x) + PHI2(t,x) ) / (2) - a2(t,x)*U(t,x)/2;

JSthetatheta_expr :=  (  (PI2(t,x) - PHI2(t,x)) / ( 2*a2(t,x) ) - U(t,x)/2  ) * b1(t,x)^2;

TS_expr := ( 3*PI2(t,x) - PHI2(t,x) ) / (2*a2(t,x) ) - 3/2*U(t,x) ;
Sx_expr := - (  rPI(t,x)*rPHI(t,x) + iPI(t,x)*iPHI(t,x)  ) / (a1(t,x))  ;


############################################################
# all_expr will be used for pre-generation of main.c routine
# this is an experimental feature in FD
############################################################
all_expr := rPHI_expr+mass+omega;

# Creating initializer routines

Gen_Eval_Code(init_alpha_expr,input="c",proc_name="init_alpha");
Gen_Eval_Code(init_beta_expr,input="c",proc_name="init_beta");


Gen_Eval_Code(init_A_expr,input="c",proc_name="init_a");
Gen_Eval_Code(init_B_expr,input="c",proc_name="init_b");

Gen_Eval_Code(init_TK_expr,input="c",proc_name="init_tk");
Gen_Eval_Code(init_Axx_expr,input="c",proc_name="init_axx");
Gen_Eval_Code(init_Athth_expr,input="c",proc_name="init_athth");
Gen_Eval_Code(init_Lamx_expr,input="c",proc_name="init_lamx");


Gen_Eval_Code(rpsi_expr,input="c",proc_name="init_rpsi");
Gen_Eval_Code(ipsi_expr,input="c",proc_name="init_ipsi");

Gen_Eval_Code(init_iPHI_expr,input="c",proc_name="init_iphi");
Gen_Eval_Code(init_rPHI_expr,input="c",proc_name="init_rphi");

Gen_Eval_Code(rpsidot_expr,input="c",proc_name="init_rpsidot");
Gen_Eval_Code(ipsidot_expr,input="c",proc_name="init_ipsidot");

Gen_Eval_Code(init_rPI_expr,input="c",proc_name="init_rpi");
Gen_Eval_Code(init_iPI_expr,input="c",proc_name="init_ipi");

Gen_Eval_Code(init_iPIb2_expr,input="c",proc_name="init_ipib2");
Gen_Eval_Code(init_rPIb2_expr,input="c",proc_name="init_rpib2");



Gen_Eval_Code(phi_expr,input="c",proc_name="compute_phi");
Gen_Eval_Code(psi_expr,input="c",proc_name="compute_psi");
Gen_Eval_Code(a2_expr,input="c",proc_name="compute_a2");
Gen_Eval_Code(a1_expr,input="c",proc_name="compute_a1");
Gen_Eval_Code(b1_expr,input="c",proc_name="compute_b1");

Gen_Eval_Code(iPI_comp_expr,input="c",proc_name="compute_ipi");
Gen_Eval_Code(rPI_comp_expr,input="c",proc_name="compute_rpi");


Gen_Eval_Code(PHI2_expr,input="c",proc_name="compute_phi2");
Gen_Eval_Code(PI2_expr,input="c",proc_name="compute_pi2");

Gen_Eval_Code(U_expr,input="c",proc_name="compute_u");

Gen_Eval_Code(rho_expr,input="c",proc_name="compute_rho");
Gen_Eval_Code(JSxx_expr,input="c",proc_name="compute_jsxx");
Gen_Eval_Code(JSthetatheta_expr,input="c",proc_name="compute_jsthth");
Gen_Eval_Code(TS_expr,input="c",proc_name="compute_ts");
Gen_Eval_Code(Sx_expr,input="c",proc_name="compute_sx");

C5psi_expr := C5s(t,x)*psi(t,x)^5;
Gen_Eval_Code(C5psi_expr,input="c",proc_name="compute_c5psi");


# Newly generated header files, will be passed to C code generator
Other_Headers := Other_Headers , "init_alpha", "init_beta", "init_A" , "init_B", "init_TK", "init_Axx", "init_Athth", 
                "init_Lamx", "compute_phi", "compute_psi", "compute_a2", "compute_a1",
                "compute_b1", "init_rpsi", "init_ipsi", "init_iphi", "init_rphi",
                "init_ipib2", "init_rpib2", "compute_ipi", "compute_rpi", 
                "compute_phi2", "compute_pi2", "compute_u", "compute_rho", "compute_jsxx",
                "compute_jsthth", "compute_ts", "compute_sx", "init_rpsidot", "init_ipsidot", "init_rpi","init_ipi";
 
# Creating DDS for work variables and calculator routines
all_work_tensors:=[ 
      [EMPH,0,dn], [DB,0,dn],
         [met,2,dn], [DDL,2,dn],
         [RR,2,dn], [U,2,dn] , [P,2,dn],
         [C1,0,dn], [C5,0,dn], [DRPSI,0,dn], [DIPSI,0,dn]
             ];

for ii from 1 to nops(all_work_tensors) do
  nm:= all_work_tensors[ii][1]:
  rank := all_work_tensors[ii][2]:
  if rank <> 0 then
     for jj from 1 to 2 do
       cn := coord[jj];
       # DDS for work variables
       spec := [
                     { i = [1,1,1] , b=xmin }     = eval_expr_tbl[nm][disorg][cn]    ,
                     { i = [2,Nx-1,1] }  = eval_expr_tbl[nm][dis][cn]  ,
                     { i = [Nx,Nx,1] , b=xmax }  = eval_expr_tbl[nm][disinf][cn]
               ]; 
      # creating calculator routines fortran routines for work variables
      A_Gen_Eval_Code(spec,input="d",proc_name=LowerCase(cat("compute_",nm,cn)));
      Other_Headers := Other_Headers , LowerCase(cat("compute_",nm,cn));
     end do;     
  else
        spec := [
                     { i = [1,1,1] , b=xmin }     = eval_expr_tbl[nm][disorg]    ,
                     { i = [2,Nx-1,1] }  = eval_expr_tbl[nm][dis]  ,
                     { i = [Nx,Nx,1] , b=xmax  }  = eval_expr_tbl[nm][disinf]
               ]; 
      A_Gen_Eval_Code(spec,input="d",proc_name=LowerCase(cat("compute_",nm))); 
      Other_Headers := Other_Headers , LowerCase(cat("compute_",nm));
  end if;
end do:

# Treating conformal connection equation first:

spec := [
  { i = [1,1,1] , b=xmin } = myzero*x(i),
  { i = [2,Nx-1,1] } = tt[DLAM][dis][x] ,
  { i = [Nx,Nx,1] , b=xmax  } = myzero*x(i)
];

A_Gen_Eval_Code(spec, input="d",proc_name="eval_dlam");

# Creating a table to hold the K/O dissipation at various points of discrete domain
dis_tbl := table([]);
for ii from 1 to nops(evol_trm) do
  nm := evol_trm[ii][1];
  rank := evol_trm[ii][2];
  if rank <> 0 then
    for jj from 1 to 2 do
       cn := coord[jj];
       fn := dis_func_tbl[nm][cn];

       dis_tbl[nm][mid][cn] := -AVGT(DISS(fn));
       dis_tbl[nm][org][cn] := -AVGT(DISS(fn));
       dis_tbl[nm][orgp1][cn] :=  -AVGT(DISS(fn));

       # Using FD's even odd facility to built the correct version of FDA
       # That can be used at the origin / infinity 
       t_mp_1:= A_FD_Even(dis_tbl[nm][org][cn],x,set_of_even_funcs,0,"forward");
       dis_tbl[nm][org][cn]:= A_FD_Odd(t_mp_1,x,set_of_odd_funcs,0,"forward");

       t_mp_1:= A_FD_Even(dis_tbl[nm][orgp1][cn],x,set_of_even_funcs,-1,"forward");
       dis_tbl[nm][orgp1][cn]:= A_FD_Odd(t_mp_1,x,set_of_odd_funcs,-1,"forward");
        
    end do
  else
      fn := dis_func_tbl[nm];

       dis_tbl[nm][mid] := -AVGT(DISS(fn));
       dis_tbl[nm][org] := -AVGT(DISS(fn));
       dis_tbl[nm][orgp1] :=  -AVGT(DISS(fn));

       t_mp_1:= A_FD_Even(dis_tbl[nm][org],x,set_of_even_funcs,0,"forward");
       dis_tbl[nm][org]:= A_FD_Odd(t_mp_1,x,set_of_odd_funcs,0,"forward");

       t_mp_1:= A_FD_Even(dis_tbl[nm][orgp1],x,set_of_even_funcs,-1,"forward");
       dis_tbl[nm][orgp1]:= A_FD_Odd(t_mp_1,x,set_of_odd_funcs,-1,"forward");

  end if;
end do:

####### NOW EVOLUTION EQUATIONS ##########

for ii from 1 to nops(evol_trm) do 
  nm := evol_trm[ii][1]:
  rank := evol_trm[ii][2]:
  if rank <> 0 then
    for jj from 1 to 2 do
      cn := coord[jj];
      fn := dis_func_tbl[nm][cn]:
      # there are 5 conditions, origin, next to origin, middle, next to infinity and infinity 
      # here is the glory of FD! (Also we are adding the right dissipation expression)
      spec := [ 
                 {i = [1,1,1] , b=xmin }    =  evol_eq_tbl[nm][disorg][cn],
                 {i = [2,2,1] , b=xmin }    =  evol_eq_tbl[nm][disorgp1][cn] + dis_tbl[nm][orgp1][cn],
                 {i = [3,Nx-2,1] } =  evol_eq_tbl[nm][dis][cn] + dis_tbl[nm][mid][cn],
                 {i = [Nx-1,Nx-1,1] } =  evol_eq_tbl[nm][dis][cn],
                 {i = [Nx,Nx,1] , b = xmax } =  evol_eq_tbl[nm][disinf][cn]
              ];
      if not( nm=DLAM and cn=theta) then
         print(cat("Solve code for:",nm,cn," Generated."));
         A_Gen_Solve_Code(spec,{solve_for_tbl[nm][cn] },input="d",proc_name=LowerCase(cat("evolve_",nm,cn)));
         Other_Headers := Other_Headers , cat("evolve_",nm,cn);
         A_Gen_Res_Code(spec,input="d",proc_name=LowerCase(cat("resid_",nm,cn)));
         Other_Headers := Other_Headers , cat("resid_",nm,cn);
         print(cat("Res code for:",nm,cn," Generated."));
      end if;
    end do;
  else
      fn := dis_func_tbl[nm];
      spec := [ 
                 {i = [1,1,1] , b=xmin }    =  evol_eq_tbl[nm][disorg],
                 {i = [2,2,1] , b=xmin }    =  evol_eq_tbl[nm][disorgp1] + dis_tbl[nm][orgp1],
                 {i = [3,Nx-2,1] } =  evol_eq_tbl[nm][dis] + dis_tbl[nm][mid],
                 {i = [Nx-1,Nx-1,1] } =  evol_eq_tbl[nm][dis],
                 {i = [Nx,Nx,1] , b=xmax  } =  evol_eq_tbl[nm][disinf]

              ];
    A_Gen_Solve_Code(spec,{solve_for_tbl[nm]},proc_name=LowerCase(cat("evolve_",nm)));
    Other_Headers := Other_Headers , cat("evolve_",nm);
    print(cat("Solve code for:",nm," Generated."));
    A_Gen_Res_Code(spec,input="d",proc_name=LowerCase(cat("resid_",nm)));
    Other_Headers := Other_Headers , cat("resid_",nm);
    print(cat("Res code for:",nm," Generated."));
  end if;
end do:

#Generating driver header:
pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl):

#############################################################
# This is where I use the Other_Headers in FD
# GDC is an experimental feature, it generates a tempalte
# for a main C driver routine.
# It uses all_expr variable to learn all necessary definition
# for a solver routine, enable it and see what it generates
#############################################################

#GDC(all_expr+myzero,input="c",output="main",other_headers={Other_Headers});


# Creating IRE for BSSN equations (notice how quick this is!)

printf("Generating IREs...\n");
pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl):
interface(warnlevel = 0):

for ii from 1 to nops(evol_trm) do
  nm := evol_trm[ii][1];
  rank := evol_trm[ii][2];
  if rank <> 0 then
    for jj from 1 to 2 do
       cn := coord[jj];
       all_expr:=all_expr + ( evol_eq_LHS_tbl[nm][jj] -  tt[nm][all][cn] );
       sten := Gen_Sten( evol_eq_LHS_tbl[nm][jj] -  tt[nm][all][cn]);
       if (sten <> 0 ) then
         Gen_Res_Code(sten,input="d",proc_name=LowerCase(cat("ire_",nm,cn)));
         Gen_Eval_Code(sten,input="d",proc_name=LowerCase(cat("irev_",nm,cn)));
         Other_Headers := Other_Headers , cat("ire_",nm,cn);
       end if;
    end do
  else
       all_expr:=all_expr + ( evol_eq_LHS_tbl[nm] -  tt[nm][all] );
       sten := Gen_Sten( evol_eq_LHS_tbl[nm] -  tt[nm][all]);
       if (sten <> 0) then
         Gen_Res_Code(sten,input="d",proc_name=LowerCase(cat("ire_",nm)));
         Gen_Eval_Code(sten,input="d",proc_name=LowerCase(cat("irev_",nm)));
         Other_Headers := Other_Headers , cat("ire_",nm);
       end if;
  end if;
end do:



