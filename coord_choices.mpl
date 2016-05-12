######################################
# Implementation of the lapse choices
# in standard G-BSSN
# Gamma drive and 1 + log
# Also some experimental other lapse
# choices
######################################

read "/home/arman/FD/FD.mpl":

CFD();
MFD();

grid_functions:= {alpha,K,beta,Lamx,DLamx,BB,ctfmp};

# A generalization of 1 + log
eq1 := diff(alpha(t,x),t) = -gma*alpha(t,x)^gmb*K(t,x);

AVGT := f -> (FD(f,[[1],[0,0,0]]) + FD(f,[[0],[0,0,0]]))/2 ;
FD_table[t] := [ [0] , [0,1] ];


spec := [ 
  { i = [1,1,1] , b=xmin } = Gen_Sten(lhs(eq1)) - AVGT(Gen_Sten(rhs(eq1))) ,
  { i = [2,Nx-1,1]       } = Gen_Sten(lhs(eq1)) - AVGT(Gen_Sten(rhs(eq1))),
  { i=[Nx,Nx,1], b=xmax  } = alpha(n+1,i) + myzero*x(i) - 1
];


A_Gen_Solve_Code(spec,{alpha(n+1,i)},proc_name="evol_alphaopl");
A_Gen_Res_Code(spec,input="d",proc_name="resid_alphaopl");

# Some kinda generalization of 1 + log
eq2:= diff(alpha(t,x),t)  + epsal*diff(K(t,x),t) =  -epsal*ck*K(t,x);



spec2 := [ 
  { i = [1,1,1] , b=xmin } = Gen_Sten(lhs(eq2)) - AVGT(Gen_Sten(rhs(eq2))) ,
  { i = [2,Nx-1,1]       } = Gen_Sten(lhs(eq2)) - AVGT(Gen_Sten(rhs(eq2))),
  { i=[Nx,Nx,1], b=xmax  } = alpha(n+1,i) + myzero*x(i) - 1
];


A_Gen_Solve_Code(spec2,{alpha(n+1,i)},proc_name="evol_alphakd");
A_Gen_Res_Code(spec2,input="d",proc_name="resid_alphakd");


# This is gamma drive

eq3 := diff(beta(t,x),t) = mus*Lamx(t,x) - eta*beta(t,x);
eq4 := diff(alpha(t,x),t) = -2*alpha(t,x)*K(t,x);

spec := [ 
  { i = [1,1,1] , b=xmin } = Gen_Sten(lhs(eq4)) - AVGT(Gen_Sten(rhs(eq4))) ,
  { i = [2,Nx-1,1]       } = Gen_Sten(lhs(eq4)) - AVGT(Gen_Sten(rhs(eq4))),
  { i=[Nx,Nx,1], b=xmax  } = alpha(n+1,i) + myzero*x(i) - 1
];


A_Gen_Solve_Code(spec,{alpha(n+1,i)},proc_name="evol_alphadyopl");
A_Gen_Res_Code(spec,input="d",proc_name="resid_alphadyopl");


spec := [ 
  { i = [1,1,1] , b=xmin } = Gen_Sten(lhs(eq3)) - AVGT(Gen_Sten(rhs(eq3))) ,
  { i = [2,Nx-1,1]       } = Gen_Sten(lhs(eq3)) - AVGT(Gen_Sten(rhs(eq3))),
  { i=[Nx,Nx,1], b=xmax  } = beta(n+1,i) + myzero*x(i) 
];

A_Gen_Solve_Code(spec,{beta(n+1,i)},proc_name="evol_betagd");
A_Gen_Res_Code(spec,input="d",proc_name="resid_betagd");


# Some strange driver, for testing
eq5 := diff(beta(t,x),t) = mus*diff(K(t,x),x) - eta*beta(t,x);

spec := [ 
  { i = [1,1,1] , b=xmin } = beta(n+1,i) + myzero*x(i) ,
  { i = [2,Nx-1,1]       } = Gen_Sten(lhs(eq5)) - AVGT(Gen_Sten(rhs(eq5))),
  { i=[Nx,Nx,1], b=xmax  } = beta(n+1,i) + myzero*x(i) 
];

A_Gen_Solve_Code(spec,{beta(n+1,i)},proc_name="evol_betakd");
A_Gen_Res_Code(spec,input="d",proc_name="resid_betakd");


# More playing around with shift condition
eq6:= diff(beta(t,x),t)  + epsal*diff(Lamx(t,x),t) =  -epsal*ck*Lamx(t,x);


spec2 := [ 
  { i = [1,1,1] , b=xmin } = Gen_Sten(lhs(eq6)) - AVGT(Gen_Sten(rhs(eq6))) ,
  { i = [2,Nx-1,1]       } = Gen_Sten(lhs(eq6)) - AVGT(Gen_Sten(rhs(eq6))),
  { i=[Nx,Nx,1], b=xmax  } = beta(n+1,i) + myzero*x(i)
];


A_Gen_Solve_Code(spec2,{beta(n+1,i)},proc_name="evol_betagd2");
A_Gen_Res_Code(spec2,input="d",proc_name="resid_betagd2");


eq7 := diff(beta(t,x),t) = advc*beta(t,x)*diff(beta(t,x),x)/ctfmp(x) + mus*BB(t,x); 
eq8 := diff(BB(t,x),t) = advc*beta(t,x)*diff(BB(t,x),x)/ctfmp(x) - eta*BB(t,x) + DLamx(t,x);
eq9 := diff(alpha(t,x),t) = advc*beta(t,x)*diff(alpha(t,x),x)/ctfmp(x) - 2*alpha(t,x)*K(t,x);

spec := [ 
  { i = [1,1,1] , b=xmin } = beta(n+1,i)+myzero*x(i) ,
  { i = [2,Nx-1,1]       } = Gen_Sten(lhs(eq7)) - AVGT(Gen_Sten(rhs(eq7))),
  { i=[Nx,Nx,1], b=xmax  } = beta(n+1,i) + myzero*x(i)
];

A_Gen_Solve_Code(spec,{beta(n+1,i)},proc_name="evol_beta_hgd_adv");
A_Gen_Res_Code(spec,input="d",proc_name="resid_beta_hgd_adv");


spec := [ 
  { i = [1,1,1] , b=xmin } = BB(n+1,i)+myzero*x(i) ,
  { i = [2,Nx-1,1]       } = Gen_Sten(lhs(eq8)) - AVGT(Gen_Sten(rhs(eq8))),
  { i=[Nx,Nx,1], b=xmax  } = BB(n+1,i) + myzero*x(i)
];

A_Gen_Solve_Code(spec,{BB(n+1,i)},proc_name="evol_BB_hgd_adv");
A_Gen_Res_Code(spec,input="d",proc_name="resid_BB_hgd_adv");

spec := [ 
  { i = [1,1,1] , b=xmin } =  Gen_Sten(lhs(eq9)) - AVGT(Gen_Sten(rhs(eq9))),
  { i = [2,Nx-1,1]       } = Gen_Sten(lhs(eq9)) - AVGT(Gen_Sten(rhs(eq9))),
  { i=[Nx,Nx,1], b=xmax  } = alpha(n+1,i) - 1 + myzero*x(i)
];

A_Gen_Solve_Code(spec,{alpha(n+1,i)},proc_name="evol_alpha_hgd_adv");
A_Gen_Res_Code(spec,input="d",proc_name="resid_alpha_hgd_adv");


