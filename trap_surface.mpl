############################################################
# A function evaulator routine that
# finds the formation of trap surface
# i.e. blackhole formation, see the paper
# for its definition: http://arxiv.org/pdf/1508.01614v2.pdf
############################################################
read "/home/arman/FD/FD.mpl":

CFD();
MFD();

grid_functions:= {a1,b1,alpha,beta,ctfm,ctfmp};

AVGT := f -> (FD(f,[[1],[0,0,0]]) + FD(f,[[0],[0,0,0]]))/2 ;
FD_table[t] := [ [0] , [0,1] ];



spec := [ 
  { i = [1,1,1] , b=xmin } = 1+myzero*x(i),
  { i = [2,Nx-1,1]       } = (AVGT(Gen_Sten(1/alpha(t,x)))*Gen_Sten(diff(b1(t,x),t)) + AVGT( Gen_Sten( (1/a1(t,x) - beta(t,x)/alpha(t,x))*( diff(b1(t,x),x)/ctfmp(x) + b1(t,x)/ctfm(x) ) ) ) )/ (AVGT(Gen_Sten(b1(t,x)))),
  { i=[Nx,Nx,1], b=xmax  } = myzero*x(i)
];


A_Gen_Eval_Code(spec,input="d",proc_name="compute_th");
