ctfm_expr :=  exp(x) - exp(lv) + ((1/(1-x/rv) - 1/(1-lv/rv)))^(rv-x+1) ;
ctfmp_expr := ((diff(ctfm_expr,x)));
octfmp_expr := -(diff(1/ctfmp_expr,x));

read ("/home/arman/FD/FD.mpl"):
MFD(); grid_functions:={};
Gen_Eval_Code(ctfm_expr,input="c",proc_name="init_ctfm2");
Gen_Eval_Code(ctfmp_expr,input="c",proc_name="init_ctfmp2");
Gen_Eval_Code(octfmp_expr,input="c",proc_name="init_octfmp2");















