ctfm_expr :=  exp(x) - exp(lv) + 1/(1-x/rv) - 1/(1-lv/rv);
ctfmp_expr := ((diff(ctfm_expr,x)));
octfmp_expr := -(diff(1/ctfmp_expr,x));

ss:=eval(diff(ctfm_expr,x),{lv=-12,rv=3});
evalf(eval(ss,x=2))/evalf(eval(ss,x=-12));
tval_ctfm:=eval((ctfm_expr),{lv=-12,rv=3});
(evalf(eval(tval_ctfm,x=2)));
read ("/home/arman/FD/FD.mpl"):
MFD(); grid_functions:={};
Gen_Eval_Code(ctfm_expr,input="c",proc_name="init_ctfm");
Gen_Eval_Code(ctfmp_expr,input="c",proc_name="init_ctfmp");
Gen_Eval_Code(octfmp_expr,input="c",proc_name="init_octfmp");

grid_functions:={};
ctfm_expr :=  x/(1-x/rv);
ctfmp_expr := simplify(diff(ctfm_expr,x));
octfmp_expr := -simplify(diff(1/ctfmp_expr,x));
Gen_Eval_Code(ctfm_expr,input="c",proc_name="init_ctfm2");
Gen_Eval_Code(ctfmp_expr,input="c",proc_name="init_ctfmp2");
Gen_Eval_Code(octfmp_expr,input="c",proc_name="init_octfmp2");

