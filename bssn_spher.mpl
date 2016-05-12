################################
# Derivation of BSSN equations
# Using GR-Tensor
# See bssn_spher.mw for interactive
# version of this script
################################
grtw();
grload(g,"metrics/spher_3.mpl");
grdisplay(g(dn,dn));
grcalc(g(up,up));
grdisplay(g(up,up));
grdef(`g_f{a b} := kdelta{a $x}*kdelta{b $x}*1 + kdelta{a $theta}*kdelta{b $theta}*x^2 + kdelta{a $phi}*kdelta{b $phi}*x^2*(sin(theta))^2`); grcalc(g_f(dn,dn)); grdisplay(g_f(dn,dn));
grdef(`g_f_inv{^a ^b} := kdelta{^a $x}*kdelta{^b $x}*1 + kdelta{^a $theta}*kdelta{^b $theta}*1/x^2 + kdelta{^a $phi}*kdelta{^b $phi}*1/(x^2*(sin(theta))^2)`); grcalc(g_f_inv(up,up)); grdisplay(g_f_inv(up,up));
grdef(`ChL_f{d a b} := 1/2*( g_f{b d ,a} + g_f{a d,b} - g_f{a b,d})`);
grdef(`MCh_f{^c a b} := g_f_inv{^c ^d} * ChL_f{d a b}`):
grcalc(MCh_f(up,dn,dn)); grdisplay(MCh_f(up,dn,dn));
grdef(`ChL{d a b} := 1/2*( g{b d ,a} + g{a d,b} - g{a b,d})`):
grdef(`MCh{^c a b} := g{^c ^d} * ChL{d a b}`):
grcalc(MCh(up,dn,dn));
gralter(MCh(up,dn,dn),simplify,expand);
grdisplay(MCh(up,dn,dn));
grdef(`TG{^i} := g{^j ^k}*MCh{^i j k}`);
grcalc(TG(up));
grdisplay(TG(up));
grdef(`TG2{^i} := g{^j ^k}*(MCh{^i j k}-MCh_f{^i j k})`); grcalc(TG2(up)); gralter(TG2(up),simplify,expand); grdisplay(TG2(up));
grdef(`Lam{^a} := [Lamx(t,x),0,0]`);
grcalc(Lam(up)); grdisplay(Lam(up));
grdef(`GA{^i} := Lam{^i} + g{^j ^k}*MCh_f{^i j k}`);  grcalc(GA(up)); grdisplay(GA(up));
grdef(`SS{i j} := 2*g{^l ^m}*MCh{^k l (i}*ChL{j) k m} + g{^l ^m}*MCh{^k i m}*ChL{k l j}`);
grcalc(SS(dn,dn));
gralter(SS(dn,dn),simplify);
grdisplay(SS(dn,dn));
grdef(`RB{i j} := -1/2*g{^l ^m}*g{i j ,l ,m} + g{k (i}*GA{^k ,j)} + GA{^k}*ChL{(i j) k} + SS{i j}`);
grcalc(RB(dn,dn));
gralter(RB(dn,dn),simplify,expand);
grdisplay(RB(dn,dn));
grdef(`CF{} := phi(t,x)`);
grcalc(CF);
grdef(`DCF{i} := CF{;i}`);
grcalc(DCF(dn));
grdisplay(DCF(dn));
grdef(`Rphi{i j} := -2*DCF{j ;i} - 2*g{i j}*DCF{^k ;k} + 4*DCF{i}*DCF{j} - 4*g{i j}*DCF{^k}*DCF{k}`);
grcalc(Rphi(dn,dn));
gralter(Rphi(dn,dn),expand,simplify,expand);
grdisplay(Rphi(dn,dn));
grcalc(R(dn,dn)); grdisplay(R(dn,dn));
grdef(`RBtest{i j} := -1/2*g{^l ^m}*g{i j ,l ,m} + g{k (i}*TG{^k ,j)} + TG{^k}*ChL{(i j) k} + SS{i j}`);
grdef(`Rresid{i j} := RBtest{i j} - R{i j}`);
grcalc(Rresid(dn,dn));
gralter(Rresid(dn,dn),simplify);
grdisplay(Rresid(dn,dn));
grdef(`Sh{^a}:=[beta(t,x),0,0]`);
grcalc(Sh(up));
grdisplay(Sh(up));
grdef(`Q{i j} := Sh{^k}*g{i j ,k} + g{i k}*Sh{^k ,j} + g{k j}*Sh{^k ,i} -2/3*g{i j}*vee*divbeta(t,x)`);
grcalc(Q(dn,dn));
grdisplay(Q(dn,dn));
grdef(`A{a b} := kdelta{a $x}*kdelta{b $x}*Axx(t,x)/1 + kdelta{a $theta}*kdelta{b $theta}*x^2*Athth(t,x) + kdelta{a $phi}*kdelta{b $phi}*x^2*(sin(theta))^2*Athth(t,x)`);
grcalc(A(dn,dn));grdisplay(A(dn,dn));
grdef(`P{i j} := Sh{^k}*A{i j,k} + A{i k}*Sh{^k ,j} + A{k j}*Sh{^k ,i} - 2/3*A{i j}*vee*divbeta(t,x)`):
grcalc(P(dn,dn));
grdisplay(P(dn,dn));
grdef(`Extc{} := K(t,x)`);
grcalc(Extc);
grdef(`Lapse{} := alpha(t,x)`);
grcalc(Lapse);
grdef(`S{i} :=[Sx(t,x),0,0]`);
grcalc(S(dn));

grdef(`DBS{} := divbeta(t,x)`);
grdef(`DGA{^i} := g{^l ^j}*Sh{^i ,l ,j} - 2*A{^i ^j}*Lapse{,j} + 2*Lapse*(MCh{^i j k}*A{^k ^j} - 2/3*g{^i ^j}*Extc{,j} - g{^i ^j}*S{j} + 6*A{^i ^j}*DCF{j}) + Sh{^j}*GA{^i ,j} - GA{^j}*Sh{^i ,j} +vee/3*(2*GA{^i}*divbeta(t,x)+g{^i ^l}*DBS{,l})`);
grcalc(DGA(up));

gralter(DGA(up),simplify,expand,trigsin,factor,expand); 
grdisplay(DGA(up));
grdef(`met{a b} := exp(4*CF)*(kdelta{a $x}*kdelta{b $x}*A(t,x) + kdelta{a $theta}*kdelta{b $theta}*B(t,x)*x^2 + kdelta{a $phi}*kdelta{b $phi}*B(t,x)*x^2*sin(theta)^2)`);
grcalc(met(dn,dn));
grdisplay(met(dn,dn));

grdef(`met_s{a b} := kdelta{a $x}*kdelta{b $x}*metxx(t,x) + kdelta{a $theta}*kdelta{b $theta}*x^2*metthetatheta(t,x) + kdelta{a $phi}*kdelta{b $phi}*x^2*(sin(theta))^2*metthetatheta(t,x)`);
grcalc(met_s(dn,dn));grdisplay(met_s(dn,dn));
grdef(`invmet_s{^a ^b}:=kdelta{^a $x}*kdelta{^b $x}*1/metxx(t,x) + kdelta{^a $theta}*kdelta{^b $theta}*1/(x^2*metthetatheta(t,x)) + kdelta{^a $phi}*kdelta{^b $phi}*1/(x^2*(sin(theta))^2*metthetatheta(t,x))`);
grcalc(invmet_s(up,up));
grdisplay(invmet_s(up,up));
grdef(`rChL{d a b} := 1/2*( met_s{b d ,a} + met_s{a d,b} - met_s{a b,d})`);
grdef(`rMCh{^c a b} := invmet_s{^c ^d} * rChL{d a b}`);
grdef(`DDL{b c} :=   Lapse{,c ,b} - rMCh{^d b c}*Lapse{,d}`);
grcalc(DDL(dn,dn));
grdisplay(DDL(dn,dn));
grdef(`rR{u p} := rMCh{^v u p ,v} - rMCh{^v v p ,u} + rMCh{^a u p}*rMCh{^v a v} - rMCh{^a v p}*rMCh{^v a u}`);
grdef(`RHO{} := rho(t,x)`);
grcalc(RHO);
grdef(`TraceS{} := TS(t,x)`);
grcalc(TraceS);
grdef(`DDL_s{a b} := kdelta{a $x}*kdelta{b $x}*DDLxx(t,x) + kdelta{a $theta}*kdelta{b $theta}*x^2*DDLthetatheta(t,x) + kdelta{a $phi}*kdelta{b $phi}*x^2*DDLthetatheta(t,x)*(sin(theta))^2`);
grcalc(DDL_s(dn,dn)); grdisplay(DDL_s(dn,dn));
grdef(`TD{} := invmet_s{^i ^j} * DDL_s{i j}`);
grcalc(TD);
grdisplay(TD);
grdef(`DDLTF{i j} := DDL_s{i j} - 1/3*TD*met_s{i j}`);
grcalc(DDLTF(dn,dn));
grdisplay(DDLTF(dn,dn));

grdef(`DTK{} := -TD +Lapse*(A{i j}*A{^i ^j} +1/3*Extc^2) +1/2*Lapse*(RHO+TraceS) +Sh{^i}*Extc{,i}`);
grcalc(DTK);
gralter(DTK,simplify,expand);
grdisplay(DTK);
grdef(`JS{a b} := kdelta{a $x}*kdelta{b $x}*JSxx(t,x) + kdelta{a $theta}*kdelta{b $theta}*x^2*JSthth(t,x) + kdelta{a $phi}*kdelta{b $phi}*x^2*(sin(theta))^2*JSthth(t,x)`);
grcalc(JS(dn,dn));
grdisplay(JS(dn,dn));
grdef(`TJS := invmet_s{^i ^j}*JS{i j}`);
grcalc(TJS);
gralter(TJS,simplify,expand);
grdisplay(TJS);
grdef(`U{i j} := (Extc*A{i j} - 2*A{i l}*A{^l j})`);
grcalc(U(dn,dn));
gralter(U(dn,dn),simplify,expand);
grdisplay(U(dn,dn));
grdef(`dg{i j} := -2*Lapse*A{i j} + Q{i j}`);
grcalc(dg(dn,dn));
gralter(dg(dn,dn),simplify,factor);
grdisplay(dg(dn,dn));
grdef(`dgu{^i ^j} := -g{^i ^l}*dg{l m}*g{^m ^j}`); grcalc(dgu(up,up)); grdisplay(dgu(up,up));
grdef(`deltaDGA{^i} := -dgu{^j ^k}*MCh_f{^i j k}`); grcalc(deltaDGA(up)); grdisplay(deltaDGA(up));
grdef(`DLAM{^i} := DGA{^i} + deltaDGA{^i}`); grcalc(DLAM(up)); gralter(DLAM(up),simplify,expand); grdisplay(DLAM(up));
grdef(`DPH{} := -1/6*Lapse*Extc + Sh{^i}*CF{,i} + 1/6*divbeta(t,x)`);
grcalc(DPH); grdisplay(DPH);
grdef(`TRB{} := g{^i ^j}*RB{i j}`);
grcalc(TRB);
gralter(TRB,simplify);
grdisplay(TRB);
grdef(`PSI{} := psi(t,x)`);
grcalc(PSI);
grdef(`C1 {} := -1/8*TRB`); grcalc(C1); gralter(C1,simplify,expand,factor,expand); grdisplay(C1);
grdef(`C5{} := 1/8*A{i j}*A{^i ^j} -1/12*Extc^2 + 1/4*rho(t,x)`);
grcalc(C5); gralter(C5,simplify,expand); grdisplay(C5);
grdef(`HS3{} := g{^i ^j}*PSI{;j ;i} + C5psi(t,x) + C1s(t,x)*PSI`);
grdef(`HS{} := g{^i ^j}*PSI{;j ;i} + C5s(t,x)*PSI^5 + C1s(t,x)*PSI`); grdef(`HSR{} := g{^i ^j}*PSI{;i ;j} - 1/8*TRB*PSI + 1/8*A{i j}*A{^i ^j}*PSI^5 -1/12*PSI^5*Extc^2 + 1/4*rho(t,x)*PSI^5`); grcalc(HSR);gralter(HSR,simplify,expand,factor,expand);grdisplay(HSR);
grcalc(HS);
gralter(HS,simplify,expand,factor,expand);
grdisplay(HS);
grcalc(HS3); grdisplay(HS3);gralter(HS3,expand);
grdef(`RR{i j} := RB{i j} + Rphi{i j}`); 
grcalc(RR(dn,dn));
gralter(RR(dn,dn),simplify,expand,factor,expand);
grdisplay(RR(dn,dn));
grdef(`RR_s{a b} := kdelta{a $x}*kdelta{b $x}*RRxx(t,x) + kdelta{a $theta}*kdelta{b $theta}*x^2*RRthetatheta(t,x) + kdelta{a $phi}*kdelta{b $phi}*x^2*(sin(theta))^2*RRthetatheta(t,x)`);
grcalc(RR_s(dn,dn)); grdisplay(RR_s(dn,dn));
grdef(`TR{} := invmet_s{^i ^j}*RR_s{i j}`);
grcalc(TR);
grdisplay(TR);

grdef(`RRTF{i j} := RR_s{i j} - 1/3*met_s{i j}*TR`); 
grcalc(RRTF(dn,dn));
grdisplay(RRTF(dn,dn));
grdef(`P_s{a b} := kdelta{a $x}*kdelta{b $x}*Pxx(t,x) + kdelta{a $theta}*kdelta{b $theta}*x^2*Pthetatheta(t,x) + kdelta{a $phi}*kdelta{b $phi}*x^2*(sin(theta))^2*Pthetatheta(t,x)`);
grcalc(P_s(dn,dn)); grdisplay(P_s(dn,dn));
grdef(`U_s{a b} := kdelta{a $x}*kdelta{b $x}*Uxx(t,x) + kdelta{a $theta}*kdelta{b $theta}*x^2*Uthetatheta(t,x) + kdelta{a $phi}*kdelta{b $phi}*x^2*(sin(theta))^2*Uthetatheta(t,x)`);
grcalc(U_s(dn,dn));
grdisplay(U_s(dn,dn));
grdef(`DAA{i j} := em4phi(t,x)*(-DDLTF{i j} +  Lapse*(RRTF{i j} - ( JS{i j} - 1/3*met_s{i j}*TJS)    )) + Lapse*U_s{i j} + P_s{i j}`);
grcalc(DAA(dn,dn));
gralter(DAA(dn,dn),simplify,expand);
grdisplay(DAA(dn,dn));
grdef(`EMPH{} := exp(-4*CF)`);grcalc(EMPH);grdisplay(EMPH);
grdef(`MKTR1{^j ^i} := PSI^6*A{^j ^i}`); 
grdef(`MK{^i} := MKTR1{^j ^i ;j} - 2/3*PSI^6*g{^i ^j}*Extc{;j} - PSI^6*g{^i ^j}*S{j}`); grcalc(MK(up)); gralter(MK(up),simplify,expand); grdisplay(MK(up));
grdef(`DB:= Sh{^i ;i}`); grcalc(DB); grdisplay(DB);gralter(DB,simplify,expand); grdisplay(DB);


