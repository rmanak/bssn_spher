DPH_expr:=-1/6*alpha(t,x)*K(t,x)+beta(t,x)*diff(phi(t,x),x)+1/6*diff(beta(t,x),x);
dgxx_expr:=-2*alpha(t,x)*Axx(t,x)+beta(t,x)*diff(A(t,x),x)+4/3*A(t,x)*diff(beta(t,x),x);
dgthetatheta_expr:=-1/3*x*(6*alpha(t,x)*x*Athth(t,x)-6*beta(t,x)*B(t,x)-3*beta(t,x)*x*diff(B(t,x),x)+2*x*B(t,x)*diff(beta(t,x),x));
metxx_expr:=exp(phi(t,x))^4*A(t,x);
metthetatheta_expr:=exp(phi(t,x))^4*B(t,x);
invmetxx_expr:=1/exp(phi(t,x))^4/A(t,x);
invmetthetatheta_expr:=1/exp(phi(t,x))^4/B(t,x);
DDLxx_expr:=-1/2*(-2*diff(diff(alpha(t,x),x),x)*metxx(t,x)+diff(metxx(t,x),x)*diff(alpha(t,x),x))/metxx(t,x);
DDLthetatheta_expr:=1/2*x*(2*metthetatheta(t,x)+x*diff(metthetatheta(t,x),x))/metxx(t,x)*diff(alpha(t,x),x);
DTK_expr:=-1/metxx(t,x)*DDLxx(t,x)-2/metthetatheta(t,x)*DDLthetatheta(t,x)+1/3*alpha(t,x)*K(t,x)^2+1/2*alpha(t,x)*rho(t,x)+1/2*alpha(t,x)*TS(t,x)+beta(t,x)*diff(K(t,x),x)+1/A(t,x)^2*alpha(t,x)*Axx(t,x)^2+2/B(t,x)^2*alpha(t,x)*Athth(t,x)^2;
RRxx_expr:=3/4/A(t,x)^2*diff(A(t,x),x)^2-2/x^2-2/x/B(t,x)*diff(B(t,x),x)-1/2/B(t,x)^2*diff(B(t,x),x)^2+A(t,x)*diff(Lamx(t,x),x)+2/x^2*A(t,x)/B(t,x)+2/x*A(t,x)/B(t,x)^2*diff(B(t,x),x)+1/2*diff(A(t,x),x)*Lamx(t,x)-1/x/B(t,x)*diff(A(t,x),x)-1/2/A(t,x)*diff(diff(A(t,x),x),x)-4*diff(diff(phi(t,x),x),x)+2/A(t,x)*diff(phi(t,x),x)*diff(A(t,x),x)-4/x*diff(phi(t,x),x)-2/B(t,x)*diff(phi(t,x),x)*diff(B(t,x),x);
RRthetatheta_expr:=1/A(t,x)*B(t,x)+1/2/A(t,x)/B(t,x)*diff(B(t,x),x)^2*x^2+Lamx(t,x)*x*B(t,x)+1/2*x^2*Lamx(t,x)*diff(B(t,x),x)-1/B(t,x)*x*diff(B(t,x),x)-1/2/A(t,x)*x^2*diff(diff(B(t,x),x),x)-1-6/A(t,x)*B(t,x)*x*diff(phi(t,x),x)-3/A(t,x)*x^2*diff(phi(t,x),x)*diff(B(t,x),x)-2/A(t,x)*B(t,x)*x^2*diff(diff(phi(t,x),x),x)+1/A(t,x)^2*B(t,x)*x^2*diff(phi(t,x),x)*diff(A(t,x),x)-4/A(t,x)*B(t,x)*x^2*diff(phi(t,x),x)^2;
Uxx_expr:=alpha(t,x)*Axx(t,x)*K(t,x)-2*alpha(t,x)*Axx(t,x)^2/A(t,x);
Uthetatheta_expr:=alpha(t,x)*x^2*Athth(t,x)*K(t,x)-2*alpha(t,x)*x^2*Athth(t,x)^2/B(t,x);
Pxx_expr:=beta(t,x)*diff(Axx(t,x),x)+4/3*Axx(t,x)*diff(beta(t,x),x);
Pthetatheta_expr:=2*beta(t,x)*x*Athth(t,x)+beta(t,x)*x^2*diff(Athth(t,x),x)-2/3*x^2*Athth(t,x)*diff(beta(t,x),x);
DAAxx_expr:=-2/3*em4phi(t,x)*DDLxx(t,x)+2/3/metthetatheta(t,x)*em4phi(t,x)*DDLthetatheta(t,x)*metxx(t,x)+2/3*em4phi(t,x)*alpha(t,x)*RRxx(t,x)-2/3/metthetatheta(t,x)*em4phi(t,x)*alpha(t,x)*RRthetatheta(t,x)*metxx(t,x)-2/3*em4phi(t,x)*alpha(t,x)*JSxx(t,x)+2/3/metthetatheta(t,x)*em4phi(t,x)*alpha(t,x)*JSthth(t,x)*metxx(t,x)+alpha(t,x)*Uxx(t,x)+Pxx(t,x);
DAAthetatheta_expr:=1/3*x^2/metxx(t,x)*em4phi(t,x)*DDLxx(t,x)*metthetatheta(t,x)-1/3*x^2*em4phi(t,x)*DDLthetatheta(t,x)-1/3*x^2/metxx(t,x)*em4phi(t,x)*alpha(t,x)*RRxx(t,x)*metthetatheta(t,x)+1/3*x^2*em4phi(t,x)*alpha(t,x)*RRthetatheta(t,x)-1/3*x^2*em4phi(t,x)*alpha(t,x)*JSthth(t,x)+1/3*x^2/metxx(t,x)*em4phi(t,x)*alpha(t,x)*JSxx(t,x)*metthetatheta(t,x)+x^2*alpha(t,x)*Uthetatheta(t,x)+x^2*Pthetatheta(t,x);
DLAMx_expr:=-2/A(t,x)^2*Axx(t,x)*diff(alpha(t,x),x)-4/3/A(t,x)*alpha(t,x)*diff(K(t,x),x)-2/A(t,x)*alpha(t,x)*Sx(t,x)+12/A(t,x)^2*alpha(t,x)*Axx(t,x)*diff(phi(t,x),x)+beta(t,x)*diff(Lamx(t,x),x)-2/x^2/B(t,x)*beta(t,x)-1/3*diff(beta(t,x),x)*Lamx(t,x)+2/x/B(t,x)*diff(beta(t,x),x)+1/A(t,x)^3*alpha(t,x)*diff(A(t,x),x)*Axx(t,x)+4/3/A(t,x)*diff(diff(beta(t,x),x),x)-4/A(t,x)/x/B(t,x)*alpha(t,x)*Athth(t,x)-2/A(t,x)/B(t,x)^2*alpha(t,x)*Athth(t,x)*diff(B(t,x),x)+4/x/B(t,x)^2*alpha(t,x)*Athth(t,x);
DLAMtheta_expr:=0;
C1_expr:=-3/32/A(t,x)^3*diff(A(t,x),x)^2+1/4/x/A(t,x)/B(t,x)*diff(B(t,x),x)-1/16/A(t,x)/B(t,x)^2*diff(B(t,x),x)^2-1/8*diff(Lamx(t,x),x)-1/16/A(t,x)*diff(A(t,x),x)*Lamx(t,x)+1/8/x/A(t,x)/B(t,x)*diff(A(t,x),x)+1/16/A(t,x)^2*diff(diff(A(t,x),x),x)-1/4/x*Lamx(t,x)-1/8/B(t,x)*Lamx(t,x)*diff(B(t,x),x)+1/8/A(t,x)/B(t,x)*diff(diff(B(t,x),x),x);
C5_expr:=-1/12*K(t,x)^2+1/4*rho(t,x)+1/8/A(t,x)^2*Axx(t,x)^2+1/4/B(t,x)^2*Athth(t,x)^2;
HS_expr:=C5s(t,x)*psi(t,x)^5+C1s(t,x)*psi(t,x)+1/A(t,x)*diff(diff(psi(t,x),x),x)-1/2/A(t,x)^2*diff(psi(t,x),x)*diff(A(t,x),x)+2/A(t,x)/x*diff(psi(t,x),x)+1/A(t,x)/B(t,x)*diff(psi(t,x),x)*diff(B(t,x),x);

