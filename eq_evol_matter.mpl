# Generated using bssn_4d_spher_matter.mw
eq_evol_rPIb2:= 
diff(rPIb2(t,x) ,t) = 
expand( simplify (
(-mass*rpsi(t,x) + 2*(diff(rpsi(t, x), x))/(a1(t, x)^2*x)+(diff(diff(rpsi(t, x), x), x))/a1(t, x)^2-(diff(rpsi(t, x), x))*(diff(a1(t, x), x))/a1(t, x)^3+2*rPIb2(t, x)*beta(t, x)/(alpha(t, x)*b1(t, x)^2*a1(t, x)*x)+beta(t, x)*(diff(rPIb2(t, x), x))/(alpha(t, x)*b1(t, x)^2*a1(t, x))+rPIb2(t, x)*(diff(beta(t, x), x))/(alpha(t, x)*b1(t, x)^2*a1(t, x))+2*(diff(rpsi(t, x), x))*(diff(b1(t, x), x))/(b1(t, x)*a1(t, x)^2)+(diff(rpsi(t, x), x))*(diff(alpha(t, x), x))/(alpha(t, x)*a1(t, x)^2) ) * (alpha(t,x)*b1(t,x)^2*a1(t,x))
) );

eq_evol_iPIb2:=
(diff(iPIb2(t, x), t)) = expand( simplify (

( -mass*ipsi(t,x) + 2*(diff(ipsi(t, x), x))/(a1(t, x)^2*x)+(diff(diff(ipsi(t, x), x), x))/a1(t, x)^2-(diff(ipsi(t, x), x))*(diff(a1(t, x), x))/a1(t, x)^3+2*iPIb2(t, x)*beta(t, x)/(alpha(t, x)*b1(t, x)^2*a1(t, x)*x)+iPIb2(t, x)*(diff(beta(t, x), x))/(alpha(t, x)*b1(t, x)^2*a1(t, x))+2*(diff(ipsi(t, x), x))*(diff(b1(t, x), x))/(b1(t, x)*a1(t, x)^2)+(diff(ipsi(t, x), x))*(diff(alpha(t, x), x))/(alpha(t, x)*a1(t, x)^2)+beta(t, x)*(diff(iPIb2(t, x), x))/(alpha(t, x)*b1(t, x)^2*a1(t, x))) * (alpha(t,x)*b1(t,x)^2*a1(t,x))
) );


eq_evol_rpsi := diff(rpsi(t,x),t) = alpha(t,x)/a1(t,x)*rPI(t,x)+beta(t,x)*rPHI(t,x);

eq_evol_ipsi := diff(ipsi(t,x),t) = alpha(t,x)/a1(t,x)*iPI(t,x)+beta(t,x)*iPHI(t,x);
