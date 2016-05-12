# real 3 metric for BSSN in spherical symmetry.
Ndim_ := 3:
x1_ := x:
x2_ := theta:
x3_ := phi:
complex_ := {}:
g11_ := A(x)*exp(4*phi(x)):
g12_ := 0:
g13_ := 0:
g22_ := x^2*B(x)*exp(4*phi(x)):
g23_ := 0:
g33_ := x^2*B(x)*(sin(theta))^2*exp(4*phi(x)):
