############################################################
# Generic conformal 3 metric for BSSN in spherical symmetry
# See: BSSN_Spher.mw or BSSN_Spher.mpl
############################################################

Ndim_ := 3:
x1_ := x:
x2_ := theta:
x3_ := phi:
complex_ := {}:
g11_ := (A(t,x)):
g12_ := 0:
g13_ := 0:
g22_ := x^2*(B(t,x)):
g23_ := 0:
g33_ := x^2*(B(t,x))*(sin(theta))^2:
