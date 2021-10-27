function g = Gengrad(P)
x = P(1);
y = P(2);
z = P(3);
 g = [y*exp(x*y),x*exp(x*y),2*z];