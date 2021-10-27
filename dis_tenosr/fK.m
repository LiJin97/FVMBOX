function f = fK(P)
x = P(1);
y = P(2);
z = P(3);


f = ((x-0.5)*(2*sin(y)+cos(y))+sin(y))*(x<=0.5)+(-2*(exp(x-0.5))*cos(y))*(x>0.5);
 