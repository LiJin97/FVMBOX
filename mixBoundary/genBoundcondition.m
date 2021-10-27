function g = genBoundcondition(P,id)
alpha = 1;
beta = 1;
kappa = Genk_Robin(P);
a = kappa(1,1);
b = kappa(2,2);
c = kappa(1,2);
d = kappa(3,3);
x = P(1);
y = P(2);
z = P(3);
if id == 2
    g = (c*x+a*y+1)*exp(x*y)+z^2;
%    outv  = - cross([0,0,1],[0,1,0])/norm(cross([0,0,1],[0,1,0]));
elseif id == 1
    g = GenReal_Robin([x,y,z]);
   % outv = cross([0,0,1],[0,1,0])/norm(cross([0,0,1],[0,1,0]));
elseif id == 4
    g = (b*x+c*y+1)*exp(x*y)+z^2;
   % outv = cross([-1,0,0],[0,0,1])/norm(cross([-1,0,0],[0,0,1]));
elseif id == 3
    g = (-b*x-c*y+1)*exp(x*y)+z^2;
   % outv = -cross([-1,0,0],[0,0,1])/norm(cross([-1,0,0],[0,0,1]));
elseif id == 6
    g = exp(x*y)+z^2+2*d*z;
   % outv = cross([-1,0,0],[0,-1,0])/norm(cross([-1,0,0],[0,-1,0]));
elseif id == 5
    g = exp(x*y)+z^2-2*d*z;
   % outv = -cross([-1,0,0],[0,-1,0])/norm(cross([-1,0,0],[0,-1,0]));
end
%outv = outv';