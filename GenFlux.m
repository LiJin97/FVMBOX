function  flux = GenFlux(Pk,Pa,Pb,Pc,Pd)


Po = (Pa+Pb+Pc+Pd)/4;
tAC = Pc-Pa;
tBD = Pd-Pb;
n = cross(tAC,tBD);

if dot(Pk-Po,n)>0
    dd = 1;
else
    dd = -1;
end
t =1;
kappa = Genk(Po*(1-t)+Pk*t);
g = Gengrad(Po*(1-t)+Pk*t);
flux = dd*g*kappa*n/2;
