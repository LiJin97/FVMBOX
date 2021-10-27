function [Cpk,Cpa,Cpb,Cpc,Cpd,temp,mK] = CalCoeffb(Pk,Pa,Pb,Pc,Pd)

%%
%This is a function to obtain the discretization of flux on a face of
%boundary cells.
%%
Po = (Pa+Pb+Pc+Pd)/4;

tOK = Pk-Po;

nK = cross(Pc-Pa,Pd-Pb);
Ak = norm(nK)/2;
nK = nK/norm(nK);

tAC = (Pc-Pa) - (nK)*dot(Pc-Pa,nK);
tBD = (Pd-Pb) - (nK)*dot(Pd-Pb,nK);

nAC = cross(tAC,nK);
nBD = cross(tBD,nK);

kappk = Genk(Pk);


a1 = dot(-kappk*nK,nK)/dot(tOK,nK);

a2 = dot(-kappk*nK-a1*tOK,nBD)/dot(tAC,nBD);

a3 = dot(-kappk*nK-a1*tOK,nAC)/dot(tBD,nAC);


Cpk = Ak*a1;
Cpa = -Ak*a2;
Cpb = -Ak*a3;
Cpc = -Cpa;
Cpd = -Cpb;

temp = -Cpk*GenReal(Po);

mK = abs(dot(cross(Pc-Pa,Pd-Pb),Pk-Po))/6;
