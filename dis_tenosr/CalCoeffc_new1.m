function [Cpk,Cpa,Cpb,Cpc,Cpd,mK] = CalCoeffc_new1(Pk,Pl,Pa,Pb,Pc,Pd)
%%
%This is a function to obtain the dicretization of flux on an interior
%face.
%%
Po = (Pa+Pb+Pc+Pd)/4;

tOK = Pk-Po;
tOL = Pl-Po;

nK = cross(Pc-Pa,Pd-Pb);
Ak = norm(nK)/2;
nK = nK/norm(nK);
nL = -nK;

tAC = (Pc-Pa) - (nK)*dot(Pc-Pa,nK);
tBD = (Pd-Pb) - (nK)*dot(Pd-Pb,nK);

nAC = cross(tAC,nK);
nBD = cross(tBD,nK);

kappk = Genk(Pk);
kappl = Genk(Pl);



a1 = dot(-kappk*nK,nK)/dot(tOK,nK);
b1 = dot(-kappl*nL,nL)/dot(tOL,nL);

a2 = dot(-kappk*nK-a1*tOK,nBD)/dot(tAC,nBD);
b2 = dot(-kappl*nL-b1*tOL,nBD)/dot(tAC,nBD);

a3 = dot(-kappk*nK-a1*tOK,nAC)/dot(tBD,nAC);
b3 = dot(-kappl*nL-b1*tOL,nAC)/dot(tBD,nAC);



Cpk = Ak*a1*b1/(a1+b1);
Cpa = -Ak*(a2*b1-a1*b2)/(a1+b1);
Cpb = -Ak*(a3*b1-a1*b3)/(a1+b1);
Cpc = -Cpa;
Cpd = -Cpb;

mK = abs(dot(cross(Pc-Pa,Pd-Pb),Pk-Po))/6;


