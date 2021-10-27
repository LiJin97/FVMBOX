function [Cpk,Cpa,Cpb,Cpc,Cpd,temp,mK] = CalCoeffb(Pk,Pa,Pb,Pc,Pd)


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

kappk = Genk_Robin(Pk);

% -kappk*nK = a1*tOK+a2*tAC+a3*tBD

a1 = dot(-kappk*nK,nK)/dot(tOK,nK);

a2 = dot(-kappk*nK-a1*tOK,nBD)/dot(tAC,nBD);

a3 = dot(-kappk*nK-a1*tOK,nAC)/dot(tBD,nAC);


%a1*(uk-uo)+a2*(uc-ua)+a3*(ud-ub) + b1*(ul-uo)+b2*(uc-ua)+b3*(ud-ub) = 0
%uo = (a1*uk+b1*ul)/(a1+b1) + (a2+b2)*(uc-ua)/(a1+b1) + (a3+b3)*(ud-ub)/(a1+b1) 
%Fk = a1*b1*(uk-ul)/(a1+b1) + (a2*b1-a1*b2)*(uc-ua)/(a1+b1) + (a3*b1-a1*b3)*(ud-ub)/(a1+b1)
Cpk = Ak*a1;
Cpa = -Ak*a2;
Cpb = -Ak*a3;
Cpc = -Cpa;
Cpd = -Cpb;

temp = -Cpk*GenReal_Robin(Po) + Cpa*GenReal_Robin(Pa) + Cpb*GenReal_Robin(Pb)+...
    Cpc*GenReal_Robin(Pc)+Cpd*GenReal_Robin(Pd);
Cpa = 0;
Cpb = 0;
Cpc = 0;
Cpd = 0;

mK = abs(dot(cross(Pc-Pa,Pd-Pb),Pk-Po))/6;
