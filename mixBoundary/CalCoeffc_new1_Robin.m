function [Cpk,Cpa,Cpb,Cpc,Cpd,temp,mK] = CalCoeffc_new1_Robin(Pk,Pl,Pa,Pb,Pc,Pd)
temp = 0;
Po = (Pa+Pb+Pc+Pd)/4;
% Po=(Pk+Pl)/2;

tOK = Pk-Po;
tOL = Pl-Po;

nK = cross(Pc-Pa,Pd-Pb);
Ak = norm(nK)/2;
nK = nK/norm(nK);
nL = -nK;

tAC = (Pc-Pa) - (nK)*dot(Pc-Pa,nK);
tBD = (Pd-Pb) - (nK)*dot(Pd-Pb,nK);
% Ak = norm(cross(Pa-Pb,Pb-Pc))/2+norm(cross(Pa-Pd,Pd-Pc))/2;

nAC = cross(tAC,nK);
nBD = cross(tBD,nK);

kappk = Genk_Robin(Po);
kappl = Genk_Robin(Po);%£¿

% -kappk*nK = a1*tOK+a2*tAC+a3*tBD
% -kappl*nL = b1*tOL+b2*tAC+b3*tBD

a1 = dot(-kappk*nK,nK)/dot(tOK,nK);
b1 = dot(-kappl*nL,nL)/dot(tOL,nL);

a2 = dot(-kappk*nK-a1*tOK,nBD)/dot(tAC,nBD);
b2 = dot(-kappl*nL-b1*tOL,nBD)/dot(tAC,nBD);

a3 = dot(-kappk*nK-a1*tOK,nAC)/dot(tBD,nAC);
b3 = dot(-kappl*nL-b1*tOL,nAC)/dot(tBD,nAC);


%a1*(uk-uo)+a2*(uc-ua)+a3*(ud-ub) + b1*(ul-uo)+b2*(uc-ua)+b3*(ud-ub) = 0
%uo = (a1*uk+b1*ul)/(a1+b1) + (a2+b2)*(uc-ua)/(a1+b1) + (a3+b3)*(ud-ub)/(a1+b1) 
%Fk = a1*b1*(uk-ul)/(a1+b1) + (a2*b1-a1*b2)*(uc-ua)/(a1+b1) + (a3*b1-a1*b3)*(ud-ub)/(a1+b1)
Cpk = Ak*a1*b1/(a1+b1);
Cpa = -Ak*(a2*b1-a1*b2)/(a1+b1);
Cpb = -Ak*(a3*b1-a1*b3)/(a1+b1);
Cpc = -Cpa;
Cpd = -Cpb;

if Pa(1) <=1e-6
    temp = temp+Cpa*GenReal_Robin(Pa);
    Cpa = 0;
end

if Pb(1) <=1e-6
    temp = temp+Cpb*GenReal_Robin(Pb);
    Cpb = 0;
end

if Pc(1) <=1e-6
    temp = temp+Cpc*GenReal_Robin(Pc);
    Cpc = 0;
end

if Pd(1)  <=1e-6
    temp = temp+Cpd*GenReal_Robin(Pd);
    Cpd = 0;
end

mK = abs(dot(cross(Pc-Pa,Pd-Pb),Pk-Po))/6;


