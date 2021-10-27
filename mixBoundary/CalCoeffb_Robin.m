function [Cpk,Cpa,Cpb,Cpc,Cpd,temp,mK] = CalCoeffb_Robin(Pk,Pa,Pb,Pc,Pd)
Po = (Pa+Pb+Pc+Pd)/4;
eps = 1e-10;
if abs(Po(1))<eps
    id = 1;
elseif abs(Po(1)-1)<eps
    id = 2;
elseif abs(Po(2))<eps
    id = 3;
elseif abs(Po(2)-1)<eps
    id = 4;
elseif abs(Po(3))<eps
    id = 5;
elseif abs(Po(3)-1)<eps
    id = 6;
end
gPo = genBoundcondition(Po,id);
s = norm(cross(Pc-Pa,Pd-Pb))/2;
Cpa = 1/4*s;
Cpb = 1/4*s;
Cpc = 1/4*s;
Cpd = 1/4*s;
Cpk = 0;
temp = -gPo*s;

if Pa(1) == 0
    temp = temp+Cpa*GenReal_Robin(Pa);
    Cpa = 0;
end

if Pb(1) == 0
    temp = temp+Cpb*GenReal_Robin(Pb);
    Cpb = 0;
end

if Pc(1) == 0
    temp = temp+Cpc*GenReal_Robin(Pc);
    Cpc = 0;
end

if Pd(1) == 0
    temp = temp+Cpd*GenReal_Robin(Pd);
    Cpd = 0;
end

mK = abs(dot(cross(Pc-Pa,Pd-Pb),Pk-Po))/6;



% n = cross(Pc-Pa,Pd-Pb);
% nK = cross(Pc-Pa,Pd-Pb);
% Ak = norm(nK)/2;
% nK = nK/norm(nK);
% x = dot(Po-Pk,nK)/dot(k_K'*nK,nK);
% Pk1 = Pk+x*k_K'*nK;
% %Pk1 = findcrossP(Pk,k_K*n);
% eps = 0.000000001;
% if abs(Pk1(1))<eps
%     id = 1;
% elseif abs(Pk1(1)-1)<eps
%     id = 2;
% elseif abs(Pk1(2))<eps
%     id = 3;
% elseif abs(Pk1(2)-1)<eps
%     id = 4;
% elseif abs(Pk1(3))<eps
%     id = 5;
% elseif abs(Pk1(3)-1)<eps
%     id = 6;
% end
% gPk1 = genBoundcondition(Pk1,id);
% 
% alphaPk1 = 1;
% betaPk1 = 1;
% c1 = norm(k_K*n)*alphaPk1/(norm(k_K*n)*alphaPk1+norm(Pk1-Pk)*betaPk1);
% c2 = norm(Pk1-Pk)*gPk1/(norm(k_K*n)*alphaPk1+norm(Pk1-Pk)*betaPk1);
% 
% Cpk = Ak*norm(k_K*n)*(c1-1)/norm(Pk1-Pk);
% temp = Ak*norm(k_K*n)*c2/norm(Pk1-Pk);
% 
% mK = abs(dot(cross(Pc-Pa,Pd-Pb),Pk-Po))/6;
% Cpa = 0;
% Cpb = 0;
% Cpc = 0;
% Cpd = 0;


