function kappa = Genk(P)


%discontinusou tensor
kappa = (P(1)<=0.5)*eye(3)+(P(1)>0.5)*[2,1,0;1,2,0;0,0,1];
