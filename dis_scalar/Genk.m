function kappa = Genk(P)


if P(1)>= 0.5
    kappa = 1*eye(3);
else
    kappa = 5*eye(3);
end
