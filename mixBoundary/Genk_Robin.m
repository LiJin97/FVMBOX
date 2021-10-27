function kappa = Genk_Robin(P)
theta = pi/6;
kappa = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1]*...
    [1,0,0;0,0.1,0;0,0,1]*[cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1];