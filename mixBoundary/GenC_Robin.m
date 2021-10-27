function Ch = GenC_Robin(Pxh,Pyh,Pzh,xxh,yyh,zzh)
Nx = size(Pxh,1);
Ny = Nx;
Nz = Nx;
%%
xid = 1;
yid = 1;
zid = 1;
Pa = [xxh(xid,yid,zid),yyh(xid,yid,zid),zzh(xid,yid,zid)]';
P8c = [Pxh(xid,yid,zid),Pyh(xid,yid,zid),Pzh(xid,yid,zid)]';
nK1 = [-1;0;0];
nK2 = [0;-1;0];
nK3 = [0;0;-1];
Kapa = Genk_Robin(Pa);
vK1 = Kapa'*nK1;
vK2 = Kapa'*nK2;
vK3 = Kapa'*nK3;

% % P8c-Pa = omega1*Kapa'*nK1+omega2*Kapa'*nK2+omega3*Kapa'*nK3
% % u(P8c)-u(Pa) = dua (P8c-Pa) = dua (omega1*Kapa'*nK1+omega2*Kapa'*nK2+omega3*Kapa'*nK3)
% % u(P8c)-u(Pa) = omega1*Kapa*dua*nK1+omega2*Kapa*dua*nK2+omega3*Kapa*dua*nK3
% % = omega1*(g1-ua)+omega2*(g2-ua)+omega3*(g3-ua) +O(h^2)
% 
omega1 = dot(P8c-Pa,cross(vK2,vK3))/dot(vK1,cross(vK2,vK3));
omega2 = dot(P8c-Pa,cross(vK3,vK1))/dot(vK2,cross(vK3,vK1));
omega3 = dot(P8c-Pa,cross(vK1,vK2))/dot(vK3,cross(vK1,vK2));

id = 1;
g1 = genBoundcondition(Pa,id);
id = 3;
g2 = genBoundcondition(Pa,id);
id = 5;
g3 = genBoundcondition(Pa,id);
Ch(9,xid,yid,zid) = 0;
% Ch(8,xid,yid,zid) = -1/(omega1+omega2+omega3-1);
% Ch(9,xid,yid,zid) = (omega1*g1+omega2*g2+omega3*g3)/(omega1+omega2+omega3-1);
%  %err = abs(vh(xid,yid,zid) - (Ch(8,xid,yid,zid)*uh(xid,yid,zid)+Ch(9,xid,yid,zid)));
% % err = abs(( uh(xid,yid,zid) - vh(xid,yid,zid) ) - ( omega1*(g1-vh(xid,yid,zid))+omega2*(g2-vh(xid,yid,zid))+omega3*(g3-vh(xid,yid,zid)) ));
% %disp(err)
% 
%%
xid = 1;
yid = Nx+1;
zid = 1;
Pa = [xxh(xid,yid,zid),yyh(xid,yid,zid),zzh(xid,yid,zid)]';
P5c = [Pxh(xid,yid-1,zid),Pyh(xid,yid-1,zid),Pzh(xid,yid-1,zid)]';
nK1 = [-1;0;0];
nK2 = [0;1;0];
nK3 = [0;0;-1];
Kapa = Genk_Robin(Pa);
vK1 = Kapa'*nK1;
vK2 = Kapa'*nK2;
vK3 = Kapa'*nK3;
omega1 = dot(P5c-Pa,cross(vK2,vK3))/dot(vK1,cross(vK2,vK3));
omega2 = dot(P5c-Pa,cross(vK3,vK1))/dot(vK2,cross(vK3,vK1));
omega3 = dot(P5c-Pa,cross(vK1,vK2))/dot(vK3,cross(vK1,vK2));
id = 1;
g1 = genBoundcondition(Pa,id);
id = 4;
g2 = genBoundcondition(Pa,id);
id = 5;
g3 = genBoundcondition(Pa,id);
% Ch(5,xid,yid,zid) = -1/(omega1+omega2+omega3-1);
% Ch(9,xid,yid,zid) = (omega1*g1+omega2*g2+omega3*g3)/(omega1+omega2+omega3-1);
Ch(9,xid,yid,zid) = 0;
% % err = abs(vh(xid,yid,zid) - (Ch(5,xid,yid,zid)*uh(xid,yid-1,zid)+Ch(9,xid,yid,zid)));
% % err = abs(( uh(xid,yid-1,zid) - vh(xid,yid,zid) ) - ( omega1*(g1-vh(xid,yid,zid))+omega2*(g2-vh(xid,yid,zid))+omega3*(g3-vh(xid,yid,zid)) ));
% % disp(err)
%%
xid = 1;
yid = 1;
zid = Nz+1;
Pa = [xxh(xid,yid,zid),yyh(xid,yid,zid),zzh(xid,yid,zid)]';
P4c = [Pxh(xid,yid,zid-1),Pyh(xid,yid,zid-1),Pzh(xid,yid,zid-1)]';
nK1 = [-1;0;0];
nK2 = [0;-1;0];
nK3 = [0;0;1];
Kapa = Genk_Robin(Pa);
vK1 = Kapa'*nK1;
vK2 = Kapa'*nK2;
vK3 = Kapa'*nK3;
omega1 = dot(P4c-Pa,cross(vK2,vK3))/dot(vK1,cross(vK2,vK3));
omega2 = dot(P4c-Pa,cross(vK3,vK1))/dot(vK2,cross(vK3,vK1));
omega3 = dot(P4c-Pa,cross(vK1,vK2))/dot(vK3,cross(vK1,vK2));
id = 1;
g1 = genBoundcondition(Pa,id);
id = 3;
g2 = genBoundcondition(Pa,id);
id = 6;
g3 = genBoundcondition(Pa,id);
% Ch(4,xid,yid,zid) = -1/(omega1+omega2+omega3-1);
% Ch(9,xid,yid,zid) = (omega1*g1+omega2*g2+omega3*g3)/(omega1+omega2+omega3-1);
Ch(9,xid,yid,zid) = 0;
% % err = abs(vh(xid,yid,zid) - (Ch(4,xid,yid,zid)*uh(xid,yid,zid-1)+Ch(9,xid,yid,zid)));
% % err = abs(( uh(xid,yid,zid-1) - vh(xid,yid,zid) ) - ( omega1*(g1-vh(xid,yid,zid))+omega2*(g2-vh(xid,yid,zid))+omega3*(g3-vh(xid,yid,zid)) ));
% % disp(err)
%%
xid = 1;
yid = Ny+1;
zid = Nz+1;
Pa = [xxh(xid,yid,zid),yyh(xid,yid,zid),zzh(xid,yid,zid)]';
P1c = [Pxh(xid,yid-1,zid-1),Pyh(xid,yid-1,zid-1),Pzh(xid,yid-1,zid-1)]';
nK1 = [-1;0;0];
nK2 = [0;1;0];
nK3 = [0;0;1];
Kapa = Genk_Robin(Pa);
vK1 = Kapa'*nK1;
vK2 = Kapa'*nK2;
vK3 = Kapa'*nK3;
omega1 = dot(P1c-Pa,cross(vK2,vK3))/dot(vK1,cross(vK2,vK3));
omega2 = dot(P1c-Pa,cross(vK3,vK1))/dot(vK2,cross(vK3,vK1));
omega3 = dot(P1c-Pa,cross(vK1,vK2))/dot(vK3,cross(vK1,vK2));
id = 1;
g1 = genBoundcondition(Pa,id);
id = 4;
g2 = genBoundcondition(Pa,id);
id = 6;
g3 = genBoundcondition(Pa,id);
% Ch(1,xid,yid,zid) = -1/(omega1+omega2+omega3-1);
% Ch(9,xid,yid,zid) = (omega1*g1+omega2*g2+omega3*g3)/(omega1+omega2+omega3-1);
Ch(9,xid,yid,zid) = 0;
% % err = abs(vh(xid,yid,zid) - (Ch(1,xid,yid,zid)*uh(xid,yid-1,zid-1)+Ch(9,xid,yid,zid)));
% % err = abs(( uh(xid,yid,zid-1) - vh(xid,yid,zid) ) - ( omega1*(g1-vh(xid,yid,zid))+omega2*(g2-vh(xid,yid,zid))+omega3*(g3-vh(xid,yid,zid)) ));
% % disp(err)
%%
xid = Nx+1;
yid = 1;
zid = 1;
Pa = [xxh(xid,yid,zid),yyh(xid,yid,zid),zzh(xid,yid,zid)]';
P7c = [Pxh(xid-1,yid,zid),Pyh(xid-1,yid,zid),Pzh(xid-1,yid,zid)]';
nK1 = [1;0;0];
nK2 = [0;-1;0];
nK3 = [0;0;-1];
Kapa = Genk_Robin(Pa);
vK1 = Kapa'*nK1;
vK2 = Kapa'*nK2;
vK3 = Kapa'*nK3;
omega1 = dot(P7c-Pa,cross(vK2,vK3))/dot(vK1,cross(vK2,vK3));
omega2 = dot(P7c-Pa,cross(vK3,vK1))/dot(vK2,cross(vK3,vK1));
omega3 = dot(P7c-Pa,cross(vK1,vK2))/dot(vK3,cross(vK1,vK2));
id = 2;
g1 = genBoundcondition(Pa,id);
id = 3;
g2 = genBoundcondition(Pa,id);
id = 5;
g3 = genBoundcondition(Pa,id);
Ch(7,xid,yid,zid) = -1/(omega1+omega2+omega3-1);
Ch(9,xid,yid,zid) = (omega1*g1+omega2*g2+omega3*g3)/(omega1+omega2+omega3-1);

% %err = abs(vh(xid,yid,zid) - (Ch(7,xid,yid,zid)*uh(xid-1,yid,zid)+Ch(9,xid,yid,zid)));
% %err = abs(( uh(xid-1,yid,zid) - vh(xid,yid,zid) ) - ( omega1*(g1-vh(xid,yid,zid))+omega2*(g2-vh(xid,yid,zid))+omega3*(g3-vh(xid,yid,zid)) ));
% %disp(err)
% 
%%
xid = Nx+1;
yid = Ny+1;
zid = 1;
Pa = [xxh(xid,yid,zid),yyh(xid,yid,zid),zzh(xid,yid,zid)]';
P6c = [Pxh(xid-1,yid-1,zid),Pyh(xid-1,yid-1,zid),Pzh(xid-1,yid-1,zid)]';
nK1 = [1;0;0];
nK2 = [0;1;0];
nK3 = [0;0;-1];
Kapa = Genk_Robin(Pa);
vK1 = Kapa'*nK1;
vK2 = Kapa'*nK2;
vK3 = Kapa'*nK3;
omega1 = dot(P6c-Pa,cross(vK2,vK3))/dot(vK1,cross(vK2,vK3));
omega2 = dot(P6c-Pa,cross(vK3,vK1))/dot(vK2,cross(vK3,vK1));
omega3 = dot(P6c-Pa,cross(vK1,vK2))/dot(vK3,cross(vK1,vK2));
id = 2;
g1 = genBoundcondition(Pa,id);
id = 4;
g2 = genBoundcondition(Pa,id);
id = 5;
g3 = genBoundcondition(Pa,id);
Ch(6,xid,yid,zid) = -1/(omega1+omega2+omega3-1);
Ch(9,xid,yid,zid) = (omega1*g1+omega2*g2+omega3*g3)/(omega1+omega2+omega3-1);

% % err = abs(vh(xid,yid,zid) - (Ch(6,xid,yid,zid)*uh(xid-1,yid-1,zid)+Ch(9,xid,yid,zid)));
% % err = abs(( uh(xid-1,yid-1,zid) - vh(xid,yid,zid) ) - ( omega1*(g1-vh(xid,yid,zid))+omega2*(g2-vh(xid,yid,zid))+omega3*(g3-vh(xid,yid,zid)) ));
% % disp(err)
% 
%%
xid = Nx+1;
yid = 1;
zid = Nz+1;
Pa = [xxh(xid,yid,zid),yyh(xid,yid,zid),zzh(xid,yid,zid)]';
P3c = [Pxh(xid-1,yid,zid-1),Pyh(xid-1,yid,zid-1),Pzh(xid-1,yid,zid-1)]';
nK1 = [1;0;0];
nK2 = [0;-1;0];
nK3 = [0;0;1];
Kapa = Genk_Robin(Pa);
vK1 = Kapa'*nK1;
vK2 = Kapa'*nK2;
vK3 = Kapa'*nK3;
omega1 = dot(P3c-Pa,cross(vK2,vK3))/dot(vK1,cross(vK2,vK3));
omega2 = dot(P3c-Pa,cross(vK3,vK1))/dot(vK2,cross(vK3,vK1));
omega3 = dot(P3c-Pa,cross(vK1,vK2))/dot(vK3,cross(vK1,vK2));
id = 2;
g1 = genBoundcondition(Pa,id);
id = 3;
g2 = genBoundcondition(Pa,id);
id = 6;
g3 = genBoundcondition(Pa,id);
Ch(3,xid,yid,zid) = -1/(omega1+omega2+omega3-1);
Ch(9,xid,yid,zid) = (omega1*g1+omega2*g2+omega3*g3)/(omega1+omega2+omega3-1);

% % err = abs(vh(xid,yid,zid) - (Ch(3,xid,yid,zid)*uh(xid-1,yid,zid-1)+Ch(9,xid,yid,zid)));
% % err = abs(( uh(xid-1,yid,zid-1) - vh(xid,yid,zid) ) - ( omega1*(g1-vh(xid,yid,zid))+omega2*(g2-vh(xid,yid,zid))+omega3*(g3-vh(xid,yid,zid)) ));
% % disp(err)
%%
xid = Nx+1;
yid = Ny+1;
zid = Nz+1;
Pa = [xxh(xid,yid,zid),yyh(xid,yid,zid),zzh(xid,yid,zid)]';
P2c = [Pxh(xid-1,yid-1,zid-1),Pyh(xid-1,yid-1,zid-1),Pzh(xid-1,yid-1,zid-1)]';
nK1 = [1;0;0];
nK2 = [0;1;0];
nK3 = [0;0;1];
Kapa = Genk_Robin(Pa);
vK1 = Kapa'*nK1;
vK2 = Kapa'*nK2;
vK3 = Kapa'*nK3;
omega1 = dot(P2c-Pa,cross(vK2,vK3))/dot(vK1,cross(vK2,vK3));
omega2 = dot(P2c-Pa,cross(vK3,vK1))/dot(vK2,cross(vK3,vK1));
omega3 = dot(P2c-Pa,cross(vK1,vK2))/dot(vK3,cross(vK1,vK2));
id = 2;
g1 = genBoundcondition(Pa,id);
id = 4;
g2 = genBoundcondition(Pa,id);
id = 6;
g3 = genBoundcondition(Pa,id);
Ch(2,xid,yid,zid) = -1/(omega1+omega2+omega3-1);
Ch(9,xid,yid,zid) = (omega1*g1+omega2*g2+omega3*g3)/(omega1+omega2+omega3-1);

% % err = abs(vh(xid,yid,zid) - (Ch(2,xid,yid,zid)*uh(xid-1,yid-1,zid-1)+Ch(9,xid,yid,zid)));
% % % err = abs(( uh(xid-1,yid-1,zid-1) - vh(xid,yid,zid) ) - ( omega1*(g1-vh(xid,yid,zid))+omega2*(g2-vh(xid,yid,zid))+omega3*(g3-vh(xid,yid,zid)) ));
% % disp(err)
% 
%%
% 棱
%%
%epsilon = 0.0001;
%(2:Nx,2,2)
%c = 1;
for i = 2:Nx%:Nx
    for j = 1
        for k = 1
            W = zeros(4,6);
            Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
            Qu = [xxh(i,j,k+1),yyh(i,j,k+1),zzh(i,j,k+1)];
            Qr = [xxh(i,j+1,k),yyh(i,j+1,k),zzh(i,j+1,k)];
%             Q1 = Q + [0,epsilon,0];
%             Q2 = Q + [0,0,epsilon];
            tu = Qu - Q;
            tr = Qr - Q;
            K1 = [Pxh(i-1,j,k),Pyh(i-1,j,k),Pzh(i-1,j,k)];
            K2 = [Pxh(i,j,k),Pyh(i,j,k),Pzh(i,j,k)];
            k_K1 = Genk_Robin(K1);
            k_K2 = Genk_Robin(K2);
            T11 = [1,0,0;0,1,0;0,0,1];
            T12 = [k_K1*cross(tu,tr)',tu',tr']/[k_K2*cross(tu,tr)',tu',tr'];
            alphaQ = 1;
            betaQ = 1;
            v1 = [0;-1;0];
            v2 = [0;0;-1];
            id = 3;
            g1 = genBoundcondition(Q,id);
            id = 5;
            g2 = genBoundcondition(Q,id);
            W(1,1:2) = 1;
            W(1,3:6) = betaQ;
            W(2:4,1) = T11*(K1-Q)';
            W(2:4,2) = T12*(K2-Q)';
            W(2:4,3) = alphaQ*k_K1'*v1;
            W(2:4,4) = alphaQ*k_K1'*v2;
            W(2:4,5) = alphaQ*T12*k_K2'*v1;
            W(2:4,6) = alphaQ*T12*k_K2'*v2;
            ww = W'*((W*W')\[1;0;0;0]);
            Ch(7:8,i,j,k) = ww(1:2);
            Ch(9,i,j,k) = ww(3)*g1+ww(4)*g2+ww(5)*g1+ww(6)*g2;
%             err(c) = abs(vh(i,j,k) - (Ch(7,i,j,k)*uh(i-1,j,k)+...
%                 Ch(8,i,j,k)*uh(i,j,k)+Ch(9,i,j,k)));
%             c =  c+1;
        end
    end
end
% disp(max(err(:)))
% %%
% % (2:Nx,Ny,2)
% c = 1;
for i =2:Nx
    for j = Ny+1
        for k = 1
            W = zeros(4,6);
            Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
            Qu = [xxh(i,j,k+1),yyh(i,j,k+1),zzh(i,j,k+1)];
            Ql = [xxh(i,j-1,k),yyh(i,j-1,k),zzh(i,j-1,k)];
            tu = Qu - Q;
            tl = Ql - Q;
            K1 = [Pxh(i-1,j-1,k),Pyh(i-1,j-1,k),Pzh(i-1,j-1,k)];
            K2 = [Pxh(i,j-1,k),Pyh(i,j-1,k),Pzh(i,j-1,k)];
            k_K1 = Genk_Robin(K1);
            k_K2 = Genk_Robin(K2);
            T11 = [1,0,0;0,1,0;0,0,1];
            T12 = [k_K1*cross(tl,tu)',tl',tu']/[k_K2*cross(tl,tu)',tl',tu'];
            alphaQ = 1;
            betaQ = 1;
            v1 = [0;1;0];
            v2 = [0;0;-1];
            id = 4;
            g1 = genBoundcondition(Q,id);
            id = 5;
            g2 = genBoundcondition(Q,id);
            W(1,1:2) = 1;
            W(1,3:6) = betaQ;
            W(2:4,1) = T11*(K1-Q)';
            W(2:4,2) = T12*(K2-Q)';
            W(2:4,3) = alphaQ*k_K1'*v1;
            W(2:4,4) = alphaQ*k_K1'*v2;
            W(2:4,5) = alphaQ*T12*k_K2'*v1;
            W(2:4,6) = alphaQ*T12*k_K2'*v2;
            ww = W'*((W*W')\[1;0;0;0]);
            Ch(5:6,i,j,k) = ww(2:-1:1);
            Ch(9,i,j,k) = ww(3)*g1+ww(4)*g2+ww(5)*g1+ww(6)*g2;
%             err(c) = abs(vh(i,j,k) - (Ch(5,i,j,k)*uh(i,j-1,k)+...
%                 Ch(6,i,j,k)*uh(i-1,j-1,k)+Ch(9,i,j,k)));
%             c =  c+1;
        end
    end
end
% disp(max(err(:)))
% %%
% % (2:Nx,Ny,Nz)
% c = 1;
for i =2:Nx
    for j = 1
        for k = Nz+1
            W = zeros(4,6);
            Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
            Qd = [xxh(i,j,k-1),yyh(i,j,k-1),zzh(i,j,k-1)];
            Qr = [xxh(i,j+1,k),yyh(i,j+1,k),zzh(i,j+1,k)];
%             Q1 = Q + [0,epsilon,0];
%             Q2 = Q + [0,0,-epsilon];
            td = Qd - Q;
            tr = Qr - Q;
            K1 = [Pxh(i-1,j,k-1),Pyh(i-1,j,k-1),Pzh(i-1,j,k-1)];
            K2 = [Pxh(i,j,k-1),Pyh(i,j,k-1),Pzh(i,j,k-1)];
            k_K1 = Genk_Robin(K1);
            k_K2 = Genk_Robin(K2);
            T11 = [1,0,0;0,1,0;0,0,1];
            T12 = [k_K1*cross(td,tr)',td',tr']/[k_K2*cross(td,tr)',td',tr'];
            alphaQ = 1;
            betaQ = 1;
            id = 3;
            g1 = genBoundcondition(Q,id);
            id = 6;
            g2 = genBoundcondition(Q,id);
            v1 = [0;-1;0];
            v2 = [0;0;1];
            W(1,1:2) = 1;
            W(1,3:6) = betaQ;
            W(2:4,1) = T11*(K1-Q)';
            W(2:4,2) = T12*(K2-Q)';
            W(2:4,3) = alphaQ*k_K1'*v1;
            W(2:4,4) = alphaQ*k_K1'*v2;
            W(2:4,5) = alphaQ*T12*k_K2'*v1;
            W(2:4,6) = alphaQ*T12*k_K2'*v2;
            ww = W'*((W*W')\[1;0;0;0]);
            Ch(3:4,i,j,k) = ww(1:2);
            Ch(9,i,j,k) = ww(3)*g1+ww(4)*g2+ww(5)*g1+ww(6)*g2;
%             err(c) = abs(vh(i,j,k) - (Ch(4,i,j,k)*uh(i,j,k-1)+...
%                 Ch(3,i,j,k)*uh(i-1,j,k-1)+Ch(9,i,j,k)));
%             c =  c+1;
        end
    end
end
% disp(max(err(:)))
% % 
% %%
% % (2:Nx,Ny,Nz)
% c = 1;
for i =2:Nx
    for j = Ny+1
        for k = Nz+1
            W = zeros(4,6);
            Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
            Qd = [xxh(i,j,k-1),yyh(i,j,k-1),zzh(i,j,k-1)];
            Ql = [xxh(i,j-1,k),yyh(i,j-1,k),zzh(i,j-1,k)];
%             Q1 = Q + [0,-epsilon,0];
%             Q2 = Q + [0,0,-epsilon];
            td = Qd - Q;
            tl = Ql - Q;
            K1 = [Pxh(i-1,j-1,k-1),Pyh(i-1,j-1,k-1),Pzh(i-1,j-1,k-1)];
            K2 = [Pxh(i,j-1,k-1),Pyh(i,j-1,k-1),Pzh(i,j-1,k-1)];
            k_K1 = Genk_Robin(K1);
            k_K2 = Genk_Robin(K2);
            T11 = [1,0,0;0,1,0;0,0,1];
            T12 = [k_K1*cross(td,tl)',td',tl']/[k_K2*cross(td,tl)',td',tl'];
            alphaQ = 1;
            betaQ = 1;
            v1 = [0;1;0]; 
            v2 = [0;0;1];
            id = 4;
            g1 = genBoundcondition(Q,id);
            id = 6;
            g2 = genBoundcondition(Q,id);
            W(1,1:2) = 1;
            W(1,3:6) = betaQ;
            W(2:4,1) = T11*(K1-Q)';
            W(2:4,2) = T12*(K2-Q)';
            W(2:4,3) = alphaQ*k_K1'*v1;
            W(2:4,4) = alphaQ*k_K1'*v2;
            W(2:4,5) = alphaQ*T12*k_K2'*v1;
            W(2:4,6) = alphaQ*T12*k_K2'*v2;
            ww = W'*((W*W')\[1;0;0;0]);
            Ch(1:2,i,j,k) = ww(2:-1:1);
            Ch(9,i,j,k) = ww(3)*g1+ww(4)*g2+ww(5)*g1+ww(6)*g2;
%             err(c) = abs(vh(i,j,k) - (Ch(1,i,j,k)*uh(i,j-1,k-1)+...
%                 Ch(2,i,j,k)*uh(i-1,j-1,k-1)+Ch(9,i,j,k)));
%             c =  c+1;
        end
    end
end
% disp(max(err(:)))
% % 
% %%
% % (2,2,2:Nz)
% c = 1;
for i =1
    for j = 1
        for k = 2:Nz
            W = zeros(4,6);
            Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
            Qi = [xxh(i+1,j,k),yyh(i+1,j,k),zzh(i+1,j,k)];
            Qr = [xxh(i,j+1,k),yyh(i,j+1,k),zzh(i,j+1,k)];
%             Q1 = Qi;
%             Q2 = Qr;
            ti = Qi - Q;
            tr = Qr - Q;
            K1 = [Pxh(i,j,k-1),Pyh(i,j,k-1),Pzh(i,j,k-1)];
            K2 = [Pxh(i,j,k),Pyh(i,j,k),Pzh(i,j,k)];
            k_K1 = Genk_Robin(K1);
            k_K2 = Genk_Robin(K2);
            T11 = [1,0,0;0,1,0;0,0,1];
            T12 = [k_K1*cross(tr,ti)',ti',tr']/[k_K2*cross(tr,ti)',ti',tr'];
            alphaQ = 1;
            betaQ = 1;
            v1 = [-1;0;0];
            v2 = [0;-1;0];
            id = 1;
            g1 = genBoundcondition(Q,id);
            id = 3;
            g2 = genBoundcondition(Q,id);
            W(1,1:2) = 1;
            W(1,3:6) = betaQ;
            W(2:4,1) = T11*(K1-Q)';
            W(2:4,2) = T12*(K2-Q)';
            W(2:4,3) = alphaQ*k_K1'*v1;
            W(2:4,4) = alphaQ*k_K1'*v2;
            W(2:4,5) = alphaQ*T12*k_K2'*v1;
            W(2:4,6) = alphaQ*T12*k_K2'*v2;
            ww = W'*((W*W')\[1;0;0;0]);
            Ch(9,xid,yid,zid) = 0;
%             Ch(4,i,j,k) = ww(1);
%             Ch(8,i,j,k) = ww(2);
%             Ch(9,i,j,k) = ww(3)*g1+ww(4)*g2+ww(5)*g1+ww(6)*g2;
%             err(c) = abs(vh(i,j,k) - (Ch(4,i,j,k)*uh(i,j,k-1)+...
%                 Ch(8,i,j,k)*uh(i,j,k)+Ch(9,i,j,k)));
%             c =  c+1;
        end
    end
end
% disp(max(err(:)))
% %%
% % (2,Ny,2:Nz)
% c = 1;
for i = 1
    for j = Ny+1
        for k = 2:Nz
            W = zeros(4,6);
            Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
            Qi = [xxh(i+1,j,k),yyh(i+1,j,k),zzh(i+1,j,k)];
            Ql = [xxh(i,j-1,k),yyh(i,j-1,k),zzh(i,j-1,k)];
%             Q1 = Qi;
%             Q2 = Ql;
            ti = Qi - Q;
            tl = Ql - Q;
            K1 = [Pxh(i,j-1,k-1),Pyh(i,j-1,k-1),Pzh(i,j-1,k-1)];
            K2 = [Pxh(i,j-1,k),Pyh(i,j-1,k),Pzh(i,j-1,k)];
            k_K1 = Genk_Robin(K1)';
            k_K2 = Genk_Robin(K2)';
            T11 = [1,0,0;0,1,0;0,0,1];
            T12 = [k_K1*cross(tl,ti)',ti',tl']/[k_K2*cross(tl,ti)',ti',tl'];
            alphaQ = 1;
            betaQ = 1;
            v1 = [-1;0;0];
            v2 = [0;1;0];
            id = 1;
            g1 = genBoundcondition(Q,id);
            id = 4;
            g2 = genBoundcondition(Q,id);
            W(1,1:2) = 1;
            W(1,3:6) = betaQ;
            W(2:4,1) = T11*(K1-Q)';
            W(2:4,2) = T12*(K2-Q)';
            W(2:4,3) = alphaQ*k_K1'*v1;
            W(2:4,4) = alphaQ*k_K1'*v2;
            W(2:4,5) = alphaQ*T12*k_K2'*v1;
            W(2:4,6) = alphaQ*T12*k_K2'*v2;
            ww = W'*((W*W')\[1;0;0;0]);
            Ch(9,xid,yid,zid) = 0;
%             Ch(1,i,j,k) = ww(1);
%             Ch(5,i,j,k) = ww(2);
%             Ch(9,i,j,k) = ww(3)*g1+ww(4)*g2+ww(5)*g1+ww(6)*g2;
%             err(c) = abs(vh(i,j,k) - (Ch(1,i,j,k)*uh(i,j-1,k-1)+...
%                 Ch(5,i,j,k)*uh(i,j-1,k)+Ch(9,i,j,k)));
%             c =  c+1;
        end
    end
end
% disp(max(err(:)))
% %%
% % (Nx,2,2:Nz)
% c = 1;
for i = Nx+1
    for j = 1
        for k = 2:Nz
            W = zeros(4,6);
            Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
            Qo = [xxh(i-1,j,k),yyh(i-1,j,k),zzh(i-1,j,k)];
            Qr = [xxh(i,j+1,k),yyh(i,j+1,k),zzh(i,j+1,k)];
%             Q1 = Qo;
%             Q2 = Qr;
            to = Qo - Q;
            tr = Qr - Q;
            K1 = [Pxh(i-1,j,k-1),Pyh(i-1,j,k-1),Pzh(i-1,j,k-1)];
            K2 = [Pxh(i-1,j,k),Pyh(i-1,j,k),Pzh(i-1,j,k)];
            k_K1 = Genk_Robin(K1);
            k_K2 = Genk_Robin(K2);
            T11 = [1,0,0;0,1,0;0,0,1];
            T12 = [k_K1*cross(tr,to)',to',tr']/[k_K2*cross(tr,to)',to',tr'];
            alphaQ = 1;
            betaQ = 1;
            v1 = [1;0;0];
            v2 = [0;-1;0];
            id = 2;
            g1 = genBoundcondition(Q,id);
            id = 3;
            g2 = genBoundcondition(Q,id);
            W(1,1:2) = 1;
            W(1,3:6) = betaQ;
            W(2:4,1) = T11*(K1-Q)';
            W(2:4,2) = T12*(K2-Q)';
            W(2:4,3) = alphaQ*k_K1'*v1;
            W(2:4,4) = alphaQ*k_K1'*v2;
            W(2:4,5) = alphaQ*T12*k_K2'*v1;
            W(2:4,6) = alphaQ*T12*k_K2'*v2;
            ww = W'*((W*W')\[1;0;0;0]);
            Ch(3,i,j,k) = ww(1);
            Ch(7,i,j,k) = ww(2);
             Ch(9,i,j,k) = ww(3)*g1+ww(4)*g2+ww(5)*g1+ww(6)*g2;
%              err(c) = (vh(i,j,k) - (Ch(3,i,j,k)*uh(i-1,j,k-1)+...
%                 Ch(7,i,j,k)*uh(i-1,j,k)+Ch(9,i,j,k)));
%             c =  c+1;
        end
    end
end
% disp(max(err(:)))
% %%
% % (Nx,Ny,2:Nz)
% c = 1;
for i = Nx+1
    for j = Ny+1
        for k = 2:Nz
            W = zeros(4,6);
            Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
            Qo = [xxh(i-1,j,k),yyh(i-1,j,k),zzh(i-1,j,k)];
            Ql = [xxh(i,j-1,k),yyh(i,j-1,k),zzh(i,j-1,k)];
%             Q1 = Qo;
%             Q2 = Ql;
            to = Qo - Q;
            tl = Ql - Q;
            K1 = [Pxh(i-1,j-1,k-1),Pyh(i-1,j-1,k-1),Pzh(i-1,j-1,k-1)];
            K2 = [Pxh(i-1,j-1,k),Pyh(i-1,j-1,k),Pzh(i-1,j-1,k)];
            k_K1 = Genk_Robin(K1);
            k_K2 = Genk_Robin(K2);
            T11 = [1,0,0;0,1,0;0,0,1];
            T12 = [k_K1*cross(tl,to)',to',tl']/[k_K2*cross(tl,to)',to',tl'];
            alphaQ = 1;
            betaQ = 1;
            v1 = [1;0;0];
            v2 = [0;1;0];
            id = 2;
            g1 = genBoundcondition(Q,id);
            id = 4;
            g2 = genBoundcondition(Q,id);
            W(1,1:2) = 1;
            W(1,3:6) = betaQ;
            W(2:4,1) = T11*(K1-Q)';
            W(2:4,2) = T12*(K2-Q)';
            W(2:4,3) = alphaQ*k_K1'*v1;
            W(2:4,4) = alphaQ*k_K1'*v2;
            W(2:4,5) = alphaQ*T12*k_K2'*v1;
            W(2:4,6) = alphaQ*T12*k_K2'*v2;
            ww = W'*((W*W')\[1;0;0;0]);
            Ch(2,i,j,k) = ww(1);
            Ch(6,i,j,k) = ww(2);
            Ch(9,i,j,k) = ww(3)*g1+ww(4)*g2+ww(5)*g1+ww(6)*g2;
%             err(c) = (vh(i,j,k) - (Ch(2,i,j,k)*uh(i-1,j-1,k-1)+...
%                 Ch(6,i,j,k)*uh(i-1,j-1,k)+Ch(9,i,j,k)));
%             c =  c+1;
        end
    end
end
% disp(max(abs(err(:))))
% % %
% % %(2,2:Ny,2)
% c = 1;
for i = 1
    for j = 2:Ny
        for k = 1
            W = zeros(4,6);
            Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
            Qu = [xxh(i,j,k+1),yyh(i,j,k+1),zzh(i,j,k+1)];
            Qi = [xxh(i+1,j,k),yyh(i+1,j,k),zzh(i+1,j,k)];
%             Q1 = Qu;
%             Q2 = Qi;
            tu = Qu - Q;
            ti = Qi - Q;
            K1 = [Pxh(i,j-1,k),Pyh(i,j-1,k),Pzh(i,j-1,k)];
            K2 = [Pxh(i,j,k),Pyh(i,j,k),Pzh(i,j,k)];
            k_K1 = Genk_Robin(K1);
            k_K2 = Genk_Robin(K2);
            T11 = [1,0,0;0,1,0;0,0,1];
            T12 = [k_K1*cross(tu,ti)',tu',ti']/[k_K2*cross(tu,ti)',tu',ti'];
            alphaQ = 1;
            betaQ = 1;
            v1 = [-1;0;0];
            v2 = [0;0;-1];
            id = 1;
            g1 = genBoundcondition(Q,id);
            id = 5;
            g2 = genBoundcondition(Q,id);
            W(1,1:2) = 1;
            W(1,3:6) = betaQ;
            W(2:4,1) = T11*(K1-Q)';
            W(2:4,2) = T12*(K2-Q)';
            W(2:4,3) = alphaQ*k_K1'*v1;
            W(2:4,4) = alphaQ*k_K1'*v2;
            W(2:4,5) = alphaQ*T12*k_K2'*v1;
            W(2:4,6) = alphaQ*T12*k_K2'*v2;
            ww = W'*((W*W')\[1;0;0;0]);
            Ch(9,xid,yid,zid) = 0;
%             Ch(5,i,j,k) = ww(1);
%             Ch(8,i,j,k) = ww(2);
%             Ch(9,i,j,k) = ww(3)*g1+ww(4)*g2+ww(5)*g1+ww(6)*g2;
%             err(c) = (vh(i,j,k) - (Ch(5,i,j,k)*uh(i,j-1,k)+...
%                 Ch(8,i,j,k)*uh(i,j,k)+Ch(9,i,j,k)));
%             c =  c+1;
        end
    end
end
% disp(max(abs(err(:))))
% %%
% % (2,2:Ny,Nz)
% c = 1;
for i = 1
    for j = 2:Ny
        for k = Nz+1
            W = zeros(4,6);
            Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
            Qd = [xxh(i,j,k-1),yyh(i,j,k-1),zzh(i,j,k-1)];
            Qi = [xxh(i+1,j,k),yyh(i+1,j,k),zzh(i+1,j,k)];
%             Q1 = Qd;
%             Q2 = Qi;
            td = Qd - Q;
            ti = Qi - Q;
            K1 = [Pxh(i,j-1,k-1),Pyh(i,j-1,k-1),Pzh(i,j-1,k-1)];
            K2 = [Pxh(i,j,k-1),Pyh(i,j,k-1),Pzh(i,j,k-1)];
            k_K1 = Genk_Robin(K1);
            k_K2 = Genk_Robin(K2);
            T11 = [1,0,0;0,1,0;0,0,1];
            T12 = [k_K1*cross(td,ti)',td',ti']/[k_K2*cross(td,ti)',td',ti'];
            alphaQ = 1;
            betaQ = 1;
            v1 = [-1;0;0];
            v2 = [0;0;1];
            id = 1;
            g1 = genBoundcondition(Q,id);
            id = 6;
            g2 = genBoundcondition(Q,id);
            W(1,1:2) = 1;
            W(1,3:6) = betaQ;
            W(2:4,1) = T11*(K1-Q)';
            W(2:4,2) = T12*(K2-Q)';
            W(2:4,3) = alphaQ*k_K1'*v1;
            W(2:4,4) = alphaQ*k_K1'*v2;
            W(2:4,5) = alphaQ*T12*k_K2'*v1;
            W(2:4,6) = alphaQ*T12*k_K2'*v2;
            ww = W'*((W*W')\[1;0;0;0]);
            Ch(9,xid,yid,zid) = 0;
%             Ch(1,i,j,k) = ww(1);
%             Ch(4,i,j,k) = ww(2);
%             Ch(9,i,j,k) = ww(3)*g1+ww(4)*g2+ww(5)*g1+ww(6)*g2;
%             err(c) = (vh(i,j,k) - (Ch(1,i,j,k)*uh(i,j-1,k-1)+...
%                 Ch(4,i,j,k)*uh(i,j,k-1)+Ch(9,i,j,k)));
%             c =  c+1;
        end
    end
end
% disp(max(abs(err(:))))
% %%
% (Nx,2:Ny,2)
% c = 1;
for i = Nx+1
    for j = 2:Ny
        for k = 1
            W = zeros(4,6);
            Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
            Qu = [xxh(i,j,k+1),yyh(i,j,k+1),zzh(i,j,k+1)];
            Qo = [xxh(i-1,j,k),yyh(i-1,j,k),zzh(i-1,j,k)];
%             Q1 = Qu;
%             Q2 = Qo;
            tu = Qu - Q;
            to = Qo - Q;
            K1 = [Pxh(i-1,j-1,k),Pyh(i-1,j-1,k),Pzh(i-1,j-1,k)];
            K2 = [Pxh(i-1,j,k),Pyh(i-1,j,k),Pzh(i-1,j,k)];
            k_K1 = Genk_Robin(K1)';
            k_K2 = Genk_Robin(K2)';
            T11 = [1,0,0;0,1,0;0,0,1];
            T12 = [k_K1*cross(tu,to)',tu',to']/[k_K2*cross(tu,to)',tu',to'];
            alphaQ = 1;
            betaQ = 1;
            v1 = [1;0;0];
            v2 = [0;0;-1];
            id = 2;
            g1 = genBoundcondition(Q,id);
            id = 5;
            g2 = genBoundcondition(Q,id);
            W(1,1:2) = 1;
            W(1,3:6) = betaQ;
            W(2:4,1) = T11*(K1-Q)';
            W(2:4,2) = T12*(K2-Q)';
            W(2:4,3) = alphaQ*k_K1'*v1;
            W(2:4,4) = alphaQ*k_K1'*v2;
            W(2:4,5) = alphaQ*T12*k_K2'*v1;
            W(2:4,6) = alphaQ*T12*k_K2'*v2;
            ww = W'*((W*W')\[1;0;0;0]);
            Ch(6,i,j,k) = ww(1);
            Ch(7,i,j,k) = ww(2);
            Ch(9,i,j,k) = ww(3)*g1+ww(4)*g2+ww(5)*g1+ww(6)*g2;
%             err(c) = (vh(i,j,k) - (Ch(6,i,j,k)*uh(i-1,j-1,k)+...
%                 Ch(7,i,j,k)*uh(i-1,j,k)+Ch(9,i,j,k)));
%             c =  c+1;
        end
    end
end
% disp(max(abs(err(:))))
%%
% %(Nx,2:Ny,Nz)
% c = 1;
for i =Nx+1
    for j = 2:Ny
        for k = Nz+1
            W = zeros(4,6);
            Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
            Qd = [xxh(i,j,k-1),yyh(i,j,k-1),zzh(i,j,k-1)];
            Qo = [xxh(i-1,j,k),yyh(i-1,j,k),zzh(i-1,j,k)];
%             Q1 = Qd;
%             Q2 = Qo;
            td = Qd - Q;
            to = Qo - Q;
            K1 = [Pxh(i-1,j-1,k-1),Pyh(i-1,j-1,k-1),Pzh(i-1,j-1,k-1)];
            K2 = [Pxh(i-1,j,k-1),Pyh(i-1,j,k-1),Pzh(i-1,j,k-1)];
            k_K1 = Genk_Robin(K1)';
            k_K2 = Genk_Robin(K2)';
            T11 = [1,0,0;0,1,0;0,0,1];
            T12 = [k_K1*cross(td,to)',td',to']/[k_K2*cross(td,to)',td',to'];
            alphaQ = 1;
            betaQ = 1;
            v1 = [1;0;0];
            v2 = [0;0;1];
            id = 2;
            g1 = genBoundcondition(Q,id);
            id = 6;
            g2 = genBoundcondition(Q,id);
            W(1,1:2) = 1;
            W(1,3:6) = betaQ;
            W(2:4,1) = T11*(K1-Q)';
            W(2:4,2) = T12*(K2-Q)';
            W(2:4,3) = alphaQ*k_K1'*v1;
            W(2:4,4) = alphaQ*k_K1'*v2;
            W(2:4,5) = alphaQ*T12*k_K2'*v1;
            W(2:4,6) = alphaQ*T12*k_K2'*v2;
            ww = W'*((W*W')\[1;0;0;0]);
            Ch(2,i,j,k) = ww(1);
            Ch(3,i,j,k) = ww(2);
             Ch(9,i,j,k) = ww(3)*g1+ww(4)*g2+ww(5)*g1+ww(6)*g2;
%             err(c) = (vh(i,j,k) - (Ch(2,i,j,k)*uh(i-1,j-1,k-1)+...
%                 Ch(3,i,j,k)*uh(i-1,j,k-1)+Ch(9,i,j,k)));
%             c =  c+1;
         end
    end
 end

% disp(max(abs(err(:))))
% %
% 面
%%
% %(1,,)
% c = 1;
for i = 1
    for j = 2:Ny
        for k = 2:Nz
            W = zeros(4,8);
            Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
            Qu = [xxh(i,j,k+1),yyh(i,j,k+1),zzh(i,j,k+1)];
            Qd = [xxh(i,j,k-1),yyh(i,j,k-1),zzh(i,j,k-1)];
            Ql = [xxh(i,j-1,k),yyh(i,j-1,k),zzh(i,j-1,k)];
%             Qr = [xxh(i,j+1,k),yyh(i,j+1,k),zzh(i,j+1,k)];
            Qi = [xxh(i+1,j,k),yyh(i+1,j,k),zzh(i+1,j,k)];
            tu = Qu-Q;
            td = Qd-Q;
            tl = Ql-Q;
%             tr = Qr-Q;
            ti = Qi-Q;
            Kul = [Pxh(i,j-1,k),Pyh(i,j-1,k),Pzh(i,j-1,k)];
            Kur = [Pxh(i,j,k),Pyh(i,j,k),Pzh(i,j,k)];
            Kdl = [Pxh(i,j-1,k-1),Pyh(i,j-1,k-1),Pzh(i,j-1,k-1)];
            Kdr = [Pxh(i,j,k-1),Pyh(i,j,k-1),Pzh(i,j,k-1)];
            k_Kul = Genk_Robin(Kul);
            k_Kur = Genk_Robin(Kur);
            k_Kdl = Genk_Robin(Kdl);
            k_Kdr = Genk_Robin(Kdr);
            Tul = [1,0,0;0,1,0;0,0,1];
            Tur = [k_Kul*cross(tu,ti)',tu',ti']/[k_Kur*cross(tu,ti)',tu',ti'];
            Tdl = [k_Kul*cross(tl,ti)',tl',ti']/[k_Kdl*cross(tl,ti)',tl',ti'];
            Tdr =Tdl*[k_Kdl*cross(td,ti)',td',ti']/[k_Kdr*cross(td,ti)',td',ti'];
            id = 1;
            g = genBoundcondition(Q,id);
            v = [-1;0;0];
            W(1,1:4) = 1;
            W(1,5:8) = 1;
            W(2:4,1) = Tul*(Kul-Q)';
            W(2:4,2) = Tur*(Kur-Q)';
            W(2:4,3) = Tdl*(Kdl-Q)';
            W(2:4,4) = Tdr*(Kdr-Q)';
            W(2:4,5) = Tul*k_Kul*v;
            W(2:4,6) = Tur*k_Kur*v;
            W(2:4,7) = Tdl*k_Kdl*v;
            W(2:4,8) = Tdr*k_Kdr*v;
            ww = W'*((W*W')\[1;0;0;0]);
%             Ch(5,i,j,k) = ww(1);
%             Ch(8,i,j,k) = ww(2);
%             Ch(1,i,j,k) = ww(3);
%             Ch(4,i,j,k) = ww(4);
%             Ch(9,i,j,k) = sum(ww(5:8))*g;
            Ch(9,i,j,k) = 0;
%              err(c) = abs(vh(i,j,k) - (Ch(8,i,j,k)*uh(i,j,k)+Ch(5,i,j,k)*uh(i,j-1,k)...
%                  +Ch(1,i,j,k)*uh(i,j-1,k-1)+Ch(4,i,j,k)*uh(i,j,k-1)+Ch(9,i,j,k)));
%              c = c+1;
        end
    end
end
%     disp(max(abs(err(:))))

%%
% (Nx+1,,)
% err = 0;
% c = 1;
for i = Nx+1
   for j = 2:Ny
       for k = 2:Nz
           W = zeros(4,8);
           Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
           Qu = [xxh(i,j,k+1),yyh(i,j,k+1),zzh(i,j,k+1)];
           Qd = [xxh(i,j,k-1),yyh(i,j,k-1),zzh(i,j,k-1)];
           Ql = [xxh(i,j-1,k),yyh(i,j-1,k),zzh(i,j-1,k)];
%            Qr = [xxh(i,j+1,k),yyh(i,j+1,k),zzh(i,j+1,k)];
           Qo = [xxh(i-1,j,k),yyh(i-1,j,k),zzh(i-1,j,k)];
           tu = Qu-Q;
           td = Qd-Q;
           tl = Ql-Q;
%            tr = Qr-Q;
           to = Qo-Q;
           Kul = [Pxh(i-1,j-1,k),Pyh(i-1,j-1,k),Pzh(i-1,j-1,k)];
           Kur = [Pxh(i-1,j,k),Pyh(i-1,j,k),Pzh(i-1,j,k)];
           Kdl = [Pxh(i-1,j-1,k-1),Pyh(i-1,j-1,k-1),Pzh(i-1,j-1,k-1)];
           Kdr = [Pxh(i-1,j,k-1),Pyh(i-1,j,k-1),Pzh(i-1,j,k-1)];
           k_Kul = Genk_Robin(Kul);
           k_Kur = Genk_Robin(Kur);
           k_Kdl = Genk_Robin(Kdl);
           k_Kdr = Genk_Robin(Kdr);
           Tul = [1,0,0;0,1,0;0,0,1];
           Tur = [k_Kul*cross(tu,to)',tu',to']/[k_Kur*cross(tu,to)',tu',to'];
           Tdl = [k_Kul*cross(tl,to)',tl',to']/[k_Kdl*cross(tl,to)',tl',to'];
           Tdr =Tdl*[k_Kdl*cross(td,to)',td',to']/[k_Kdr*cross(td,to)',td',to'];
           id = 2;
           g = genBoundcondition(Q,id);
           v = [1;0;0];
           alphaQ = 1;
           betaQ = 1;
           W(1,1:4) = 1;
           W(1,5:8) = betaQ;
           W(2:4,1) = Tul*(Kul-Q)';
           W(2:4,2) = Tur*(Kur-Q)';
           W(2:4,3) = Tdl*(Kdl-Q)';
           W(2:4,4) = Tdr*(Kdr-Q)';
           W(2:4,5) = alphaQ*Tul*k_Kul*v;
           W(2:4,6) = alphaQ*Tur*k_Kur*v;
           W(2:4,7) = alphaQ*Tdl*k_Kdl*v;
           W(2:4,8) = alphaQ*Tdr*k_Kdr*v;
           ww = W'*((W*W')\[1;0;0;0]);
           Ch(6,i,j,k) = ww(1);
           Ch(7,i,j,k) = ww(2);
           Ch(2,i,j,k) = ww(3);
           Ch(3,i,j,k) = ww(4);
           Ch(9,i,j,k) = sum(ww(5:8))*g;
%            err(c) = abs(vh(i,j,k) - (Ch(7,i,j,k)*uh(i-1,j,k)+Ch(6,i,j,k)*uh(i-1,j-1,k)...
%                  +Ch(2,i,j,k)*uh(i-1,j-1,k-1)+Ch(3,i,j,k)*uh(i-1,j,k-1)+Ch(9,i,j,k)));
%             c=c+1;
       end    
   end
end
% disp(max(abs(err(:))))
% %%
% % (,1,)
% % % err = 0;
% c = 1 ;
for i = 2:Nx
   for j = 1
       for k = 2:Nx
           W = zeros(4,8);
           Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
           Qu = [xxh(i,j,k+1),yyh(i,j,k+1),zzh(i,j,k+1)];
           Qd = [xxh(i,j,k-1),yyh(i,j,k-1),zzh(i,j,k-1)];
%            Qi = [xxh(i+1,j,k),yyh(i+1,j,k),zzh(i+1,j,k)];
           Qo = [xxh(i-1,j,k),yyh(i-1,j,k),zzh(i-1,j,k)];
           Qr = [xxh(i,j+1,k),yyh(i,j+1,k),zzh(i,j+1,k)];
           tu = Qu-Q;
           td = Qd-Q;
%            ti = Qi-Q;
           to = Qo-Q;
           tr = Qr-Q;
           Kuo = [Pxh(i-1,j,k),Pyh(i-1,j,k),Pzh(i-1,j,k)];
           Kui = [Pxh(i,j,k),Pyh(i,j,k),Pzh(i,j,k)];
           Kdo = [Pxh(i-1,j,k-1),Pyh(i-1,j,k-1),Pzh(i-1,j,k-1)];
           Kdi = [Pxh(i,j,k-1),Pyh(i,j,k-1),Pzh(i,j,k-1)];
           k_Kuo = Genk_Robin(Kuo);
           k_Kui = Genk_Robin(Kui);
           k_Kdo = Genk_Robin(Kdo);
           k_Kdi = Genk_Robin(Kdi);
           Tuo = [1,0,0;0,1,0;0,0,1];
           Tui = [k_Kuo*cross(tu,tr)',tu',tr']/[k_Kui*cross(tu,tr)',tu',tr'];
           Tdo = [k_Kuo*cross(tr,to)',tr',to']/[k_Kdo*cross(tr,to)',tr',to'];
           Tdi =Tdo*[k_Kdo*cross(td,tr)',td',tr']/[k_Kdi*cross(td,tr)',td',tr'];
           id = 3;
           g = genBoundcondition(Q,id);
           v = [0;-1;0];
           alphaQ = 1;
           betaQ = 1;
           W(1,1:4) = 1;
           W(1,5:8) = betaQ;
           W(2:4,1) = Tuo*(Kuo-Q)';
           W(2:4,2) = Tui*(Kui-Q)';
           W(2:4,3) = Tdo*(Kdo-Q)';
           W(2:4,4) = Tdi*(Kdi-Q)';
           W(2:4,5) = alphaQ*Tuo*k_Kuo*v;
           W(2:4,6) = alphaQ*Tui*k_Kui*v;
           W(2:4,7) = alphaQ*Tdo*k_Kdo*v;
           W(2:4,8) = alphaQ*Tdi*k_Kdi*v;
           ww = W'*((W*W')\[1;0;0;0]);
           Ch(7,i,j,k) = ww(1);
           Ch(8,i,j,k) = ww(2);
           Ch(3,i,j,k) = ww(3);
           Ch(4,i,j,k) = ww(4);
           Ch(9,i,j,k) = sum(ww(5:8))*g;
%            err(c) = abs(vh(i,j,k) - (Ch(8,i,j,k)*uh(i,j,k)+Ch(7,i,j,k)*uh(i-1,j,k)...
%                +Ch(3,i,j,k)*uh(i-1,j,k-1)+Ch(4,i,j,k)*uh(i,j,k-1)+Ch(9,i,j,k)));
%           c = c+1;
       end
    end
end
% disp(max(abs(err(:))))
% % 
% %%
% % (,Ny+1,)
% c = 1;
for i = 2:Nx
   for j = Ny+1
       for k = 2:Nz
           W = zeros(4,8);
           Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
           Qu = [xxh(i,j,k+1),yyh(i,j,k+1),zzh(i,j,k+1)];
           Qd = [xxh(i,j,k-1),yyh(i,j,k-1),zzh(i,j,k-1)];
%            Qi = [xxh(i+1,j,k),yyh(i+1,j,k),zzh(i+1,j,k)];
           Qo = [xxh(i-1,j,k),yyh(i-1,j,k),zzh(i-1,j,k)];
           Ql = [xxh(i,j-1,k),yyh(i,j-1,k),zzh(i,j-1,k)];
           tu = Qu-Q;
           td = Qd-Q;
%            ti = Qi-Q;
           to = Qo-Q;
           tl = Ql-Q;
           Kuo = [Pxh(i-1,j-1,k),Pyh(i-1,j-1,k),Pzh(i-1,j-1,k)];
           Kui = [Pxh(i,j-1,k),Pyh(i,j-1,k),Pzh(i,j-1,k)];
           Kdo = [Pxh(i-1,j-1,k-1),Pyh(i-1,j-1,k-1),Pzh(i-1,j-1,k-1)];
           Kdi = [Pxh(i,j-1,k-1),Pyh(i,j-1,k-1),Pzh(i,j-1,k-1)];
           k_Kuo = Genk_Robin(Kuo);
           k_Kui = Genk_Robin(Kui);
           k_Kdo = Genk_Robin(Kdo);
           k_Kdi = Genk_Robin(Kdi);
           Tuo = [1,0,0;0,1,0;0,0,1];
           Tui = [k_Kuo*cross(tu,tl)',tu',tl']/[k_Kui*cross(tu,tl)',tu',tl'];
           Tdo = [k_Kuo*cross(tl,to)',tl',to']/[k_Kdo*cross(tl,to)',tl',to'];
           Tdi =Tdo*[k_Kdo*cross(td,tl)',td',tl']/[k_Kdi*cross(td,tl)',td',tl'];
           v = [0;1;0];
           id = 4;
           g = genBoundcondition(Q,id);
           alphaQ = 1;
           betaQ = 1;
           W(1,1:4) = 1;
           W(1,5:8) = betaQ;
           W(2:4,1) = Tuo*(Kuo-Q)';
           W(2:4,2) = Tui*(Kui-Q)';
           W(2:4,3) = Tdo*(Kdo-Q)';
           W(2:4,4) = Tdi*(Kdi-Q)';
           W(2:4,5) = alphaQ*Tuo*k_Kuo*v;
           W(2:4,6) = alphaQ*Tui*k_Kui*v;
           W(2:4,7) = alphaQ*Tdo*k_Kdo*v;
           W(2:4,8) = alphaQ*Tdi*k_Kdi*v;
           ww = W'*((W*W')\[1;0;0;0]);
           Ch(5,i,j,k) = ww(2);
           Ch(6,i,j,k) = ww(1);
           Ch(1,i,j,k) = ww(4);
           Ch(2,i,j,k) = ww(3);
           Ch(9,i,j,k) = sum(ww(5:8))*g;
%            err(c) = abs(vh(i,j,k) - (Ch(5,i,j,k)*uh(i,j-1,k)+Ch(6,i,j,k)*uh(i-1,j-1,k)...
%                +Ch(1,i,j,k)*uh(i,j-1,k-1)+Ch(2,i,j,k)*uh(i-1,j-1,k-1)+Ch(9,i,j,k)));
%            c = c+1;
       end
    end
end
% disp(max(abs(err(:))))
% 
%%
% (,,1)
% c = 1;
for i = 2:Nx
    for j = 2:Ny
        for k = 1
            W = zeros(4,8);
            Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
            Qi = [xxh(i+1,j,k),yyh(i+1,j,k),zzh(i+1,j,k)];
            Qo = [xxh(i-1,j,k),yyh(i-1,j,k),zzh(i-1,j,k)];
            Ql = [xxh(i,j-1,k),yyh(i,j-1,k),zzh(i,j-1,k)];
%             Qr = [xxh(i,j+1,k),yyh(i,j+1,k),zzh(i,j+1,k)];
            Qu = [xxh(i,j,k+1),yyh(i,j,k+1),zzh(i,j,k+1)];
            ti = Qi-Q;
            to = Qo-Q;
            tl = Ql-Q;
%             tr = Qr-Q;
            tu = Qu-Q;
            Kil = [Pxh(i,j-1,k),Pyh(i,j-1,k),Pzh(i,j-1,k)];
            Kir = [Pxh(i,j,k),Pyh(i,j,k),Pzh(i,j,k)];
            Kol = [Pxh(i-1,j-1,k),Pyh(i-1,j-1,k),Pzh(i-1,j-1,k)];
            Kor = [Pxh(i-1,j,k),Pyh(i-1,j,k),Pzh(i-1,j,k)];
            k_Kil = Genk_Robin(Kil);
            k_Kir = Genk_Robin(Kir);
            k_Kol = Genk_Robin(Kol);
            k_Kor = Genk_Robin(Kor);
            Til = [1,0,0;0,1,0;0,0,1];
            Tir = [k_Kil*cross(tu,ti)',tu',ti']/[k_Kir*cross(tu,ti)',tu',ti'];
            Tol = [k_Kil*cross(tl,tu)',tl',tu']/[k_Kol*cross(tl,tu)',tl',tu'];
            Tor = Tol*[k_Kol*cross(to,tu)',to',tu']/[k_Kor*cross(to,tu)',to',tu'];
            v = [0;0;-1];
            alphaQ = 1;
            betaQ = 1;
            id = 5;
            g = genBoundcondition(Q,id);
            W(1,1:4) = 1;
            W(1,5:8) = betaQ;
            W(2:4,1) = Til*(Kil-Q)';
            W(2:4,2) = Tir*(Kir-Q)';
            W(2:4,3) = Tol*(Kol-Q)';
            W(2:4,4) = Tor*(Kor-Q)';
            W(2:4,5) = alphaQ*Til*k_Kil*v;
            W(2:4,6) = alphaQ*Tir*k_Kir*v;
            W(2:4,7) = alphaQ*Tol*k_Kol*v;
            W(2:4,8) = alphaQ*Tor*k_Kor*v;
            ww = W'*((W*W')\[1;0;0;0]);
            Ch(5,i,j,k) = ww(1);
            Ch(8,i,j,k) = ww(2);
            Ch(6,i,j,k) = ww(3);
            Ch(7,i,j,k) = ww(4);
            Ch(9,i,j,k) = sum(ww(5:8))*g;
%             err(c) = abs(vh(i,j,k) - (Ch(8,i,j,k)*uh(i,j,k)+Ch(5,i,j,k)*uh(i,j-1,k)...
%                +Ch(6,i,j,k)*uh(i-1,j-1,k)+Ch(7,i,j,k)*uh(i-1,j,k)+Ch(9,i,j,k)));
%             c = c+1;
        end
    end
end
% disp(max(abs(err(:))))
%%
% (,,Nz+1)
% c = 1;
for i = 2:Nx
    for j = 2:Ny
        for k = Nz+1
            W = zeros(4,8);
            Q = [xxh(i,j,k),yyh(i,j,k),zzh(i,j,k)];
            Qi = [xxh(i+1,j,k),yyh(i+1,j,k),zzh(i+1,j,k)];
            Qo = [xxh(i-1,j,k),yyh(i-1,j,k),zzh(i-1,j,k)];
            Ql = [xxh(i,j-1,k),yyh(i,j-1,k),zzh(i,j-1,k)];
%             Qr = [xxh(i,j+1,k),yyh(i,j+1,k),zzh(i,j+1,k)];
            Qd = [xxh(i,j,k-1),yyh(i,j,k-1),zzh(i,j,k-1)];
            ti = Qi-Q;
            to = Qo-Q;
            tl = Ql-Q;
%             tr = Qr-Q;
            td = Qd-Q;
            Kil = [Pxh(i,j-1,k-1),Pyh(i,j-1,k-1),Pzh(i,j-1,k-1)];
            Kir = [Pxh(i,j,k-1),Pyh(i,j,k-1),Pzh(i,j,k-1)];
            Kol = [Pxh(i-1,j-1,k-1),Pyh(i-1,j-1,k-1),Pzh(i-1,j-1,k-1)];
            Kor = [Pxh(i-1,j,k-1),Pyh(i-1,j,k-1),Pzh(i-1,j,k-1)];
            k_Kil = Genk_Robin(Kil);
            k_Kir = Genk_Robin(Kir);
            k_Kol = Genk_Robin(Kol);
            k_Kor = Genk_Robin(Kor);
            Til = [1,0,0;0,1,0;0,0,1];
            Tir = [k_Kil*cross(td,ti)',td',ti']/[k_Kir*cross(td,ti)',td',ti'];
            Tol = [k_Kil*cross(tl,td)',tl',td']/[k_Kol*cross(tl,td)',tl',td'];
            Tor = Tol*[k_Kol*cross(to,td)',to',td']/[k_Kor*cross(to,td)',to',td'];
            v = [0;0;1];
            alphaQ = 1;
            betaQ = 1;
            id = 6;
            g = genBoundcondition(Q,id);
            W(1,1:4) = 1;
            W(1,5:8) = betaQ;
            W(2:4,1) = Til*(Kil-Q)';
            W(2:4,2) = Tir*(Kir-Q)';
            W(2:4,3) = Tol*(Kol-Q)';
            W(2:4,4) = Tor*(Kor-Q)';
            W(2:4,5) = alphaQ*Til*k_Kil*v;
            W(2:4,6) = alphaQ*Tir*k_Kir*v;
            W(2:4,7) = alphaQ*Tol*k_Kol*v;
            W(2:4,8) = alphaQ*Tor*k_Kor*v;
            ww = W'*((W*W')\[1;0;0;0]);
            Ch(1,i,j,k) = ww(1);
            Ch(4,i,j,k) = ww(2);
            Ch(2,i,j,k) = ww(3);
            Ch(3,i,j,k) = ww(4);
            Ch(9,i,j,k) = sum(ww(5:8))*g;
%              err(c) = abs(vh(i,j,k) - (Ch(4,i,j,k)*uh(i,j,k-1)+Ch(1,i,j,k)*uh(i,j-1,k-1)...
%                +Ch(2,i,j,k)*uh(i-1,j-1,k-1)+Ch(3,i,j,k)*uh(i-1,j,k-1)+Ch(9,i,j,k)));
%             c = c+1;
        end
    end
end
% disp(max(abs(err(:))))
% 
%%
% %内部节点
% %%内部单元节点
% Vh = zeros(Nx-1,Ny-1,Nz-1);

for i = 1:Nx-1
    for j = 1:Ny-1
        for k = 1:Nz-1
       
            W = zeros(4,8);
            W(1,:) = 1;
            A = [xxh(i+1,j+1,k+1),yyh(i+1,j+1,k+1),zzh(i+1,j+1,k+1)];
            
%             tu = [xxh(i+1,j+1,k+2),yyh(i+1,j+1,k+2),zzh(i+1,j+1,k+2)]-A;
            td = [xxh(i+1,j+1,k),yyh(i+1,j+1,k),zzh(i+1,j+1,k)]-A;
            tf = [xxh(i+1,j,k+1),yyh(i+1,j,k+1),zzh(i+1,j,k+1)]-A;
            tb = [xxh(i+1,j+2,k+1),yyh(i+1,j+2,k+1),zzh(i+1,j+2,k+1)]-A;
            tl = [xxh(i,j+1,k+1),yyh(i,j+1,k+1),zzh(i,j+1,k+1)]-A;
            tr = [xxh(i+2,j+1,k+1),yyh(i+2,j+1,k+1),zzh(i+2,j+1,k+1)]-A;
            %过A点的6个切向量
            
            da = [Pxh(i,j,k),Pyh(i,j,k),Pzh(i,j,k)];
            db = [Pxh(i+1,j,k),Pyh(i+1,j,k),Pzh(i+1,j,k)];
            dc = [Pxh(i,j+1,k),Pyh(i,j+1,k),Pzh(i,j+1,k)];
            dd = [Pxh(i+1,j+1,k),Pyh(i+1,j+1,k),Pzh(i+1,j+1,k)];
            ua = [Pxh(i,j,k+1),Pyh(i,j,k+1),Pzh(i,j,k+1)];
            ub = [Pxh(i+1,j,k+1),Pyh(i+1,j,k+1),Pzh(i+1,j,k+1)];
            uc = [Pxh(i,j+1,k+1),Pyh(i,j+1,k+1),Pzh(i,j+1,k+1)];
            ud = [Pxh(i+1,j+1,k+1),Pyh(i+1,j+1,k+1),Pzh(i+1,j+1,k+1)];
            %8个绕A的单元中心
            
            kda = Genk_Robin(da);
            kdb = Genk_Robin(db);
            kdc = Genk_Robin(dc);
            kdd = Genk_Robin(dd);
            
            kua = Genk_Robin(ua);
            kub = Genk_Robin(ub);
            kuc = Genk_Robin(uc);
            kud = Genk_Robin(ud);
            
            %Dda*[kda*(ndadb),td,tf] = Ddb*[kdb*(ndadb),td,tf]
            Tdadb = [kdb*cross(td,tf)',td',tf']/[kda*cross(td,tf)',td',tf'];
            %Dda*[kda*(ndadc),td,tl] = Ddc*[kdc*(ndadc),td,tl]
            Tdadc = [kdc*cross(td,tl)',td',tl']/[kda*cross(td,tl)',td',tl'];
            %Ddd*[kdd*(ndddc),td,tb] = Ddc*[kdc*(ndddc),td,tb]
            Tdddc = [kdc*cross(td,tb)',td',tb']/[kdd*cross(td,tb)',td',tb'];
            
            
            %Dda*[kda*(ndaua),tf,tl] = Dua*[kua*(ndaua),tf,tl]
            Tdaua = [kua*cross(tf,tl)',tf',tl']/[kda*cross(tf,tl)',tf',tl'];
            
            %Ddb*[kdb*(ndbub),tf,tr] = Dub*[kub*(ndbub),tf,tr]
            Tdbub = [kub*cross(tf,tr)',tf',tr']/[kdb*cross(tf,tr)',tf',tr'];
            
            %Ddc*[kdc*(ndcuc),tb,tl] = Duc*[kuc*(ndcuc),tb,tl]
            Tdcuc = [kuc*cross(tb,tl)',tb',tl']/[kdc*cross(tb,tl)',tb',tl'];
            
            %Ddd*[kdd*(nddud),tb,tr] = Dud*[kud*(nddud),tb,tr]
            Tddud = [kud*cross(tb,tr)',tb',tr']/[kdd*cross(tb,tr)',tb',tr'];
            
            
            
            
            W(2:4,1) = (da-A)';
            W(2:4,2) = Tdadb\(db-A)';
            W(2:4,3) = Tdadc\(dc-A)';
            W(2:4,4) = Tdddc*((Tdadc)\(dd-A)');
            W(2:4,5) = Tdaua\(ua-A)';
            W(2:4,6) = (Tdbub)\Tdadb\(ub-A)';
            W(2:4,7) = Tdcuc\Tdadc\(uc-A)';
            W(2:4,8) = (Tddud)\Tdddc*((Tdadc)\(ud-A)');
            
%             W(5,:) = W
            
%             W(2:4,1) = (da-A)';
%             W(2:4,2) = (db-A)';
%             W(2:4,3) = (dc-A)';
%             W(2:4,4) = (dd-A)';
%             W(2:4,5) = (ua-A)';
%             W(2:4,6) = (ub-A)';
%             W(2:4,7) = (uc-A)';
%             W(2:4,8) = (ud-A)';
            
            ww = W'*((W*W')\[1;0;0;0]);
            
            
           
%                 Up1 = uh(i+1,j,k);
%             
%                 Up2 = uh(i,j,k);
%            
%                 Up3 = uh(i,j+1,k);
%           
%                 Up4 = uh(i+1,j+1,k);
%          
%                 Up5 = uh(i+1,j,k+1);
%         
%                 Up6 = uh(i,j,k+1);
%             
%                 Up7 = uh(i,j+1,k+1);
%             
%             
%                 Up8 = uh(i+1,j+1,k+1);
            
         %   Ch(1:8,i+1,j+1,k+1) = ww;
            Ch(1,i+1,j+1,k+1) = ww(2);
            Ch(2,i+1,j+1,k+1) = ww(1);
            Ch(3,i+1,j+1,k+1) = ww(3);
            Ch(4,i+1,j+1,k+1) = ww(4);
            Ch(5,i+1,j+1,k+1) = ww(6);
            Ch(6,i+1,j+1,k+1) = ww(5);
            Ch(7,i+1,j+1,k+1) = ww(7);
            Ch(8,i+1,j+1,k+1) = ww(8);
     
%             Vh(i,j,k) = Ch(9,i+1,j+1,k+1)+Ch(1,i+1,j+1,k+1)*Up1+Ch(2,i+1,j+1,k+1)*Up2+Ch(3,i+1,j+1,k+1)*Up3+...
%                 Ch(4,i+1,j+1,k+1)*Up4+Ch(5,i+1,j+1,k+1)*Up5+Ch(6,i+1,j+1,k+1)*Up6+Ch(7,i+1,j+1,k+1)*Up7+Ch(8,i+1,j+1,k+1)*Up8;
            
        end
    end
end
% vvh = vh(2:Nx,2:Ny,2:Nz);
% disp(max(abs(vvh(:)-Vh(:))))
% % 
% % 

% for i = 1:Nx+1
%     for j = 1:Ny+1
%         for k = 1:Nz+1
%             if Ch(1,i,j,k) ~= 0
%                 Up1 = uh(i,j-1,k-1);
%             else
%                 Up1 = 0;
%             end
%             if Ch(2,i,j,k) ~= 0
%                 Up2 = uh(i-1,j-1,k-1);
%             else
%                 Up2 = 0;
%             end
%             if Ch(3,i,j,k) ~= 0
%                 Up3 = uh(i-1,j,k-1);
%             else
%                 Up3 = 0;
%             end
%             if Ch(4,i,j,k) ~= 0
%                 Up4 = uh(i,j,k-1);
%             else
%                 Up4 = 0;
%             end
%             if Ch(5,i,j,k) ~= 0
%                 Up5 = uh(i,j-1,k);
%             else
%                 Up5 = 0;
%             end
%             if Ch(6,i,j,k) ~= 0
%                 Up6 = uh(i-1,j-1,k);
%             else
%                 Up6 = 0;
%             end
%             if Ch(7,i,j,k) ~= 0
%                 Up7 = uh(i-1,j,k);
%             else
%                 Up7 = 0;
%             end
%             if Ch(8,i,j,k) ~= 0
%                 Up8 = uh(i,j,k);
%             else
%                 Up8 = 0;
%             end
%             Vh(i,j,k) = Ch(9,i,j,k)+Ch(1,i,j,k)*Up1+Ch(2,i,j,k)*Up2+Ch(3,i,j,k)*Up3+...
%                 Ch(4,i,j,k)*Up4+Ch(5,i,j,k)*Up5+Ch(6,i,j,k)*Up6+Ch(7,i,j,k)*Up7+Ch(8,i,j,k)*Up8;
%             
%         end
%     end
% end
% disp(max(abs((Vh(:)-vh(:)))))
% end
% % % 
