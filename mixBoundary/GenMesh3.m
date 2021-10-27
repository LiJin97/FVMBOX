function [xx,yy,zz] = GenMesh3(Nx,Ny,Nz,delta)
% clear 
% clc
% Nx = 10;
% Ny = 10;
% Nz = 10;
% delta = 0.2;
xx = zeros(Nx+1,Ny+1,Nz+1);
yy = xx;
zz = yy;

for xid = 1:Nx+1
    for yid = 1:Ny+1
        for zid = 1:Nz+1
            if xid==1 || xid==(Nx+1) || yid==1 || yid==(Ny+1) ||zid==1 || zid==(Nz+1)
            xx(xid,yid,zid) = (xid-1)/Nx;
            yy(xid,yid,zid) = (yid-1)/Ny;
            zz(xid,yid,zid) = (zid-1)/Nz;
            else
            xx(xid,yid,zid) = (xid-1)/Nx+delta*(rand-0.5)/Nx;
            yy(xid,yid,zid) = (yid-1)/Ny+delta*(rand-0.5)/Ny;
            zz(xid,yid,zid) = (zid-1)/Nz;%+delta*(rand-0.5)/Nz;
            end
        end
    end
end

% xx(Nx/2+1,:,:) = 0.5;

% %%
% %out
% mm = 1;
% surf(xx(:,:,mm),yy(:,:,mm),zz(:,:,mm),zz(:,:,mm));hold on
% 
% mm = 1;
% surf(reshape(xx(:,mm,:),Ny+1,Nz+1),reshape(yy(:,mm,:),Ny+1,Nz+1),reshape(zz(:,mm,:),Ny+1,Nz+1));hold on
% 
% mm = 1;
% surf(reshape(xx(mm,:,:),Nx+1,Nz+1),reshape(yy(mm,:,:),Nx+1,Nz+1),reshape(zz(mm,:,:),Nx+1,Nz+1));hold on
% 
% %%
% %in
% mm = 3;
% surf(xx(3:end,3:end,mm),yy(3:end,3:end,mm),zz(3:end,3:end,mm),zz(3:end,3:end,mm));hold on
% 
% mm = 3;
% surf(reshape(xx(3:end,mm,3:end),Ny-1,Nz-1),reshape(yy(3:end,mm,3:end),Ny-1,Nz-1),reshape(zz(3:end,mm,3:end),Ny-1,Nz-1));hold on
% 
% mm = 3;
% surf(reshape(xx(mm,3:end,3:end),Nx-1,Nz-1),reshape(yy(mm,3:end,3:end),Nx-1,Nz-1),reshape(zz(mm,3:end,3:end),Nx-1,Nz-1));hold on
% 
% %%
% %small piece
% mm = Nz+1;
% surf(xx(3:end,1:3,mm),yy(3:end,1:3,mm),zz(3:end,1:3,mm),zz(3:end,1:3,mm));hold on
% mm = Nz+1;
% surf(xx(1:3,3:end,mm),yy(1:3,3:end,mm),zz(1:3,3:end,mm),zz(1:3,3:end,mm));hold on
% 
% mm = Nx+1;
% surf(reshape(xx(3:end,mm,1:3),Ny-1,3),reshape(yy(3:end,mm,1:3),Ny-1,3),reshape(zz(3:end,mm,1:3),Ny-1,3));hold on
% mm = Nx+1;
% surf(reshape(xx(1:3,mm,3:end),3,Nz-1),reshape(yy(1:3,mm,3:end),3,Nz-1),reshape(zz(1:3,mm,3:end),3,Nz-1));hold on
% 
% mm = Ny+1;
% surf(reshape(xx(mm,3:end,1:3),Nx-1,3),reshape(yy(mm,3:end,1:3),Nx-1,3),reshape(zz(mm,3:end,1:3),Nx-1,3));hold on
% mm = Ny+1;
% surf(reshape(xx(mm,1:3,3:end),3,Nz-1),reshape(yy(mm,1:3,3:end),3,Nz-1),reshape(zz(mm,1:3,3:end),3,Nz-1));hold on
% 
% %%
% %little piece
% mm = Nz+1;
% surf(xx(1:3,1:3,mm),yy(1:3,1:3,mm),zz(1:3,1:3,mm),zz(1:3,1:3,mm));hold on
% 
% mm = Nx+1;
% surf(reshape(xx(1:3,mm,1:3),3,3),reshape(yy(1:3,mm,1:3),3,3),reshape(zz(1:3,mm,1:3),3,3));hold on
% 
% mm = Ny+1;
% surf(reshape(xx(mm,1:3,1:3),3,3),reshape(yy(mm,1:3,1:3),3,3),reshape(zz(mm,1:3,1:3),3,3));hold on
% 
% 
% 
% 
% 
% axis equal
% view([1,1,1])