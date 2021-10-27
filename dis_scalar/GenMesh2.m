function [xx,yy,zz] = GenMesh2(Nx,Ny,Nz,delta)

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
            zz(xid,yid,zid) = (zid-1)/Nz+delta*(rand-0.5)/Nz;
            end
        end
    end
end

xx(Nx/2+1,:,:) = 0.5;

