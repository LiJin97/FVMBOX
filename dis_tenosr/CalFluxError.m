function err = CalFluxError(xx,yy,zz,Uh,Ch)
%%
%This is a function to obtain the error of fluxes.
%%

Pxh = (xx(1:end-1,1:end-1,1:end-1)+xx(1:end-1,2:end,1:end-1)+xx(2:end,1:end-1,1:end-1)+xx(2:end,2:end,1:end-1)+...
    xx(1:end-1,1:end-1,2:end)+xx(1:end-1,2:end,2:end)+xx(2:end,1:end-1,2:end)+xx(2:end,2:end,2:end))/8;
Pyh = (yy(1:end-1,1:end-1,1:end-1)+yy(1:end-1,2:end,1:end-1)+yy(2:end,1:end-1,1:end-1)+yy(2:end,2:end,1:end-1)+...
    yy(1:end-1,1:end-1,2:end)+yy(1:end-1,2:end,2:end)+yy(2:end,1:end-1,2:end)+yy(2:end,2:end,2:end))/8;
Pzh = (zz(1:end-1,1:end-1,1:end-1)+zz(1:end-1,2:end,1:end-1)+zz(2:end,1:end-1,1:end-1)+zz(2:end,2:end,1:end-1)+...
    zz(1:end-1,1:end-1,2:end)+zz(1:end-1,2:end,2:end)+zz(2:end,1:end-1,2:end)+zz(2:end,2:end,2:end))/8;
N = size(Pxh,1);
Vh = zeros(N+1,N+1,N+1);
for xid = 1:N+1
    for yid = 1:N+1
        for zid = 1:N+1
           if xid == 1 || xid == N+1 || yid == 1 || yid == N+1|| zid == 1 || zid == N+1 
            Vh(xid,yid,zid) = GenReal([xx(xid,yid,zid),yy(xid,yid,zid),zz(xid,yid,zid)]);
           else     
            Uda = Uh((zid-1-1)*N*N+(yid-1-1)*N+xid-1);
            Udb = Uh((zid-1-1)*N*N+(yid-1-1)*N+xid);
            Udc = Uh((zid-1-1)*N*N+(yid-1)*N+xid-1);
            Udd = Uh((zid-1-1)*N*N+(yid-1)*N+xid);
            Uua = Uh((zid-1)*N*N+(yid-1-1)*N+xid-1);
            Uub = Uh((zid-1)*N*N+(yid-1-1)*N+xid);
            Uuc = Uh((zid-1)*N*N+(yid-1)*N+xid-1); 
            Uud = Uh((zid-1)*N*N+(yid-1)*N+xid);
            u = [Uda,Udb,Udc,Udd,Uua,Uub,Uuc,Uud]';
            Vh(xid,yid,zid) = Ch(:,(zid-2)*(N-1)*(N-1)+(yid-2)*(N-1)+xid-1)'*u;
            end
        end
    end
end



err = 0;
for zid = 1:N
    for yid = 1:N
        for xid = 1:N
            BoxP = [reshape(xx(xid:xid+1, yid:yid+1, zid:zid+1),[],1),reshape(yy(xid:xid+1,yid:yid+1,zid:zid+1),[],1),reshape(zz(xid:xid+1,yid:yid+1,zid:zid+1),[],1)];
            VoxP = reshape(Vh(xid:xid+1, yid:yid+1, zid:zid+1),[],1);
            
            Pcc = [Pxh(xid,yid,zid),Pyh(xid,yid,zid),Pzh(xid,yid,zid)];
            if xid==1
          
                
                indt = [1,5,7,3];
                
                Prc = [Pxh(xid+1,yid,zid),Pyh(xid+1,yid,zid),Pzh(xid+1,yid,zid)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffb(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');    
                err = err + (Cpk*Uh(xid+(yid-1)*N+(zid-1)*N*N) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)+tmp - GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;
                indt = [2,4,8,6];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmK] = CalCoeffc_new1(Pcc',Prc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                err = err + (Cpk*(Uh(xid+(yid-1)*N+(zid-1)*N*N)-Uh(xid+1+(yid-1)*N+(zid-1)*N*N)) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)-GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;
                
                
            elseif xid==N
                indt = [1,5,7,3];
                
                Plc = [Pxh(xid-1,yid,zid),Pyh(xid-1,yid,zid),Pzh(xid-1,yid,zid)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmK] = CalCoeffc_new1(Pcc',Plc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                err = err + (Cpk*(Uh(xid+(yid-1)*N+(zid-1)*N*N)-Uh(xid-1+(yid-1)*N+(zid-1)*N*N)) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)-GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;

                
                indt = [2,4,8,6];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffb(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                err = err + (Cpk*Uh(xid+(yid-1)*N+(zid-1)*N*N) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)+tmp - GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;

            else
                
                indt = [1,5,7,3];
                Plc = [Pxh(xid-1,yid,zid),Pyh(xid-1,yid,zid),Pzh(xid-1,yid,zid)];
                Prc = [Pxh(xid+1,yid,zid),Pyh(xid+1,yid,zid),Pzh(xid+1,yid,zid)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmK] = CalCoeffc_new1(Pcc',Plc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                err = err + (Cpk*(Uh(xid+(yid-1)*N+(zid-1)*N*N)-Uh(xid-1+(yid-1)*N+(zid-1)*N*N)) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)-GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;

                indt = [2,4,8,6];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmK] = CalCoeffc_new1(Pcc',Prc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                err = err + (Cpk*(Uh(xid+(yid-1)*N+(zid-1)*N*N)-Uh(xid+1+(yid-1)*N+(zid-1)*N*N)) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)-GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;
                
            end
            
            if yid==1
                indt = [1,2,6,5];
                
                Pbc = [Pxh(xid,yid+1,zid),Pyh(xid,yid+1,zid),Pzh(xid,yid+1,zid)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffb(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                err = err + (Cpk*Uh(xid+(yid-1)*N+(zid-1)*N*N) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)+tmp - GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;

                indt = [3,7,8,4];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmK] = CalCoeffc_new1(Pcc',Pbc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                err = err + (Cpk*(Uh(xid+(yid-1)*N+(zid-1)*N*N)-Uh(xid+(yid-1+1)*N+(zid-1)*N*N)) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)-GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;
                
            elseif yid==N
                indt = [1,2,6,5];
                
                Pfc = [Pxh(xid,yid-1,zid),Pyh(xid,yid-1,zid),Pzh(xid,yid-1,zid)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmK] = CalCoeffc_new1(Pcc',Pfc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                err = err + (Cpk*(Uh(xid+(yid-1)*N+(zid-1)*N*N)-Uh(xid+(yid-1-1)*N+(zid-1)*N*N)) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)-GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;
                
                
                indt = [3,7,8,4];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffb(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                err = err + (Cpk*Uh(xid+(yid-1)*N+(zid-1)*N*N) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)+tmp - GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;

            else
                indt = [1,2,6,5];
                
                Pfc = [Pxh(xid,yid-1,zid),Pyh(xid,yid-1,zid),Pzh(xid,yid-1,zid)];
                Pbc = [Pxh(xid,yid+1,zid),Pyh(xid,yid+1,zid),Pzh(xid,yid+1,zid)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmK] = CalCoeffc_new1(Pcc',Pfc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                 err = err + (Cpk*(Uh(xid+(yid-1)*N+(zid-1)*N*N)-Uh(xid+(yid-1-1)*N+(zid-1)*N*N)) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)-GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;
                
                indt = [3,7,8,4];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmK] = CalCoeffc_new1(Pcc',Pbc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                 err = err + (Cpk*(Uh(xid+(yid-1)*N+(zid-1)*N*N)-Uh(xid+(yid-1+1)*N+(zid-1)*N*N)) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)-GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;
                
            end
            
            if zid==1
                indt = [1,3,4,2];
                
                Puc = [Pxh(xid,yid,zid+1),Pyh(xid,yid,zid+1),Pzh(xid,yid,zid+1)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffb(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                err = err + (Cpk*Uh(xid+(yid-1)*N+(zid-1)*N*N) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)+tmp - GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;

                indt = [5,6,8,7];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmK] = CalCoeffc_new1(Pcc',Puc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                err = err + (Cpk*(Uh(xid+(yid-1)*N+(zid-1)*N*N)-Uh(xid+(yid-1)*N+(zid-1+1)*N*N)) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)-GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;
                
            elseif zid==N
                indt = [1,3,4,2];
                
                Pdc = [Pxh(xid,yid,zid-1),Pyh(xid,yid,zid-1),Pzh(xid,yid,zid-1)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmK] = CalCoeffc_new1(Pcc',Pdc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                err = err + (Cpk*(Uh(xid+(yid-1)*N+(zid-1)*N*N)-Uh(xid+(yid-1)*N+(zid-1-1)*N*N)) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)-GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;
                
                
                indt = [5,6,8,7];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffb(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                err = err + (Cpk*Uh(xid+(yid-1)*N+(zid-1)*N*N) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)+tmp - GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;

            else
                indt = [1,3,4,2];
                Pdc = [Pxh(xid,yid,zid-1),Pyh(xid,yid,zid-1),Pzh(xid,yid,zid-1)];
                Puc = [Pxh(xid,yid,zid+1),Pyh(xid,yid,zid+1),Pzh(xid,yid,zid+1)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmK] = CalCoeffc_new1(Pcc',Pdc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                err = err + (Cpk*(Uh(xid+(yid-1)*N+(zid-1)*N*N)-Uh(xid+(yid-1)*N+(zid-1-1)*N*N)) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)-GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;
                
                indt = [5,6,8,7];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmK] = CalCoeffc_new1(Pcc',Puc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                err = err + (Cpk*(Uh(xid+(yid-1)*N+(zid-1)*N*N)-Uh(xid+(yid-1)*N+(zid-1+1)*N*N)) + [Cp1,Cp2,Cp3,Cp4]*VoxP(indt)-GenFlux(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)'))^2;
               
            end
            
            
            
        end
    end
end
err = sqrt(err);
