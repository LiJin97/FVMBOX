function [A,d,u,w,mKt] = GenAd_new_Robin(xx,yy,zz,Pxh,Pyh,Pzh,Ch)
N = size(Pxh,1);
d = zeros(N^3,1);
u = zeros(N^3,1);
mKt = u;
nnz = 8*1*8+12*(N-2)*12+6*((N-2)^2)*18+1*((N-2)^3)*27;
%8

val_a = zeros(N*N*N+1,1);
val_a(2) = 8;
for col = 2:N-1
    val_a(col+1) = val_a(col) + 12;
end
val_a(N+1) = val_a(N)+8;
for row = 2:N-1
    val_a((row-1)*N+2) = val_a((row-1)*N+1)+12;
    for col = 2:N-1
        val_a((row-1)*N+col+1) = val_a((row-1)*N+col)+18;
    end
    val_a((row)*N+1) = val_a((row)*N)+12;
end
val_a((N-1)*N+2) = val_a((N-1)*N+1)+8;
for col = 2:N-1
    val_a((N-1)*N+col+1) = val_a((N-1)*N+col) + 12;
end
val_a((N-1)*N+N+1) = val_a((N-1)*N+N)+8;

for pag = 2:N-1
    val_a((pag-1)*N^2+2) = val_a((pag-1)*N^2+1)+12;
    for col = 2:N-1
        val_a((pag-1)*N^2+col+1) = val_a((pag-1)*N^2+col) + 18;
    end
    val_a((pag-1)*N^2+N+1) = val_a((pag-1)*N^2+N)+12;
    
    for row = 2:N-1
        val_a((pag-1)*N^2+(row-1)*N+2) = val_a((pag-1)*N^2+(row-1)*N+1)+18;
        for col = 2:N-1
            val_a((pag-1)*N^2+(row-1)*N+col+1) = val_a((pag-1)*N^2+(row-1)*N+col)+27;
        end
        val_a((pag-1)*N^2+(row)*N+1) = val_a((pag-1)*N^2+(row)*N)+18;
    end
    val_a((pag-1)*N^2+(N-1)*N+2) = val_a((pag-1)*N^2+(N-1)*N+1)+12;
    for col = 2:N-1
        val_a((pag-1)*N^2+(N-1)*N+col+1) = val_a((pag-1)*N^2+(N-1)*N+col) + 18;
    end
    val_a((pag-1)*N^2+(N-1)*N+N+1) = val_a((pag-1)*N^2+(N-1)*N+N)+12;
    
end

pag = N;

val_a((pag-1)*N^2+2) = val_a((pag-1)*N^2+1)+8;
for col = 2:N-1
    val_a((pag-1)*N^2+col+1) = val_a((pag-1)*N^2+col) + 12;
end
val_a((pag-1)*N^2+N+1) = val_a((pag-1)*N^2+N)+8;

for row = 2:N-1
    val_a((pag-1)*N^2+(row-1)*N+2) = val_a((pag-1)*N^2+(row-1)*N+1)+12;
    for col = 2:N-1
        val_a((pag-1)*N^2+(row-1)*N+col+1) = val_a((pag-1)*N^2+(row-1)*N+col)+18;
    end
    val_a((pag-1)*N^2+(row)*N+1) = val_a((pag-1)*N^2+(row)*N)+12;
end
val_a((pag-1)*N^2+(N-1)*N+2) = val_a((pag-1)*N^2+(N-1)*N+1)+8;
for col = 2:N-1
    val_a((pag-1)*N^2+(N-1)*N+col+1) = val_a((pag-1)*N^2+(N-1)*N+col) + 12;
end
val_a((pag-1)*N^2+(N-1)*N+N+1) = val_a((pag-1)*N^2+(N-1)*N+N)+8;

CsrInd = zeros(nnz,2);


% ind27 = zeros(3,3,3);
ind27 = [[-N-1, -N, -N+1, -1, 0, 1, N-1, N, N+1]'-N*N;[-N-1, -N, -N+1, -1, 0, 1, N-1, N, N+1]';[-N-1, -N, -N+1, -1, 0, 1, N-1, N, N+1]'+N*N];
% Coeff27 = zeros(3,3,3);
volum = 0;
count = 0;
for zid = 1:N
    for yid = 1:N
        for xid = 1:N
            BoxP = [reshape(xx(xid:xid+1,yid:yid+1,zid:zid+1),[],1),reshape(yy(xid:xid+1,yid:yid+1,zid:zid+1),[],1),reshape(zz(xid:xid+1,yid:yid+1,zid:zid+1),[],1)];
            %BoxP 每一行对应着节点的x,y,z坐标
            Pcc = [Pxh(xid,yid,zid),Pyh(xid,yid,zid),Pzh(xid,yid,zid)]; 
            
            Coeff27 = zeros(3,3,3);
            Coeff8 = zeros(8,1);
            CenId = ones(3,3,3);
            VerId = ones(8,1);
            temp = 0;
            mK = 0;
            %27 18 12 8
            if xid==1
                CenId(1,:,:) = 0;
                
                indt = [1,5,7,3];
                VerId(indt) = 0;

                 Prc = [Pxh(xid+1,yid,zid),Pyh(xid+1,yid,zid),Pzh(xid+1,yid,zid)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffb(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                temp = temp+tmp;
                Coeff27([1,2],2,2) = Coeff27([1,2],2,2)+[-Cpk,Cpk]';
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                mK = mK+tmK;
                
                indt = [2,4,8,6]; 
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffc_new1_Robin(Pcc',Prc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
%                 [Cpk,Cp1,Cp2,Cp3,Cp4]
                mK = mK+tmK;
                Coeff27([3,2],2,2) = Coeff27([3,2],2,2)+[-Cpk,Cpk]';
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                temp = temp+tmp;
            elseif xid==N
                CenId(3,:,:) = 0;  
                indt = [1,5,7,3];  
                
                Plc = [Pxh(xid-1,yid,zid),Pyh(xid-1,yid,zid),Pzh(xid-1,yid,zid)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffc_new1_Robin(Pcc',Plc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                temp = temp+tmp;
                mK = mK+tmK;
                Coeff27([1,2],2,2) = Coeff27([1,2],2,2)+[-Cpk,Cpk]';
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                
                indt = [2,4,8,6];
                VerId(indt) = 0;
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffb_Robin(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                Coeff27([3,2],2,2) = Coeff27([3,2],2,2)+[-Cpk,Cpk]';
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                temp = temp+tmp;
                mK = mK+tmK;
            else
                
                indt = [1,5,7,3];                
                Plc = [Pxh(xid-1,yid,zid),Pyh(xid-1,yid,zid),Pzh(xid-1,yid,zid)];
                Prc = [Pxh(xid+1,yid,zid),Pyh(xid+1,yid,zid),Pzh(xid+1,yid,zid)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffc_new1_Robin(Pcc',Plc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                temp = temp+tmp;
                mK = mK+tmK;
                Coeff27([1,2],2,2) = Coeff27([1,2],2,2)+[-Cpk,Cpk]';
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                
                indt = [2,4,8,6];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffc_new1_Robin(Pcc',Prc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                Coeff27([3,2],2,2) = Coeff27([3,2],2,2)+[-Cpk,Cpk]';
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                temp = temp+tmp;
                mK = mK+tmK;
            end
            
            if yid==1
                CenId(:,1,:) = 0;
                indt = [1,2,6,5];                
                VerId(indt) = 0;

                Pbc = [Pxh(xid,yid+1,zid),Pyh(xid,yid+1,zid),Pzh(xid,yid+1,zid)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffb_Robin(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                temp = temp+tmp;
                Coeff27(2,[1,2],2) = Coeff27(2,[1,2],2)+[-Cpk,Cpk];
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                mK = mK+tmK;
                indt = [3,7,8,4]; 
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffc_new1_Robin(Pcc',Pbc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                Coeff27(2,[3,2],2) = Coeff27(2,[3,2],2)+[-Cpk,Cpk];
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                temp = temp+tmp;
                mK = mK+tmK;
            elseif yid==N
                CenId(:,3,:) = 0;
                indt = [1,2,6,5];                
                
                Pfc = [Pxh(xid,yid-1,zid),Pyh(xid,yid-1,zid),Pzh(xid,yid-1,zid)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffc_new1_Robin(Pcc',Pfc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                temp = temp+tmp;
                mK = mK+tmK;
                Coeff27(2,[1,2],2) = Coeff27(2,[1,2],2)+[-Cpk,Cpk];
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                
                indt = [3,7,8,4];
                VerId(indt) = 0;
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffb_Robin(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                Coeff27(2,[3,2],2) = Coeff27(2,[3,2],2)+[-Cpk,Cpk];
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                temp = temp+tmp;
                mK = mK+tmK;
            else
                indt = [1,2,6,5];                
                
                Pfc = [Pxh(xid,yid-1,zid),Pyh(xid,yid-1,zid),Pzh(xid,yid-1,zid)];
                Pbc = [Pxh(xid,yid+1,zid),Pyh(xid,yid+1,zid),Pzh(xid,yid+1,zid)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffc_new1_Robin(Pcc',Pfc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                Coeff27(2,[1,2],2) = Coeff27(2,[1,2],2)+[-Cpk,Cpk];
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                temp = temp+tmp;
                mK = mK+tmK;
                indt = [3,7,8,4]; 
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffc_new1_Robin(Pcc',Pbc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                Coeff27(2,[3,2],2) = Coeff27(2,[3,2],2)+[-Cpk,Cpk];
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                temp = temp+tmp;
                mK = mK+tmK;
            end
            
            if zid==1
                CenId(:,:,1) = 0;
                indt = [1,3,4,2];                
                VerId(indt) = 0;

                 Puc = [Pxh(xid,yid,zid+1),Pyh(xid,yid,zid+1),Pzh(xid,yid,zid+1)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffb_Robin(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                temp = temp+tmp;
                Coeff27(2,2,[1,2]) = Coeff27(2,2,[1,2])+reshape([-Cpk,Cpk],1,1,2);
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                mK = mK+tmK;
                indt = [5,6,8,7]; 
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffc_new1_Robin(Pcc',Puc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                Coeff27(2,2,[3,2]) = Coeff27(2,2,[3,2])+reshape([-Cpk,Cpk],1,1,2);
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                temp = temp+tmp;
                mK = mK+tmK;
            elseif zid==N
                CenId(:,:,3) = 0;
                indt = [1,3,4,2];                

                Pdc = [Pxh(xid,yid,zid-1),Pyh(xid,yid,zid-1),Pzh(xid,yid,zid-1)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffc_new1_Robin(Pcc',Pdc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                temp = temp+tmp;
                mK = mK+tmK;
                Coeff27(2,2,[1,2]) = Coeff27(2,2,[1,2])+reshape([-Cpk,Cpk],1,1,2);
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                
                indt = [5,6,8,7];   
                VerId(indt) = 0;
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffb_Robin(Pcc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                Coeff27(2,2,[3,2]) = Coeff27(2,2,[3,2])+reshape([-Cpk,Cpk],1,1,2);
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                temp = temp+tmp;
                mK = mK+tmK;
                
            else
                indt = [1,3,4,2];                
                Pdc = [Pxh(xid,yid,zid-1),Pyh(xid,yid,zid-1),Pzh(xid,yid,zid-1)];
                Puc = [Pxh(xid,yid,zid+1),Pyh(xid,yid,zid+1),Pzh(xid,yid,zid+1)];
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffc_new1_Robin(Pcc',Pdc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                Coeff27(2,2,[1,2]) = Coeff27(2,2,[1,2])+reshape([-Cpk,Cpk],1,1,2);
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                temp = temp+tmp;
                mK = mK+tmK;
                indt = [5,6,8,7]; 
                [Cpk,Cp1,Cp2,Cp3,Cp4,tmp,tmK] = CalCoeffc_new1_Robin(Pcc',Puc',BoxP(indt(1),:)',BoxP(indt(2),:)',BoxP(indt(3),:)',BoxP(indt(4),:)');
                Coeff27(2,2,[3,2]) = Coeff27(2,2,[3,2])+reshape([-Cpk,Cpk],1,1,2);
                Coeff8(indt) = Coeff8(indt) + [Cp1,Cp2,Cp3,Cp4]';
                temp = temp+tmp;
                mK = mK+tmK;
            end
            VerId = reshape(VerId,2,2,2);
            for pz = 1:2
                for py = 1:2
                    for px = 1:2
                        pp = (pz-1)*4+(py-1)*2+px;                          
                            
                            ww = [Ch(2,xid+px-1,yid+py-1,zid+pz-1);Ch(1,xid+px-1,yid+py-1,zid+pz-1);Ch(3,xid+px-1,yid+py-1,zid+pz-1);Ch(4,xid+px-1,yid+py-1,zid+pz-1);...
                                Ch(6,xid+px-1,yid+py-1,zid+pz-1);Ch(5,xid+px-1,yid+py-1,zid+pz-1);Ch(7,xid+px-1,yid+py-1,zid+pz-1);Ch(8,xid+px-1,yid+py-1,zid+pz-1)];
                           
                            
                            Coeff27(px:px+1,py:py+1,pz:pz+1) = Coeff27(px:px+1,py:py+1,pz:pz+1)+Coeff8(pp)*reshape(ww,2,2,2);
                            temp = temp + Coeff8(pp)*Ch(9,xid+px-1,yid+py-1,zid+pz-1);
                         
%                         
                    end
                end
            end
            
            
            
            
            Cid = sum(CenId(:));
            ind = (zid-1)*N*N + (yid-1)*N + xid;
            
            
            CsrInd(count+1:count+Cid,1) = ind + ind27(CenId~=0)';
            CsrInd(count+1:count+Cid,2) = Coeff27(CenId~=0)';

            d(ind) = mK*fK_Robin(Pcc')-temp;%QuadA(BoxP(1,:),BoxP(2,:),BoxP(3,:),BoxP(4,:),BoxP(5,:),BoxP(6,:),BoxP(7,:),BoxP(8,:),Pcc)*fK(Pcc')-temp;
            u(ind) = GenReal_Robin(Pcc');
            count = count+Cid;
            mKt(xid+(yid-1)*N+(zid-1)*N*N) = mK;
        end
    end
end
% volum
% CsrInd
w = zeros((N+1)^3,1);
% for i = 1:N+1
%     for j = 1:N+1
%         for k = 1:N+1
% %             d((i-1)*(N)+j+(k-1)*N*N) =  (1/N/N)*fK((j-1+0.5)/N,(i-1+0.5)/N,(k-1+0.5)/N);
% %             u((i-1)*(N)+j+(k-1)*N*N) = GenReal((j-1+0.5)/N,(i-1+0.5)/N,(k-1+0.5)/N);
%             w((i-1)*(N+1)+j+(k-1)*(N+1)*(N+1)) = GenReal([xx(j,i,k),yy(j,i,k),zz(j,i,k)]);
%         end
%     end
% end

% reshape(d(1:end-1),3,[])'
% reshape(CsrInd(1:102,2),3,[])'
A = csr2sparse(CsrInd(:,2),val_a,CsrInd(:,1));