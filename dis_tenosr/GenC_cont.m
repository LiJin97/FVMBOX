function Ch = GenC_cont(Px,Py,Pz,XX,YY,ZZ)
Ny = size(Px,1);
Nx = size(Px,2);
Nz = size(Px,3);
Ch = zeros(8,(Ny-1)*(Nx-1)*(Nz-1));


for i = 1:Ny-1
    for j = 1:Nx-1
        for k = 1:Nz-1
            W = zeros(4,8);
            W(1,:) = 1;
            A = [XX(i+1,j+1,k+1),YY(i+1,j+1,k+1),ZZ(i+1,j+1,k+1)];
            
            tu = [XX(i+1,j+1,k+2),YY(i+1,j+1,k+2),ZZ(i+1,j+1,k+2)]-A;
            td = [XX(i+1,j+1,k),YY(i+1,j+1,k),ZZ(i+1,j+1,k)]-A;
            tf = [XX(i+1,j,k+1),YY(i+1,j,k+1),ZZ(i+1,j,k+1)]-A;
            tb = [XX(i+1,j+2,k+1),YY(i+1,j+2,k+1),ZZ(i+1,j+2,k+1)]-A;
            tl = [XX(i,j+1,k+1),YY(i,j+1,k+1),ZZ(i,j+1,k+1)]-A;
            tr = [XX(i+2,j+1,k+1),YY(i+2,j+1,k+1),ZZ(i+2,j+1,k+1)]-A;
            
            
            da = [Px(i,j,k),Py(i,j,k),Pz(i,j,k)];
            db = [Px(i+1,j,k),Py(i+1,j,k),Pz(i+1,j,k)];
            dc = [Px(i,j+1,k),Py(i,j+1,k),Pz(i,j+1,k)];
            dd = [Px(i+1,j+1,k),Py(i+1,j+1,k),Pz(i+1,j+1,k)];
            ua = [Px(i,j,k+1),Py(i,j,k+1),Pz(i,j,k+1)];
            ub = [Px(i+1,j,k+1),Py(i+1,j,k+1),Pz(i+1,j,k+1)];
            uc = [Px(i,j+1,k+1),Py(i,j+1,k+1),Pz(i,j+1,k+1)];
            ud = [Px(i+1,j+1,k+1),Py(i+1,j+1,k+1),Pz(i+1,j+1,k+1)];
            
            kda = Genk(da);
            kdb = Genk(db);
            kdc = Genk(dc);
            kdd = Genk(dd);
            
            kua = Genk(ua);
            kub = Genk(ub);
            kuc = Genk(uc);
            kud = Genk(ud);
            
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
            
%             W(2:4,1) = kda(1)\(da-A)';
%             W(2:4,2) = kdb(1)\(db-A)';
%             W(2:4,3) = kdc(1)\(dc-A)';
%             W(2:4,4) = kdd(1)\(dd-A)';
%             W(2:4,5) = kua(1)\(ua-A)';
%             W(2:4,6) = kub(1)\(ub-A)';
%             W(2:4,7) = kuc(1)\(uc-A)';
%             W(2:4,8) = kud(1)\(ud-A)';
            
            ww = W'*((W*W')\[1;0;0;0]);
            
            
            Ch(:,(k-1)*(Nx-1)*(Ny-1)+(j-1)*(Nx-1)+i) = ww;
        end
    end
end

% reshape(W(1:end-2),3,[])'