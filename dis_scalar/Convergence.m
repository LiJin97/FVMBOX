clear 
close all
%  clc


delta = 0.2;
N = 4;

for level = 1:3
tic
%% to generate mesh 
[xxh,yyh,zzh] = GenMesh1(N,N,N,delta);

%% to obtain the location of cell centers
Pxh = (xxh(1:end-1,1:end-1,1:end-1)+xxh(1:end-1,2:end,1:end-1)+xxh(2:end,1:end-1,1:end-1)+xxh(2:end,2:end,1:end-1)+...
    xxh(1:end-1,1:end-1,2:end)+xxh(1:end-1,2:end,2:end)+xxh(2:end,1:end-1,2:end)+xxh(2:end,2:end,2:end))/8;
Pyh = (yyh(1:end-1,1:end-1,1:end-1)+yyh(1:end-1,2:end,1:end-1)+yyh(2:end,1:end-1,1:end-1)+yyh(2:end,2:end,1:end-1)+...
    yyh(1:end-1,1:end-1,2:end)+yyh(1:end-1,2:end,2:end)+yyh(2:end,1:end-1,2:end)+yyh(2:end,2:end,2:end))/8;
Pzh = (zzh(1:end-1,1:end-1,1:end-1)+zzh(1:end-1,2:end,1:end-1)+zzh(2:end,1:end-1,1:end-1)+zzh(2:end,2:end,1:end-1)+...
    zzh(1:end-1,1:end-1,2:end)+zzh(1:end-1,2:end,2:end)+zzh(2:end,1:end-1,2:end)+zzh(2:end,2:end,2:end))/8;

%% obtain the GT transformer 
Ch = GenC_cont(Pxh,Pyh,Pzh,xxh,yyh,zzh);

%% to obatin final system 
%  Ah: coefficient matrix; dh: right tern
%  uh: exact solution; mKt: the volume of a cell
[Ah,dh,uh,wh,mKt] = GenAd_new(xxh,yyh,zzh,Pxh,Pyh,Pzh,Ch);
Uh = bicgstab(Ah,dh,10^(-13),1000);

error(2,level) = CalFluxError(xxh,yyh,zzh,Uh,Ch);%This is a function to obtain error of flux
error(1,level) = sqrt(sum(mKt.*((Uh-uh).^2))); %error of solution


N = N*2;
toc
end


ratio=log2(error(:,1:end-1)./error(:,2:end));
results=[error(1,:)',[0;ratio(1,:)'],error(2,:)',[0;ratio(2,:)']];

fprintf('mesh  ||err_u||L2     Order    ||err_flux||L2    Order\n');
for i=1:level
    fprintf('%3d    %10.2e   %8.2f    %10.2e  %8.2f\n',(N/(2^(level+1)))*2^i, results(i,:));
end

