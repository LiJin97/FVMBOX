%% Main function 
clear 
close all
clc


delta = 0.8; % This is a parameter to control distortion of mesh.
N = 4;
tic
for level = 1:2

[xxh,yyh,zzh] = GenMesh1(N,N,N,delta); %% There are two kinds of 3D meshes. GenMesh1(GenMesh2) is distorted(random) mesh.

%% To obtain the location of cell centers
Pxh = (xxh(1:end-1,1:end-1,1:end-1)+xxh(1:end-1,2:end,1:end-1)+xxh(2:end,1:end-1,1:end-1)+xxh(2:end,2:end,1:end-1)+...
    xxh(1:end-1,1:end-1,2:end)+xxh(1:end-1,2:end,2:end)+xxh(2:end,1:end-1,2:end)+xxh(2:end,2:end,2:end))/8;
Pyh = (yyh(1:end-1,1:end-1,1:end-1)+yyh(1:end-1,2:end,1:end-1)+yyh(2:end,1:end-1,1:end-1)+yyh(2:end,2:end,1:end-1)+...
    yyh(1:end-1,1:end-1,2:end)+yyh(1:end-1,2:end,2:end)+yyh(2:end,1:end-1,2:end)+yyh(2:end,2:end,2:end))/8;
Pzh = (zzh(1:end-1,1:end-1,1:end-1)+zzh(1:end-1,2:end,1:end-1)+zzh(2:end,1:end-1,1:end-1)+zzh(2:end,2:end,1:end-1)+...
    zzh(1:end-1,1:end-1,2:end)+zzh(1:end-1,2:end,2:end)+zzh(2:end,1:end-1,2:end)+zzh(2:end,2:end,2:end))/8;

%% To obtain matrix of gradient transfer
Ch = GenC_cont(Pxh,Pyh,Pzh,xxh,yyh,zzh);

%% Ah is coefficient matrix, dh is right term, uh is excat solution, mKt is volume of the cell.
[Ah,dh,uh,~,mKt] = GenAd_new(xxh,yyh,zzh,Pxh,Pyh,Pzh,Ch);

Uh = Ah\dh;
%% These are errors of slution and flux.
error(2,level) =CalFluxError(xxh,yyh,zzh,Uh,Ch);
error(1,level) = sqrt(sum(mKt.*((Uh-uh).^2)));

N = N*2;
toc
end


ratio=log2(error(:,1:end-1)./error(:,2:end));
results=[error(1,:)',[0;ratio(1,:)'],error(2,:)',[0;ratio(2,:)']];

fprintf('mesh  ||err_u||L2       Order    ||err_flux||L2    Order\n');
for i=1:level
    fprintf('%3d    %10.2e   %8.2f    %10.2e     %8.2f\n',(N/(2^(level+1)))*2^i, results(i,:));
end

