function [recon,info]=diffmap_phase_retrieval(u,opts)
if nargin<1
    diffmap_phase_retrieval_example1;
    return;
end;

if (~isfield(opts,'eps_hack')) opts.eps_hack=0; end;
if (~isfield(opts,'alpha_A')) opts.alpha_A=1; end;
if (~isfield(opts,'alpha_B')) opts.alpha_B=1; end;

%initial_phase=zeros(size(u));
initial_phase=rand(size(u))*2*pi;
f0=real(ifft2b(u.*exp(i*initial_phase)));
opts.u=u;
info=struct;

cmap='gray';

%figure; imagesc(rand(10,10)*2-1); set(gca,'clim',[-1,1]); colorbar; colormap(cmap);
fA=figure; set(fA,'position',[100,100,2600,1300]);
fB=figure; set(fB,'position',[100,1500,300,300]);
resids=[];
tA=tic;
if (isfield(opts,'reference'))
    figure; imagesc(opts.reference); set(gca,'clim',[-1,1]); colorbar; colormap(cmap);
end;

fA=figure; set(fA,'position',[100,100,1300,1300]);
for it=1:opts.max_iterations
    %f1 = pi_A(2pi_B(f0)-f0) - pi_B(f0)
    if strcmp(opts.diffmap_method,'AB')
        pi_B_f0=pi_B(f0,opts);
        proj1=pi_A(2*pi_B_f0-f0,opts);
        proj2=pi_B_f0;
        f1=f0+proj1-proj2;
        eps0=opts.eps_hack;
        f1=(abs(f1)<eps0).*(eps0*sign(f1))+(abs(f1)>=eps0).*f1;
        f0=f1;
        figure(fA);
        subplot(1,1,1); imagesc(pi_B(f0,opts)); set(gca,'clim',[-1,1]);  colormap(cmap);
        %subplot(1,2,2); imagesc(f0); set(gca,'clim',[-1,1]);  colormap(cmap);
    elseif strcmp(opts.diffmap_method,'AB_HIO')
        pi_B_f0=pi_B(f0,opts);
        %beta=min(1,0.5+0.5*it/(opts.max_iterations*0.9));
        beta=0.5;
        proj1=pi_A((1+beta)*pi_B_f0-f0,opts);
        proj2=pi_B_f0;
        f1=f0+proj1-beta*proj2;
        eps0=opts.eps_hack;
        f1=(abs(f1)<eps0).*(eps0*sign(f1))+(abs(f1)>=eps0).*f1;
        f0=f1;
        tmp=f0.*opts.support_mask;
        if mean(tmp(:)<0), f0=-f0; end;
        figure(fA);
        subplot(1,1,1); imagesc(pi_B(f0,opts)); set(gca,'clim',[-1,1]);  colormap(cmap);
        %subplot(1,2,2); imagesc(f0); set(gca,'clim',[-1,1]); colormap(cmap);
    elseif strcmp(opts.diffmap_method,'BA')
        pi_A_f0=pi_A(f0,opts);
        proj1=pi_B(2*pi_A_f0-f0,opts);
        proj2=pi_A_f0;
        f1=f0+proj1-proj2;
        eps0=opts.eps_hack;
        f1=(abs(f1)<eps0).*(eps0*sign(f1))+(abs(f1)>=eps0).*f1;
        f0=f1;
        figure(fA);
        subplot(1,1,1); imagesc(pi_A(f0,opts)); set(gca,'clim',[-1,1]);  colormap(cmap);
        %subplot(1,2,2); imagesc(f0); set(gca,'clim',[-1,1]); colormap(cmap);
    elseif strcmp(opts.diffmap_method,'BA_HIO')
        pi_A_f0=pi_A(f0,opts);
        %beta=min(1,0.5+0.5*it/(opts.max_iterations*0.9));
        beta=0.5;
        proj1=pi_B((1+beta)*pi_A_f0-f0,opts);
        proj2=pi_A_f0;
        f1=f0+proj1-beta*proj2;
        eps0=opts.eps_hack;
        f1=(abs(f1)<eps0).*(eps0*sign(f1))+(abs(f1)>=eps0).*f1;
        f0=f1;
        tmp=f0.*opts.support_mask;
        if mean(tmp(:)<0), f0=-f0; end;
        figure(fA);
        subplot(1,1,1); imagesc(pi_A(f0,opts)); set(gca,'clim',[-1,1]);  colormap(cmap);
        %subplot(1,2,2); imagesc(f0); set(gca,'clim',[-1,1]); colormap(cmap);
    elseif strcmp(opts.diffmap_method,'AP')
        pi_A_f0=pi_A(f0,opts);
        f1=pi_B(pi_A_f0,opts);
        eps0=opts.eps_hack;
        f1=(abs(f1)<eps0).*(eps0*sign(f1))+(abs(f1)>=eps0).*f1;
        f0=f1;
        figure(fA);
        subplot(1,1,1); imagesc(pi_A(f0,opts)); set(gca,'clim',[-1,1]);  colormap(cmap);
    elseif strcmp(opts.diffmap_method,'HIO')
        pi_A_f0=pi_A(f0,opts);
        beta=0.5;
        f1=pi_A_f0.*opts.support_mask + (f0-beta*pi_A_f0).*(1-opts.support_mask);
        f0=f1;
        tmp=f0.*opts.support_mask;
        if mean(tmp(:)<0), f0=-f0; end;
        figure(fA);
        subplot(1,1,1); imagesc(f0); set(gca,'clim',[-1,1]);  colormap(cmap);
    end;
    
    if ~strcmp(opts.diffmap_method,'HIO')
        resid0=compute_residual(proj1,proj2);
        title(sprintf('iteration %d - resid between pi_A and pi_B = %g',it,resid0));
        resids(it)=resid0;
        figure(fB);
        semilogy(1:it,resids);
    end;
    if (toc(tA)>0.3)
        drawnow;
        tA=tic;
    end;
end;

if strcmp(opts.diffmap_method,'AB')
    recon=pi_B(f0,opts);
elseif strcmp(opts.diffmap_method,'AB_HIO')
    recon=pi_B(f0,opts);
elseif strcmp(opts.diffmap_method,'BA')
    recon=pi_A(f0,opts);
elseif strcmp(opts.diffmap_method,'BA_HIO')
    recon=pi_A(f0,opts);
elseif strcmp(opts.diffmap_method,'AP')
    recon=pi_A(f0,opts);
end;

if (isfield(opts,'reference'))
    recon=register_to_reference(recon,opts.reference);
end;

end

function f2=pi_A(f0,opts)
f1=fft2b(f0);
f1=opts.u.*exp(i*angle(f1));
f1=real(ifft2b(f1));
f2=f0+opts.alpha_A*(f1-f0);
end

function f0=pi_B(f0,opts)
if opts.positivity
    f0=f0.*(f0>=0);
end;
if opts.support
    f0=f0.*opts.support_mask;
elseif opts.support2
    if opts.positivity
        mask0=find_optimal_mask(f0.*(f0>=0),opts.support_mask);
    else
        mask0=find_optimal_mask(f0,opts.support_mask);
    end;
    f0=f0.*mask0;
end;
end

function mask0=find_optimal_mask(f,mask)

tmp=abs(ifftn(fftn(abs(mask).^2).*conj(fftn(f))));
[~,ind0]=max(tmp(:));
[i1,i2]=ind2sub(size(tmp),ind0);

mask0=circshift(circshift(mask,-i1,1),-i2,2);

end


% function f0=find_optimal_centered_image(f,mask)
% 
% tmp=abs(ifftn(fftn(abs(f).^2).*conj(fftn(mask))));
% [~,ind0]=max(tmp(:));
% [i1,i2]=ind2sub(size(tmp),ind0);
% 
% f0=circshift(circshift(f,-i1,1),-i2,2);
% 
% end

function Y=fft2b(X)
Y=fftshift(fft2(fftshift(X)));
end

function Y=ifft2b(X)
Y=fftshift(ifft2(fftshift(X)));
end

function A=register_to_reference(A,A_ref)

N1=size(A_ref,1);
N2=size(A_ref,2);
N3=size(A_ref,3);
M1=ceil((N1+1)/2);
M2=ceil((N2+1)/2);
M3=ceil((N3+1)/2);

tmp1=abs(ifft0(fft0(A).*conj(fft0(A_ref))));
tmp2=abs(ifft0(conj(fft0(A)).*conj(fft0(A_ref))));
[max1,ii1]=max(tmp1(:));
[max2,ii2]=max(tmp2(:));

if (max1>max2)
    ii=ii1;
else
    ii=ii2;
    A=ifft0(conj(fft0(A)));
end;

[i1,i2,i3]=ind2sub(size(A),ii);
A=circshift(A,-(i1-M1),1);
A=circshift(A,-(i2-M2),2);
A=circshift(A,-(i3-M3),3);

dxrange=-2:0.2:2;
dyrange=-2:0.2:2;
dzrange=-2:0.2:2;
if (N3==1) dzrange=0; end;
best_resid=inf; best_A=A;
for dz=dzrange
for dy=dyrange
for dx=dxrange
    Atry=fracshift(A,dx,dy,dz);
    resid=Atry-A_ref; resid=sum(resid(:).^2);
    if (resid<best_resid)
        best_A=Atry;
        best_resid=resid;
    end;
end;
end;
end;

A=best_A;

% [GX,GY,GZ]=ndgrid((0:N1-1)-M1,(0:N2-1)-M2,(0:N3-1)-M3);
% A1_CX=sum(A_ref(:).*GX(:))/sum(A_ref(:));
% A1_CY=sum(A_ref(:).*GY(:))/sum(A_ref(:));
% A1_CZ=sum(A_ref(:).*GZ(:))/sum(A_ref(:));
% 
% A2_CX=sum(A(:).*GX(:))/sum(A(:));
% A2_CY=sum(A(:).*GY(:))/sum(A(:));
% A2_CZ=sum(A(:).*GZ(:))/sum(A(:));
% 
% DX=A2_CX-A1_CX;
% DY=A2_CY-A1_CY;
% DZ=A2_CZ-A1_CZ;
% 
% A2=real(ifft0(fft0(A).*exp(2*pi*i*(DX*GX/N1+DY*GY/N2+DZ*GZ/N3))));

end

function resid=compute_residual(f,ref)
resid=sqrt(sum((f(:)-ref(:)).^2))/sqrt(sum(ref(:).^2));
end

function X=fracshift(X,dx,dy,dz)
Xhat=fftn(X);
[N1,N2,N3]=size(X);
%N=7:[0,1,2, 3,-3,-2,-1] or N=8:[0,1,2,3, -4,-3,-2,-1]
M1a=floor((N1-1)/2); M1b=floor(N1/2);
kxrange=[0:M1a,-M1b:-1];
M2a=floor((N2-1)/2); M2b=floor(N2/2);
kyrange=[0:M2a,-M2b:-1];
M3a=floor((N3-1)/2); M3b=floor(N3/2);
kzrange=[0:M3a,-M3b:-1];
[kx,ky,kz]=ndgrid(kxrange,kyrange,kzrange);
theta=(kx*dx/N1+ky*dy/N2+kz*dz/N3)*2*pi;
Xhat=Xhat.*exp(i*theta);
X=real(ifftn(Xhat));
end
