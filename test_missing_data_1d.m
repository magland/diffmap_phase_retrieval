function test_missing_data_1d

close all;

N=40;
oversamp=3;

k_locations=linspace(-N/2,N/2-1/oversamp,N*oversamp);
x_locations=linspace(-0.5,0.5-1/N,N);

k_locations=k_locations(find(abs(k_locations)>=3));

[kk,xx]=ndgrid(k_locations,x_locations);

figure; plot(k_locations,'b+');

AA=exp(2*pi*i*xx.*kk);
cond(AA)
%AA(1,:)'*AA(2,:)
[U,D,V]=svd(AA,'econ');
diag(D)'
% figure; plot(real(U));
% figure; plot(real(U'));
% figure; plot(real(V));
%figure; plot(real(U(:,36:40)));

f=create_a_function(N)';
figure; plot(f);
fhat=AA*f;
D_dagger=inv(D).*(D>=max(D(:))*1e-2);
diag(D_dagger)'
AA_dagger=V*D_dagger*U';
frecon=real(AA_dagger*fhat);
figure; plot(frecon);

figure; plot(real(V(:,end-1:end)));
figure; plot(real(U(:,end-1:end)));

end

function f=create_a_function(N)


seed=1;
rng(seed);

example_opts.example_num=4;
example_opts.num_disks=50;
example_opts.N1=N;
example_opts.N2=N;
example_opts.oversamp=1;
example_opts.k=2;
example_opts.offset=0;
noise_level=0.00;
relative_support_sizes=0.6:0.025:1.4;

alg_opts.max_iterations=100;
alg_opts.method='BA_HIO';
alg_opts.positivity=0;
alg_opts.support=1;
alg_opts.beta=0.8;
alg_opts.num_restarts=20;

[f,support_mask]=create_example(example_opts);
fnoisy=f+randn(size(f))*noise_level;
figure; imagesc(f');
f_proj=sum(f,1);

f=f_proj;

end

function stuff_to_hold;
fnoisy_proj=sum(fnoisy,1);
true_support_range=[min(find(f_proj)),max(find(f_proj))];
figure; plot(fnoisy_proj);
uproj=abs(fft1b(fnoisy_proj));

resids=[];
for ii=1:length(relative_support_sizes)
    relative_support_size=relative_support_sizes(ii);
    support_range=(true_support_range-mean(true_support_range))*relative_support_size+mean(true_support_range);
    aaa=(1:length(f_proj));
    alg_opts.support_mask=(support_range(1)<=aaa).*(aaa<=support_range(2));

    [fproj_recon,info]=do_recon_1d(uproj,alg_opts);
    figure;
    plot(fnoisy_proj,'k'); hold on;
    plot(fproj_recon,'b'); hold on;
    plot(pi_A_1d(fproj_recon,uproj,alg_opts),'r');

    figure; semilogy(info.resids);
    set(gcf,'position',[50,50,800,800]);
    fprintf('rel. support size=%g, resid=%g\n',relative_support_size,info.resid);
    
    resids(ii)=info.resid;
end;
    
figure; plot(relative_support_sizes,resids);
xlabel('Relative support constraint widths');
ylabel('Residuals');
    
end

function [recon,info]=do_recon_1d(u,opts)

if opts.num_restarts>1
    best_recon=[];
    best_info.resid=inf;
    num_restarts=opts.num_restarts;
    opts.num_restarts=1;
    for ct=1:num_restarts
        fprintf('.');
        [recon0,info0]=do_recon_1d(u,opts);
        if (info0.resid<best_info.resid)
            best_recon=recon0;
            best_info=info0;
        end;
    end;
    fprintf('\n');
    recon=best_recon;
    info=best_info;
    return;
end;

initial_phase=rand(size(u))*2*pi;
f0=real(ifft1b(u.*exp(i*initial_phase)));

info.resids=[];
beta=opts.beta;

for it=1:opts.max_iterations
    if strcmp(opts.method,'AB')
        pi_B_f0=pi_B_1d(f0,opts);
        proj1=pi_A_1d(2*pi_B_f0-f0,u,opts);
        proj2=pi_B_f0;
        f1=f0+proj1-proj2;
        f0=f1;
        recon=pi_B_1d(f0,opts);
    elseif strcmp(opts.method,'AB_HIO')
        pi_B_f0=pi_B_1d(f0,opts);
        proj1=pi_A_1d((1+beta)*pi_B_f0-f0,u,opts);
        proj2=pi_B_f0;
        f1=f0+proj1-beta*proj2;
        f0=f1;
        recon=pi_B_1d(f0,opts);
    elseif strcmp(opts.method,'BA')
        pi_A_f0=pi_A_1d(f0,u,opts);
        proj1=pi_B_1d(2*pi_A_f0-f0,opts);
        proj2=pi_A_f0;
        f1=f0+proj1-proj2;
        f0=f1;
        recon=pi_A_1d(f0,u,opts);
    elseif strcmp(opts.method,'BA_HIO')
        pi_A_f0=pi_A_1d(f0,u,opts);
        proj1=pi_B_1d((1+beta)*pi_A_f0-f0,opts);
        proj2=pi_A_f0;
        f1=f0+proj1-beta*proj2;
        f0=f1;
        recon=pi_A_1d(f0,u,opts);
    end;
    info.resids(end+1)=compute_residual_1d(pi_A_1d(recon,u,opts),pi_B_1d(pi_A_1d(recon,u,opts),opts));
end;

info.resid=info.resids(end);

end

function [f,support_mask]=create_example(opts)
oversamp=opts.oversamp;
N1=opts.N1; N2=opts.N2; N3=1;

N1b=N1*oversamp;
N2b=N2*oversamp;
N3b=1;
[xx,yy,zz]=ndgrid(linspace(-oversamp,oversamp,N1b),linspace(-oversamp,oversamp,N2b),0);
support_mask=(abs(xx)<=1.1).*(abs(yy)<=1.1).*(abs(zz)<=1.2);
if (opts.example_num==1)
    f=zeros(size(xx));
    for kk=1:opts.num_disks
        cc=(rand(3,1)*2-1)*0.7;
        rr=(rand*2-1)*0.2;
        if (N3b==1) cc(3)=0; end;
        R=sqrt((xx-cc(1)).^2+(yy-cc(2)).^2+(zz-cc(3)).^2);
        f=f+...
            exp(-(R.^2)/(rr/2)^2).*...
            (R.^2<=rr^2);
    end;
elseif (opts.example_num==2)
    f=zeros(size(xx));
    for kk=1:opts.num_disks
        cc=(rand(3,1)*2-1)*0.7;
        rr=(rand*2-1)*0.2;
        if (N3b==1) cc(3)=0; end;
        R=sqrt((xx-cc(1)).^2+(yy-cc(2)).^2+(zz-cc(3)).^2);
        f=f+(R.^2<=rr^2);
    end;
elseif (opts.example_num==3)
    f=zeros(size(xx));
    for kk=1:opts.num_disks
        cc=(rand(3,1)*2-1)*0.7;
        rr=(rand*2-1)*0.2;
        if (N3b==1) cc(3)=0; end;
        R=sqrt((xx-cc(1)).^2+(yy-cc(2)).^2+(zz-cc(3)).^2);
        f=f+exp(-(R.^2)/(rr/2)^2);
    end;
elseif (opts.example_num==4)
    f=zeros(size(xx));
    for kk=1:opts.num_disks
        cc=(rand(3,1)*2-1)*0.7;
        rr=abs((rand*2-1)*0.2);
        if (N3b==1) cc(3)=0; end;
        R=sqrt((xx-cc(1)).^2+(yy-cc(2)).^2+(zz-cc(3)).^2);
        k=opts.k;
        ampl=rand;
        f=f+ampl*(R/rr<1).*(1-(R/rr).^2).^k;
    end;
    f=f+opts.offset;
    f=f/max(f(:));
end;
f=f.*support_mask;
%support_mask=(abs(f)>1e-8);
end

function Y=fft2b(X)
Y=fftshift(fft2(fftshift(X)));
end

function Y=ifft2b(X)
Y=fftshift(ifft2(fftshift(X)));
end

function Y=fft1b(X)
Y=fftshift(fft(fftshift(X)));
end

function Y=ifft1b(X)
Y=fftshift(ifft(fftshift(X)));
end

function f2=pi_A_1d(f0,u,opts)
f1=fft1b(f0);
f1=u.*exp(i*angle(f1));
f1=real(ifft2b(f1));
f2=f1;
end

function f0=pi_B_1d(f0,opts)
if opts.positivity
    f0=f0.*(f0>=0);
end;
if opts.support
    f0=f0.*opts.support_mask;
end;
end

function resid=compute_residual_1d(f,ref)
resid=sqrt(sum((f(:)-ref(:)).^2))/sqrt(sum(ref(:).^2));
end
