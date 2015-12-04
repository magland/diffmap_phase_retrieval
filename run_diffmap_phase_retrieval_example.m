function run_diffmap_phase_retrieval_example(example_opts,opts)
if (nargin<1)
    diffmap_phase_retrieval_example1;
    return;
end;

f=create_example(example_opts);
f_noise=f+randn(size(f))*example_opts.noise_level;
u=abs(fft2b(f_noise));
opts.reference=f;

recon=diffmap_phase_retrieval(u,opts);
recon=real(ifft0(u.*exp(i*angle(fft0(recon))))); %enforce the magnitude constraint on the last step -- this makes the resid calculation fair

resid=compute_residual(recon,f);

fig=figure; set(fig,'position',[50,50,1400,500]);
maxval=max(f(:));
subplot(1,3,1);
imagescb(f,[0,maxval]);
subplot(1,3,2);
imagescb(recon,[0,maxval]);
subplot(1,3,3);
imagescb(abs(f-recon),[0,maxval]);
title(sprintf('Resid = %g, %g\n',resid));

end

function resid=compute_residual(f,ref)
resid=sqrt(sum((f(:)-ref(:)).^2))/sqrt(sum(ref(:).^2));
end

function f=create_example(opts)
oversamp=opts.oversamp;
N1=opts.N1; N2=opts.N2; N3=1;

N1b=N1*oversamp;
N2b=N2*oversamp;
N3b=1;
[xx,yy,zz]=ndgrid(linspace(-oversamp,oversamp,N1b),linspace(-oversamp,oversamp,N2b),0);
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
end;
end

function imagescb(X,range)
imagesc(X); colormap('gray');
if (nargin>=2)
    caxis(range);
end;

end

function Y=fft2b(X)
Y=fftshift(fft2(fftshift(X)));
end

function Y=ifft2b(X)
Y=fftshift(ifft2(fftshift(X)));
end

