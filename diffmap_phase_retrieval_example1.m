function Y=diffmap_phase_retrieval_example1(X)

opts.max_iterations=200;
opts.diffmap_method='BA';
example_opts.example_num=1;
example_opts.N1=128;
example_opts.N2=128;
example_opts.oversamp=2;
noise_level=0.05;

f=create_example(example_opts);
f_noise=f+randn(size(f))*noise_level;
u=abs(fft2b(f_noise));
opts.reference=f;

recon=diffmap_phase_retrieval(u,opts);

fig=figure; set(fig,'position',[50,50,1000,500]);
subplot(1,3,1);
imagescb(f,[0,1]);
subplot(1,3,2);
imagescb(recon,[0,1]);
subplot(1,3,3);
imagescb(abs(f-recon),[0,1]);

resid=compute_residual(recon,f)
best_possible_resid=compute_residual(f_noise,f)

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
    for kk=1:20
        cc=(rand(3,1)*2-1)*0.7;
        rr=(rand*2-1)*0.2;
        if (N3b==1) cc(3)=0; end;
        R=sqrt((xx-cc(1)).^2+(yy-cc(2)).^2+(zz-cc(3)).^2);
        f=f+...
            exp(-(R.^2)/(rr/2)^2).*...
            (R.^2<=rr^2);
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

