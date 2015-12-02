function fhat=fft0(f)
fhat=fftshift(fftn(fftshift(f)));
end

