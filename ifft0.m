function f=ifft0(fhat)
f=fftshift(ifftn(fftshift(fhat)));
end