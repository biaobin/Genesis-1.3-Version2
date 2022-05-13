function [l,plunit]=spectrim(s,po,phase,N,l0)
%calculate the spectrum of FEL pulse
%change 
 %l0 nm


p    = po;	%[W] 
phi  = phase;	%[rad] 

c0   = 299.792458;	%[nm/fs]
w0   = 2*pi*c0/l0;
dt   = (s(end)-s(1))*1e9/c0/N;
dw   = 2*pi/N/dt;
w    = (-N/2:N/2-1)*dw;
absw = w+w0;

et   = sqrt(p).*exp(1i*phi);
pt   = et.*conj(et);
ew   = fftshift(fft(et,N));
pw   = ew.*conj(ew);

temp = find(absw>0);
zp   = temp(1);
l    = 2*pi*c0./absw(end:-1:zp);
pl   = pw(end:-1:zp);
plunit=pl/max(pl);


