interval = 0.00001;
Fs = 1/interval;
t = 0:interval:2;
y = chirp(t,100,1,150);
plot(t,y);
f=fft(y, pow2(nextpow2(length(y))));
NFFT = pow2(nextpow2(length(y)));
freq = Fs/2*linspace(0,1,NFFT/2+1);
figure
plot(freq,2*abs(f(1:NFFT/2+1)),'r');