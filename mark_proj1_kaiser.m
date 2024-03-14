%mark_proj1_kaiser.m
%Filter Characteristics
fpass = 4e3;
fstop = 4.5e3;
dp = 0.1;
ds = 0.05;
fs = 20e3;
A = [1 0]
w = 0:1/512:1-1/512;
[n, Wn, beta, ftype] = kaiserord([fpass fstop], A, [dp ds], fs)

%Ensure that the filter is type 1 by making the order even if its odd.
if mod(n, 2) == 1
    n = n + 1;
end

KFIRF = fir1(n, Wn, kaiser(n+1,beta),ftype, 'noscale');
[FAmp, FPh] = freqz(KFIRF);
figure(1)
plot(KFIRF)
title('Impulse Response')
xlabel('Sample')
ylabel('Amplitude')
figure(2)
freqz(KFIRF)
title('Frequency Response')
figure(3)
plot(w,abs(FAmp))
title('Magnitude Response')
xlabel('normalized frequency x pi / sample')
ylabel('Magnitude')