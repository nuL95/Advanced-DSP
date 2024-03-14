%mark_proj1_equiripple
%Filter Characteristics
fpass = 4e3;
fstop = 4.5e3;
ds = 0.05;
dp = 0.1;
fs = 20e3;
A = [1 0];
omega = 0:1/512:1-1/512
[n, f0, a0, w] = firpmord([fpass fstop], A, [dp ds], fs);

%Ensure that the filter is type 1 by making the order even if its odd.
if mod(n, 2) == 1
    n = n + 1;
end

FIRF_equiripple = firpm(n, f0, a0, w);
figure(1)
plot(FIRF_equiripple)
title('Impulse Response')
xlabel('Sample')
ylabel('Amplitude')
figure(2)
freqz(FIRF_equiripple)
title('Frequency Response')


[FeqAm, Feqf] = freqz(FIRF_equiripple);
FeqMag = abs(FeqAm);

fcut_index = 2+round(512*fstop/(fs/2));
ds_new = max(FeqMag(fcut_index:end));

mpf = FIRF_equiripple;
mpf(n/2+1) = mpf(n/2+1)+ds_new;
figure(3)
zplane(mpf)
title('zero/pole plot of modified equiripple filter')

mpfd = FIRF_equiripple;
mpfd(n/2+1) = mpf(n/2+1)+ds_new*1.01;
figure(4)
zplane(mpfd)
title('zero/pole plot of the equiripple with extra large ds')

zeros_d = sort(roots(mpfd));
zeros_IUC = abs(zeros_d) < 1;
zeros = sort(roots(mpf));
zeros_keep = nonzeros(zeros.*zeros_IUC);
mpf = poly(zeros_keep);
mpfFR = freqz(mpf);
m_pass = mpfFR(1:round(512*fpass/(fs/2)));
norm_factor = mean(abs(m_pass));
mpf = mpf/norm_factor;
mpfFR = freqz(mpf);

figure(5)
zplane(mpf);
title('zero/pole plot of the minimum phase filter');
figure(6)
hold on
plot(omega, abs(mpfFR));
plot(omega, FeqMag);
title('Magnitude Response of Equiripple and Minimum Phase Filter');
xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude');
legend('Blue - Minimum Phase', 'Orange - Equiripple')
hold off
figure(7)
freqz(mpf);
title('Frequency Response of Minimum Phase Filter');