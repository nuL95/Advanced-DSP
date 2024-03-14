%mark_proj1_equiripple
%Filter Characteristics
fpass = 4e3;
fstop = 4.5e3;
ds = 0.05;
dp = 0.1;
fs = 20e3;
A = [1 0];

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
freqz(FIRF_equiripple,N)
title('Frequency Response')


[FeqAm, Feqf] = freqz(FIRF_equiripple);
FeqMag = abs(FeqAm);

fcut_index = 2+round(512*fstop/(fs/2));
ds_new = max(FeqMag(fcut_index:end));

mpf = FIRF_equiripple;
mpf(n/2+1) = mpf(n/2+1)+ds_new;
mpfRoots = roots(mpf);
mpfRoots = sort(mpfRoots);
mpfRootsMag = abs(mpfRoots);

tol = 1e-2;
IUCroots = abs(mpfRoots) < 1-tol;
UCroots = abs(mpfRoots) > 1-tol & abs(mpfRoots) < 1+tol;

IUCzeros = nonzeros(IUCroots.*mpfRoots);
UCzeros = nonzeros(UCroots.*mpfRoots);
UCzeros = sort(UCzeros, 'ComparisonMethod', 'real');
if(mod(length(UCzeros),4) ~= 0)
    oddzero = UCzeros(1);
    UCzeros = UCzeros(3:end);

end

zplane(poly(UCzeros))
ics = reshape(UCzeros, 2, length(UCzeros)/2);
ics = ics(:,1:2:end);
UCzeros = reshape(ics, 1, length(UCzeros)/2)';
zerosKeep = [oddzero' UCzeros' IUCzeros'];

mpf = poly(zerosKeep);
zplane(mpf)
mpfFR = freqz(mpf);
m_pass = mpfFR(1:round(512*fpass/(fs/2)));
norm_factor = mean(abs(m_pass));
mpf = mpf/norm_factor;
mpfFR = mpfFR/norm_factor;
figure(4)
plot(abs(mpfFR))
hold on
plot(FeqMag)
hold off