%% 
% Declaration of variables
N = 50;
M = N/2;
fs = 24e3;
p = zeros(M+1,1);
Q = zeros(M+1,M+1);
F = [0 4000 5000 7000 7500 8500 9000 12000];
Fnotch = [4500 8000];
Fnotch = pi*Fnotch/(fs/2);
E = zeros(length(Fnotch), M+2);
Fn = F * pi / (fs/2);
int_grid = 1024;
W = zeros(length(Fn)/2,int_grid);
%% 
% Calculation of p, d, and Q
for ii = 1:length(Fn)/2
    W(ii,:) = linspace(Fn(2*ii-1), Fn(2*ii), int_grid);
end

for ii = 1:M+1
    for jj = 1:M+1
        for kk = 1:length(Fn)/2
            Q(ii,jj) = Q(ii,jj) + trapz(W(kk,:), cos(W(kk,:)*(ii-1)).*cos(W(kk,:)*(jj-1)));
        end
    end
end

d = trapz(W(3,:), ((12/pi)*W(3,:)-7.5).^2) + 1*Fn(2);

for ii = 1:M+1
    p(ii) = trapz(W(1,:), cos((ii-1)*W(1,:))) + trapz(W(3,:), ((12/pi)*W(3,:)-7.5).*cos((ii-1)*W(3,:)));
end

Qt = [Q p;p' d];
for ii = 1:M+1
    for jj = 1:length(Fnotch)
        E(jj,ii) = cos((ii-1)*Fnotch(jj));
    end
end
%%
% Minimization of Rayleigh quotient by applying Courant Fischer Theorem
% With linear constraints for notch filter 
B = null(E);
Qaug = B'*Qt*B;
[evec, ev] = eig(Qaug);
[mag, ind] = min(nonzeros(ev));
e0 = evec(:,ind);
a = B*e0;
a = a/-a(end);
a = a(1:end-1);
%%
%Retrieving filter coefficients from cosine polynomial coefficients
h = zeros(1,N+1);
h(M+1) = a(1);
af = fliplr(a(2:end)')/2;
h(1:M) = af(1:end);
h(M+2:N+1) = a(2:end)/2;
wplot = 0:12e3/512:12e3-1/512;
Hfr = freqz(h);
%%
%Plotting magnitude response
figure(2)
plot(wplot, abs(Hfr))
title('Magnitude Response of Eigenfilter With Notches')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
%%