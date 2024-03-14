%% 
% Declaration of variables
N = 50;
M = N/2;
fs = 24e3;
p = zeros(M+1,1);
Q = zeros(M+1,M+1);
F = [0 4000 5000 12000];
Fnotch = [4500 8000];
Fnotch = pi*Fnotch/(fs/2);
E = zeros(length(Fnotch)*2, M+2);
Fn = F * pi / (fs/2);
int_grid = 1024;
W = zeros(length(Fn)/2,int_grid);
Qtest = zeros(M+1,M+1);
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
d = 1*Fn(2);

for ii = 1:M+1
    p(ii) = trapz(W(1,:), cos((ii-1)*W(1,:))) + trapz(W(3,:));
end
%%
% Minimization of Rayleigh quotient by applying Courant Fischer Theorem
% Formulation of linear constraint matrix, and completing the square to
% form Qt
Qt = [Q p;p' d];

for ii = 1:M+1
    for jj = 1:length(Fnotch)
        E(jj,ii) = cos((ii-1)*Fnotch(jj));
%         E(jj+2,ii) = (ii-1)*sin((ii-1)*Fnotch(jj));
    end
end
%eigenfilter
B = null(E);
[evec, ev] = eig(Qt);
[mag, ind] = min(nonzeros(ev));
e0 = evec(:,ind);
a = e0/-e0(end);
a = a(1:end-1);
%eigen filter with notch
Qaug = B'*Qt*B;
[evec, ev] = eig(Qaug);
[mag, ind] = min(nonzeros(ev));
e0 = evec(:,ind);
an = B*e0;
an = an/-an(end);
an = an(1:end-1);

%%
%Retrieving filter coefficients from cosine polynomial coefficients
%eigenfilter
h = zeros(1,N+1);
h(M+1) = a(1);
af = fliplr(a(2:end)')/2;
h(1:M) = af(1:end);
h(M+2:N+1) = a(2:end)/2;
wplot = 0:12e3/512:12e3-1/512;
Hfr = freqz(h);
%eigenfilter with notch
hn = zeros(1,N+1);
hn(M+1) = an(1);
afn = fliplr(an(2:end)')/2;
hn(1:M) = afn(1:end);
hn(M+2:N+1) = an(2:end)/2;
Hfrn = freqz(hn);
%%
%Plotting magnitude response
figure(1)
plot(wplot, abs(Hfr))
title('Magnitude Response of Eigenfilter')
xlabel('Frequency (Hz)')
ylabel('Magnitude')

figure(2)
plot(wplot, abs(Hfrn))
title('Magnitude Response of Eigenfilter With Notch')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
%%

D = zeros(1,512);
D(1:171) = 1;
f = 0:42;
D(321:363) = f/42;
figure(3)
plot(wplot, D)
hold on
plot(wplot, abs(Hfr))
plot(wplot, abs(Hfrn))
title('Ideal and Eigenfilters On Same Plot')
ylabel('Magnitude')
xlabel('Frequency (Hz')
legend('desired', 'designed', 'designed with notches')