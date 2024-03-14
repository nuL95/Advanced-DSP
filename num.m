N = 50;
M = N/2;
fs = 24e3;
[c w] = cospolyvec(M);
c = c.';
C = c*c.';
Q = zeros(M+1,M+1);
F = [0 4000 5000 7000 7500 8500 9000 12000];
Fn = F * pi / (fs/2);

ramp = (12/pi)*w-7.5;



%intervals
interval = 1024;
    w1 = 0:pi/interval:pi/3-pi/interval;
    w2 = 5*pi/12:pi/interval:7*pi/12 - pi/interval;
    w3 = 7.5*pi/12:pi/interval:8.5*pi/12-pi/interval;
    w4 = 3*pi/4:pi/interval:pi-pi/interval;
for ii = 1:M+1
    for jj = 1:M+1
        Q(ii,jj) = trapz(w1, cos(w1*(ii-1)).*cos(w1*(jj-1))) + trapz(w2, cos(w2*(ii-1)).*cos(w2*(jj-1))) + trapz(w3, cos(w3*(ii-1)).*cos(w3*(jj-1))) + trapz(w4, cos(w4*(ii-1)).*cos(w4*(jj-1)));
    end
end
%this looks really sloppy, the problem is that all the intervals are of
%different length , i think if i used pointers i could kind of 'fix' this
%but i dont want to think about that right now

d = double(pi/3 + int(ramp^2, 7.5*pi/12,8.5*pi/12));
p = double(int(c, 0, pi/3) + int(ramp.*c,7.5*pi/12,8.5*pi/12));
Qt = [Q p;p' d];
e1 = double(subs(c, w, 4.5*pi/12));
e2 = double(subs(c, w, 8*pi/12));
E = [e1 e2]';
B = null(E);
[evec, ev] = eig(Qt)

e0 = evec(:,1);

a = e0/-e0(end);

a = a(1:end-1);

h = zeros(1,N+1);
h(M+1) = a(1)
af = fliplr(a(2:end)')/2;
h(1:M) = af(1:end);
h(M+2:N+1) = a(2:end)/2
wplot = 0:12e3/512:12e3-1/512;
Hfr = freqz(h);
plot(wplot, abs(Hfr))
function [cpv, sym1] = cospolyvec(M)
    syms w
    cpv = sym('x',[1 M+1]);
    for ii = 1:M+1
        cpv(ii) = [cos((ii-1)*w)];
    end
    sym1 = w;
end