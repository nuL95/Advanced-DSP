N = 50;
M = N/2;
[c w] = cospolyvec(M);
c = c.';
C = c*c.';
Q = zeros(M+1,M+1);
fs = 24e3;
ramp = (12/pi)*w-7.5;

Cint = int(C, w);
d = double(pi/3 + int(ramp^2, 7.5*pi/12,8.5*pi/12));
p = double(int(c, 0, pi/3) + int(ramp.*c,7.5*pi/12,8.5*pi/12));
Q = double(subs(Cint, w, pi/3)-subs(Cint, w, 0) + subs(Cint,w, 7*pi/12)-subs(Cint,w,5*pi/12) + subs(Cint,w, 8.5*pi/12)-subs(Cint,w,7.5*pi/12) + subs(Cint, w, pi)-subs(Cint, w,3*pi/4));
Qt = [Q p;p' d];
cd = diff(c, w);
cdd = diff(cd, w);
e1 = double(subs(c, w, 4.5*pi/12));
e2 = double(subs(c, w, 8*pi/12));

v = zeros(1,6);
E = [e1' v(1);e2' v(2)];
B = [null(E)];
Qaug = B'*Qt*B;
[evec, ev] = eig(Qaug)
[mag, ind] = min(nonzeros(ev))
e0 = evec(:,ind);
a = B*e0;

a = a/-a(end);

a = a(1:end-1);

h = zeros(1,N+1);
h(M+1) = a(1)
af = fliplr(a(2:end)')/2;
h(1:M) = af(1:end);
h(M+2:N+1) = a(2:end)/2
wplot = 0:12e3/512:12e3-1/512;
Hfr = freqz(h);
figure(3)
plot(wplot, abs(Hfr))
function [cpv, sym1] = cospolyvec(M)
    syms w
    cpv = sym('x',[1 M+1]);
    for ii = 1:M+1
        cpv(ii) = [cos((ii-1)*w)];
    end
    sym1 = w;
end