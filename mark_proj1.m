zplane(mp_firf)
x = MPFIR_roots_UC
xabs = abs(x);
abs(MPFIR_roots_UC);

for jj = 1:length(MPFIR_roots_UC)/4
    mpuc(2*jj-1:2*jj) = MPFIR_roots_UC(4*jj-3:4*jj-2);
end

R = conv(poly(mpuc), poly(MPFIR_roots_OUC))
zplane(R)