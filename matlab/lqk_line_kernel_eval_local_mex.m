function [funvals_local, integrand0_up, I_local, kval] = lqk_line_kernel_eval_local_mex(nquad, npan, xdisp, ydisp, zdisp, stangpan, sppan, dswpan, target_loc, fptr_int, kdata, funvals_local, integrand0_up, I_local, kval)
nquad = double(nquad);
npan = double(npan);
nquad3 = 3*nquad;
nup = nquad*npan;
xdisp = double(reshape(xdisp, nquad, npan));
ydisp = double(reshape(ydisp, nquad, npan));
zdisp = double(reshape(zdisp, nquad, npan));
stangpan = double(reshape(stangpan, nquad3, npan));
sppan = double(reshape(sppan, nquad, npan));
dswpan = double(reshape(dswpan, nquad, npan));
target_loc = double(target_loc(:));
fptr_int = double(fptr_int);
kdata = double(kdata(:));
funvals_local = double(reshape(funvals_local, nquad, npan));
integrand0_up = double(reshape(integrand0_up, nup, 1));
I_local = double(I_local);
kval = double(reshape(kval, nquad, npan));
mex_id_ = 'lqk_line_kernel_eval_local_mex(c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[xx], c i double[xx], c i double[xx], c i double[xx], c i double[x], c i int64_t[x], c i double[x], c io double[xx], c io double[x], c io double[x], c io double[xx])';
[funvals_local, integrand0_up, I_local, kval] = LineQuaaadrature_mex(mex_id_, nquad, npan, xdisp, ydisp, zdisp, stangpan, sppan, dswpan, target_loc, fptr_int, kdata, funvals_local, integrand0_up, I_local, kval, 1, 1, nquad, npan, nquad, npan, nquad, npan, nquad3, npan, nquad, npan, nquad, npan, 3, 1, 3, nquad, npan, nup, 1, nquad, npan);
stangpan = double(reshape(stangpan, 3, nquad, npan));
integrand0_up = double(reshape(integrand0_up, nup, 1));
end

% --------------------------------------------------------------------------
