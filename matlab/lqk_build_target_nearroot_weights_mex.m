function [funvals0, weights, troot, accepted_i, I_local] = lqk_build_target_nearroot_weights_mex(nquad, tgl, wgl, legmat, xj, yj, zj, spj, stauj, xjhat, yjhat, zjhat, rho, xtk, ytk, ztk, fptr_int, kdata, funvals0, weights, troot, accepted_i, I_local, kernel_id, adaptive_fallback_i)
nquad = double(nquad);
tgl = double(tgl(:));
wgl = double(wgl(:));
legmat = double(reshape(legmat, nquad, nquad));
xj = double(xj(:));
yj = double(yj(:));
zj = double(zj(:));
spj = double(spj(:));
stauj = double(reshape(stauj, 3, nquad));
xjhat = double(xjhat(:));
yjhat = double(yjhat(:));
zjhat = double(zjhat(:));
rho = double(rho);
xtk = double(xtk);
ytk = double(ytk);
ztk = double(ztk);
fptr_int = double(fptr_int);
kdata = double(kdata(:));
funvals0 = double(funvals0(:));
weights = double(weights(:));
troot = complex(troot);
accepted_i = double(accepted_i);
I_local = double(I_local);
kernel_id = double(kernel_id);
adaptive_fallback_i = double(adaptive_fallback_i);
mex_id_ = 'lqk_build_target_nearroot_weights_mex(c i int64_t[x], c i double[x], c i double[x], c i double[xx], c i double[x], c i double[x], c i double[x], c i double[x], c i double[xx], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i int64_t[x], c i double[x], c io double[x], c io double[x], c io dcomplex[x], c io int64_t[x], c io double[x], c i int64_t[x], c i int64_t[x])';
[funvals0, weights, troot, accepted_i, I_local] = LineQuaaadrature_mex(mex_id_, nquad, tgl, wgl, legmat, xj, yj, zj, spj, stauj, xjhat, yjhat, zjhat, rho, xtk, ytk, ztk, fptr_int, kdata, funvals0, weights, troot, accepted_i, I_local, kernel_id, adaptive_fallback_i, 1, nquad, nquad, nquad, nquad, nquad, nquad, nquad, nquad, 3, nquad, nquad, nquad, nquad, 1, 1, 1, 1, 1, 3, nquad, nquad, 1, 1, 1, 1, 1);
end

% --------------------------------------------------------------------------
