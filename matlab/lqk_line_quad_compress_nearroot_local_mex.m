function [funvals0, sxbdw] = lqk_line_quad_compress_nearroot_local_mex(root_ok, nquad, npan, tgl, wgl, legmat, xhat, yhat, zhat, t0, tpan, wpan, target_loc, fptr_int, kdata, integrand0_up, funvals0, sxbdw)
root_ok = double(root_ok);
nquad = double(nquad);
npan = double(npan);
nup = nquad*npan;
tgl = double(tgl(:));
wgl = double(wgl(:));
legmat = double(reshape(legmat, nquad, nquad));
xhat = double(xhat(:));
yhat = double(yhat(:));
zhat = double(zhat(:));
t0 = double(t0);
tpan = double(reshape(tpan, nquad, npan));
wpan = double(reshape(wpan, nquad, npan));
target_loc = double(target_loc(:));
fptr_int = double(fptr_int);
kdata = double(kdata(:));
integrand0_up = double(reshape(integrand0_up, nup, 1));
funvals0 = double(reshape(funvals0, nquad, 1));
sxbdw = double(reshape(sxbdw, nquad, 1));
mex_id_ = 'lqk_line_quad_compress_nearroot_local_mex(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[xx], c i double[x], c i double[x], c i double[x], c i double[x], c i double[xx], c i double[xx], c i double[x], c i int64_t[x], c i double[x], c i double[x], c io double[x], c io double[x])';
[funvals0, sxbdw] = LineQuaaadrature_mex(mex_id_, root_ok, nquad, npan, tgl, wgl, legmat, xhat, yhat, zhat, t0, tpan, wpan, target_loc, fptr_int, kdata, integrand0_up, funvals0, sxbdw, 1, 1, 1, nquad, nquad, nquad, nquad, nquad, nquad, nquad, 1, nquad, npan, nquad, npan, 3, 1, 3, nup, nquad, nquad);
end



% --------------------------------------------------------------------------
