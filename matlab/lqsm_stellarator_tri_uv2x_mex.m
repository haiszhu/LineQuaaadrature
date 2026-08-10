function [x, ier] = lqsm_stellarator_tri_uv2x_mex(mp, np, p, nfp, nmode, mn, rc, zs, nchart, charts, itri, nuv, uv, x, ier)
mp = double(mp); np = double(np); p = double(p);
nfp = double(nfp); nmode = double(nmode); nchart = double(nchart);
mn = double(reshape(mn, 2, nmode));
rc = double(rc(:)); zs = double(zs(:));
charts = double(reshape(charts, 6, nchart));
itri = double(itri); nuv = double(nuv); uv = double(reshape(uv, 2, nuv));
x = double(reshape(x, 3, nuv)); ier = double(ier);
mex_id_ = 'lqsm_stellarator_tri_uv2x_mex(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[xx], c i double[x], c i double[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c io double[xx], c io int64_t[x])';
[x, ier] = LineQuaaadrature_mex(mex_id_, mp, np, p, nfp, nmode, mn, rc, zs, nchart, charts, itri, nuv, uv, x, ier, 1, 1, 1, 1, 1, 2, nmode, nmode, nmode, 1, 6, nchart, 1, 1, 2, nuv, 3, nuv, 1);
end

% --------------------------------------------------------------------------
