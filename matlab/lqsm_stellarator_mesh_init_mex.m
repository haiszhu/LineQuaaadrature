function [charts, nchart, ntri, ier] = lqsm_stellarator_mesh_init_mex(mp, np, p, nfp, nmode, mn, rc, zs, restol, cap, charts, nchart, ntri, ier)
mp = double(mp); np = double(np); p = double(p);
nfp = double(nfp); nmode = double(nmode); cap = double(cap);
mn = double(reshape(mn, 2, nmode));
rc = double(rc(:)); zs = double(zs(:)); restol = double(restol);
charts = double(reshape(charts, 6, cap));
nchart = double(nchart); ntri = double(ntri); ier = double(ier);
mex_id_ = 'lqsm_stellarator_mesh_init_mex(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[xx], c i double[x], c i double[x], c i double[x], c i int64_t[x], c io double[xx], c io int64_t[x], c io int64_t[x], c io int64_t[x])';
[charts, nchart, ntri, ier] = LineQuaaadrature_mex(mex_id_, mp, np, p, nfp, nmode, mn, rc, zs, restol, cap, charts, nchart, ntri, ier, 1, 1, 1, 1, 1, 2, nmode, nmode, nmode, 1, 1, 6, cap, 1, 1, 1);
end
