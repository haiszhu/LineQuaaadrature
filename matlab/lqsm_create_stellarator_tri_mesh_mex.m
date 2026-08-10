function [sx, snx, sw, rts, rps, ier] = lqsm_create_stellarator_tri_mesh_mex(mp, np, p, nfp, nmode, mn, rc, zs, nchart, charts, ntri, sx, snx, sw, rts, rps, ier)
mp = double(mp); np = double(np); p = double(p);
nfp = double(nfp); nmode = double(nmode); nchart = double(nchart);
mn = double(reshape(mn, 2, nmode));
rc = double(rc(:)); zs = double(zs(:));
charts = double(reshape(charts, 6, nchart));
ntri = double(ntri); hdim = p*(p+1)/2; nsrc = ntri*hdim;
sx = double(reshape(sx, 3, nsrc)); snx = double(reshape(snx, 3, nsrc));
sw = double(sw(:)); rts = double(reshape(rts, 3, nsrc));
rps = double(reshape(rps, 3, nsrc)); ier = double(ier);
mex_id_ = 'lqsm_create_stellarator_tri_mesh_mex(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[xx], c i double[x], c i double[x], c i int64_t[x], c i double[xx], c i int64_t[x], c io double[xx], c io double[xx], c io double[x], c io double[xx], c io double[xx], c io int64_t[x])';
[sx, snx, sw, rts, rps, ier] = LineQuaaadrature_mex(mex_id_, mp, np, p, nfp, nmode, mn, rc, zs, nchart, charts, ntri, sx, snx, sw, rts, rps, ier, 1, 1, 1, 1, 1, 2, nmode, nmode, nmode, 1, 6, nchart, 1, 3, nsrc, 3, nsrc, nsrc, 3, nsrc, 3, nsrc, 1);
end
