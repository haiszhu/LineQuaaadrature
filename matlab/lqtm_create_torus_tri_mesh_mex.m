function [sx, snx, sw, rts, rps, ier] = lqtm_create_torus_tri_mesh_mex(mp, np, p, radii, scales, nosc, orig, ntri, sx, snx, sw, rts, rps, ier)
mp = double(mp); np = double(np); p = double(p);
radii = double(reshape(radii, 3, 1));
scales = double(reshape(scales, 3, 1));
nosc = double(nosc); orig = double(orig); ntri = double(ntri);
hdim = p*(p+1)/2; nsrc = ntri*hdim;
sx = double(reshape(sx, 3, nsrc)); snx = double(reshape(snx, 3, nsrc));
sw = double(sw(:)); rts = double(reshape(rts, 3, nsrc));
rps = double(reshape(rps, 3, nsrc)); ier = double(ier);
mex_id_ = 'lqtm_create_torus_tri_mesh_mex(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io double[xx], c io double[xx], c io double[x], c io double[xx], c io double[xx], c io int64_t[x])';
[sx, snx, sw, rts, rps, ier] = LineQuaaadrature_mex(mex_id_, mp, np, p, radii, scales, nosc, orig, ntri, sx, snx, sw, rts, rps, ier, 1, 1, 1, 3, 3, 1, 1, 1, 3, nsrc, 3, nsrc, nsrc, 3, nsrc, 3, nsrc, 1);
end

% --------------------------------------------------------------------------
