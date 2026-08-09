function [rfc, IalphaAsvestas] = lqs_evaluate_solid_angle_integral_fast_mex(m, r0, nbd, sbdnp, nquad, sxbd, stangbd, sspbd, len1, sxbd1, stangbd1, swbd1, len2, sxbd2, stangbd2, swbd2, len3, sxbd3, stangbd3, swbd3, qhat, tgl, wgl, Dgl, w_bclag, bclagmatlr, troot, xroot, yroot, zroot, rfc, IalphaAsvestas, rho_in, sxbd_raw, tx_raw, Rfr, alpha_fr, Legmat)
m = double(m);
nbd = double(nbd);
sbdnp = double(sbdnp);
nquad = double(nquad);
len1 = double(len1);
len2 = double(len2);
len3 = double(len3);
nb1 = len1*nbd;
nb2 = len2*nbd;
nb3 = len3*nbd;
r0 = double(reshape(r0, 3, m));
sxbd = double(reshape(sxbd, 3, nbd));
stangbd = double(reshape(stangbd, 3, nbd));
sspbd = double(sspbd(:));
sxbd1 = double(reshape(sxbd1, 3, nb1));
stangbd1 = double(reshape(stangbd1, 3, nb1));
swbd1 = double(swbd1(:));
sxbd2 = double(reshape(sxbd2, 3, nb2));
stangbd2 = double(reshape(stangbd2, 3, nb2));
swbd2 = double(swbd2(:));
sxbd3 = double(reshape(sxbd3, 3, nb3));
stangbd3 = double(reshape(stangbd3, 3, nb3));
swbd3 = double(swbd3(:));
qhat = double(qhat(:));
tgl = double(tgl(:));
wgl = double(wgl(:));
Dgl = double(reshape(Dgl, nquad, nquad));
w_bclag = double(w_bclag(:));
bclagmatlr = double(reshape(bclagmatlr, nquad, 2));
troot = double(reshape(troot, m, sbdnp));
xroot = double(reshape(xroot, m, sbdnp));
yroot = double(reshape(yroot, m, sbdnp));
zroot = double(reshape(zroot, m, sbdnp));
if nargin < 31 || isempty(rfc), rfc = zeros(m, sbdnp); end
if nargin < 32 || isempty(IalphaAsvestas), IalphaAsvestas = zeros(m, 1); end
rho_in = double(rho_in);
sxbd_raw = double(reshape(sxbd_raw, 3, nbd));
tx_raw = double(reshape(tx_raw, 3, m));
Rfr = double(reshape(Rfr, 3, 3));
alpha_fr = double(alpha_fr);
Legmat = double(reshape(Legmat, nquad, nquad));
mex_id_ = 'lqs_evaluate_solid_angle_integral_fast_mex(c i int64_t[x], c i double[xx], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i double[x], c i double[x], c i double[x], c i double[xx], c i double[x], c i double[xx], c i dcomplex[xx], c i double[xx], c i double[xx], c i double[xx], c io int64_t[xx], c io double[x], c i double[x], c i double[xx], c i double[xx], c i double[xx], c i double[x], c i double[xx])';
[rfc, IalphaAsvestas] = LineQuaaadrature_mex(mex_id_, m, r0, nbd, sbdnp, nquad, sxbd, stangbd, sspbd, len1, sxbd1, stangbd1, swbd1, len2, sxbd2, stangbd2, swbd2, len3, sxbd3, stangbd3, swbd3, qhat, tgl, wgl, Dgl, w_bclag, bclagmatlr, troot, xroot, yroot, zroot, rfc, IalphaAsvestas, rho_in, sxbd_raw, tx_raw, Rfr, alpha_fr, Legmat, 1, 3, m, 1, 1, 1, 3, nbd, 3, nbd, nbd, 1, 3, nb1, 3, nb1, nb1, 1, 3, nb2, 3, nb2, nb2, 1, 3, nb3, 3, nb3, nb3, 3, nquad, nquad, nquad, nquad, nquad, nquad, 2, m, sbdnp, m, sbdnp, m, sbdnp, m, sbdnp, m, sbdnp, m, 1, 3, nbd, 3, m, 3, 3, 1, nquad, nquad);
end


% --------------------------------------------------------------------------
% line_quad_compress_mex: error-driven GL bisection weight compression.
%
% [funvals, sxbdw] = line_quad_compress_mex(
%     m, r0, nbd, sbdnp, nquad,
%     sxbd, sxpbd, stangbd, sspbd,
%     tgl, wgl, Dgl, w_bclag, Legmat, bclagmatlr,
%     fptr_int, kdata, funvals, sxbdw)
%
% funvals, sxbdw: (nquad*sbdnp) x m here; Fortran sees (nquad,sbdnp,m).
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
% create_ellipsoid_tri_mesh_mex: build VR triangular mesh on an ellipsoid.
%
% [x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr] = ...
%     create_ellipsoid_tri_mesh_mex(mp, np, p, ratio, nquad_bdry, nvr, ntri, ...)
%
% Pass [] for any output array to auto-allocate.
% 3D arrays x/nx/xbd/tri_vert are passed as 2D (first dims flattened).
% Reshape outputs back after the call:
%   x        = reshape(x,        3, nvr, ntri)
%   nx       = reshape(nx,       3, nvr, ntri)
%   xbd      = reshape(xbd,      3, 3*nquad_bdry, ntri)
%   tri_vert = reshape(tri_vert, 3, 3, ntri)
% tri2face(ntri), tri2cell(2,ntri), ptr(ntri+1) returned as double.
% --------------------------------------------------------------------------
