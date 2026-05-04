function [IalphaAsvestas] = evaluate_solid_angle_integral_mex(m, tx, n, sx, snx, sw, r_vert, nbd, sxbd_in, use_nearroot, IalphaAsvestas)
m   = double(m);
n   = double(n);
nbd = double(nbd);
if nargin < 10 || isempty(use_nearroot)
  use_nearroot = 0;
end
use_nearroot = double(logical(use_nearroot));
if nargin < 11 || isempty(IalphaAsvestas)
  IalphaAsvestas = zeros(m, 1);
end
mex_id_ = 'evaluate_solid_angle_integral(c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i double[xx], c i int64_t[x], c i double[xx], c i int64_t[x], c io double[x])';
[IalphaAsvestas] = LineQuaaadrature_mex(mex_id_, m, tx, n, sx, snx, sw, r_vert, nbd, sxbd_in, use_nearroot, IalphaAsvestas, 1, 3, m, 1, 3, n, 3, n, n, 3, 3, 1, 3, nbd, 1, m);
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
