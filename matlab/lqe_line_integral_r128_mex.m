function [q_lq64] = lqe_line_integral_r128_mex(m, r0, nbd, sbdnp, nquad, sxbd, sxpbd, stangbd, sspbd, fptr_int, kdata, q_lq64)
m        = double(m);
nbd      = double(nbd);
sbdnp    = double(sbdnp);
nquad    = double(nquad);
fptr_int = double(fptr_int);
if nargin < 12 || isempty(q_lq64)
  q_lq64 = zeros(m, 1);
else
  q_lq64 = double(q_lq64(:));
end
mex_id_ = 'lqe_line_integral_r128_mex(c i int64_t[x], c i double[xx], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[xx], c i double[x], c i int64_t[x], c i double[xx], c io double[x])';
[q_lq64] = LineQuaaadrature_mex(mex_id_, m, r0, nbd, sbdnp, nquad, sxbd, sxpbd, stangbd, sspbd, fptr_int, kdata, q_lq64, 1, 3, m, 1, 1, 1, 3, nbd, 3, nbd, 3, nbd, nbd, 1, 3, m, m);
end

% --------------------------------------------------------------------------
% line_quad_compress_nearroot_mex: nearroot GL weight compression.
%
% [funvals, sxbdw] = line_quad_compress_nearroot_mex(
%     m, r0, nbd, sbdnp, nquad,
%     sxbd, sxpbd, stangbd, sspbd,
%     tgl, wgl, Dgl, w_bclag, Legmat, bclagmatlr,
%     fptr_int, kdata, funvals, sxbdw,
%     root_re, root_im, root_ok)
%
% root_re(m)  : real part of nearest singularity in parameter [-1,1]
% root_im(m)  : imaginary part (distance from real axis)
% root_ok(m)  : nonzero = use nearroot path; 0 = fall back to adaptive
% funvals, sxbdw: (nquad*sbdnp) x m here; Fortran sees (nquad,sbdnp,m).
% --------------------------------------------------------------------------
