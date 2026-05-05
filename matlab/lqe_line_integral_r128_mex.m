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

