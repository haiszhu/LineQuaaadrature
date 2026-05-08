function [funvals_pre] = lqs_eval_moments_funvals_mex(m, tx, nbd, sxbd, nquad, order, ncol, funvals_pre)
m = double(m);
nbd = double(nbd);
nquad = double(nquad);
order = double(order);
ncol = double(ncol);
nbdncol = nbd*ncol;
tx = double(reshape(tx, 3, m));
sxbd = double(reshape(sxbd, 3, nbd));
if nargin < 8 || isempty(funvals_pre)
  funvals_pre = zeros(nbdncol, m);
else
  funvals_pre = double(reshape(funvals_pre, [nbdncol, m]));
end
mex_id_ = 'lqs_eval_moments_funvals_mex(c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io double[xx])';
[funvals_pre] = LineQuaaadrature_mex(mex_id_, m, tx, nbd, sxbd, nquad, order, ncol, funvals_pre, 1, 3, m, 1, 3, nbd, 1, 1, 1, nbdncol, m);
funvals_pre = double(reshape(funvals_pre, [nbd, ncol, m]));
end

