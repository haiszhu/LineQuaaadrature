function [troot, ifconv] = lqa_root_refine_r128_mex(xhat, yhat, zhat, n_expa, tx, ty, tz, tinit, troot, ifconv)
n_expa   = double(n_expa);
xhat     = double(xhat(:));
yhat     = double(yhat(:));
zhat     = double(zhat(:));
tx       = double(tx);
ty       = double(ty);
tz       = double(tz);

if nargin < 9 || isempty(troot)
  troot = 0;
end
if nargin < 10 || isempty(ifconv)
  ifconv = 0;
end

mex_id_ = 'lqa_root_refine_r128_mex(c i double[x], c i double[x], c i double[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i dcomplex[x], c io dcomplex[x], c io int64_t[x])';
[troot, ifconv] = LineQuaaadrature_mex(mex_id_, xhat, yhat, zhat, n_expa, tx, ty, tz, tinit, troot, ifconv, n_expa, n_expa, n_expa, 1, 1, 1, 1, 1, 1, 1);
end

% --------------------------------------------------------------------------
% lqk_eval_mex: evaluate kernel at all boundary quad nodes.
%
% [funvals] = lqk_eval_mex(m, r0, nbd, sbdnp, nquad,
%                                   sxbd, sxpbd, stangbd,
%                                   fptr_int, kdim, kdata, funvals)
%
% fptr_int: integer(8) from kernel_lookup_mex
% kdata:    kdim x m  (per-target kernel data, e.g. kdata(:,j) = qhat)
% funvals:  nquad x sbdnp x m  (output)
% --------------------------------------------------------------------------
