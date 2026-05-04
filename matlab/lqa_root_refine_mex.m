function [troot, ifconv] = lqa_root_refine_mex(xhat, yhat, zhat, n_expa, tx, ty, tz, tinit, troot, ifconv)
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

mex_id_ = 'lqa_root_refine_mex(c i double[x], c i double[x], c i double[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i dcomplex[x], c io dcomplex[x], c io int64_t[x])';
[troot, ifconv] = LineQuaaadrature_mex(mex_id_, xhat, yhat, zhat, n_expa, tx, ty, tz, tinit, troot, ifconv, n_expa, n_expa, n_expa, 1, 1, 1, 1, 1, 1, 1);
end

% --------------------------------------------------------------------------
% lqa_root_refine_mex
%
% [troot_re, troot_im, ifconv] = lqa_root_refine_mex(
%     n_expa, xhat, yhat, zhat, tx, ty, tz, tinit_re, tinit_im)
% --------------------------------------------------------------------------
