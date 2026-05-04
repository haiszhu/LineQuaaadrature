function [tinit] = lqa_root_initial_guess_mex(tgl, x, y, z, n, tx, ty, tz, tinit)
n   = double(n);
tgl = double(tgl(:));
x   = double(x(:));
y   = double(y(:));
z   = double(z(:));
tx  = double(tx);
ty  = double(ty);
tz  = double(tz);

if nargin < 9 || isempty(tinit)
  tinit = 0;
end

mex_id_ = 'lqa_root_initial_guess(c i double[x], c i double[x], c i double[x], c i double[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c io dcomplex[x])';
[tinit] = LineQuaaadrature_mex(mex_id_, tgl, x, y, z, n, tx, ty, tz, tinit, n, n, n, n, 1, 1, 1, 1, 1);
end


% --------------------------------------------------------------------------
% lqa_root_refine_mex
%
% [troot_re, troot_im, ifconv] = lqa_root_refine_mex(
%     n_expa, xhat, yhat, zhat, tx, ty, tz, tinit_re, tinit_im)
% --------------------------------------------------------------------------
