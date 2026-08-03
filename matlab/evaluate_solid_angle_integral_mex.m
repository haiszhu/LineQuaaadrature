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

