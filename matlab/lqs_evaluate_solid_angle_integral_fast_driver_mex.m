function [IalphaAsvestas] = lqs_evaluate_solid_angle_integral_fast_driver_mex(m, tx, n, sx, snx, sw, r_vert, nbd, sxbd_in, IalphaAsvestas)
m   = double(m);
n   = double(n);
nbd = double(nbd);
if nargin < 10 || isempty(IalphaAsvestas)
  IalphaAsvestas = zeros(m, 1);
end
mex_id_ = 'lqs_evaluate_solid_angle_integral_fast_driver_mex(c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i double[xx], c i int64_t[x], c i double[xx], c io double[x])';
[IalphaAsvestas] = LineQuaaadrature_mex(mex_id_, m, tx, n, sx, snx, sw, r_vert, nbd, sxbd_in, IalphaAsvestas, 1, 3, m, 1, 3, n, 3, n, n, 3, 3, 1, 3, nbd, m);
end
