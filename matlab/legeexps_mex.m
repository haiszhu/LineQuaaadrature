function [x, u, v, whts] = legeexps_mex(itype, n, x, u, v, whts)
itype = double(itype);
n = double(n);
x = double(x(:));
u = double(reshape(u, n, n));
v = double(reshape(v, n, n));
whts = double(whts(:));
mex_id_ = 'legeexps_mex(c i int64_t[x], c i int64_t[x], c io double[x], c io double[xx], c io double[xx], c io double[x])';
[x, u, v, whts] = LineQuaaadrature_mex(mex_id_, itype, n, x, u, v, whts, 1, 1, n, n, n, n, n, n);
end

% --------------------------------------------------------------------------
