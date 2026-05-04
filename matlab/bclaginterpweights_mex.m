function [w] = bclaginterpweights_mex(n, x, w)
n = double(n);
x = double(x(:));
w = double(w(:));
mex_id_ = 'bclaginterpweights_mex(c i int64_t[x], c i double[x], c io double[x])';
[w] = LineQuaaadrature_mex(mex_id_, n, x, w, 1, n, n);
end

% --------------------------------------------------------------------------
