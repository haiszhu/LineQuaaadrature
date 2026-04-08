function [tgl, wgl, Dgl] = gauss_mex(n, tgl, wgl, Dgl)
n = double(n);
tgl = zeros(n, 1);
wgl = zeros(n, 1);
Dgl = zeros(n, n);
mex_id_ = 'gauss_r64(c i int64_t[x], c io double[x], c io double[x], c io double[xx])';
[tgl, wgl, Dgl] = LineQuaaadrature_mex(mex_id_, n, tgl, wgl, Dgl, 1, n, n, n, n);
end

% --------------------------------------------------------------------------
% kernel_lookup_mex: look up a kernel subroutine by name via dlsym.
%
% Usage:
%   fptr = kernel_lookup_mex('asvestas_kernel_r64');
%
% Returns an integer(8) handle (0 = not found).
% Pass to line_kernel_eval_mex / line_quad_compress_mex as fptr_int.
% --------------------------------------------------------------------------
