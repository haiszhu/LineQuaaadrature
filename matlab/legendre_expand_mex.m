function [coeffs] = legendre_expand_mex(nn, tnodes, wnodes, values, coeffs)
nn = double(nn);
tnodes = double(tnodes(:));
wnodes = double(wnodes(:));
values = double(values(:));
coeffs = double(coeffs(:));
mex_id_ = 'legendre_expand_mex(c i int64_t[x], c i double[x], c i double[x], c i double[x], c io double[x])';
[coeffs] = LineQuaaadrature_mex(mex_id_, nn, tnodes, wnodes, values, coeffs, 1, nn, nn, nn, nn);
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
