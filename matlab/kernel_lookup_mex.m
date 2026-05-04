function [fptr_int] = kernel_lookup_mex(sym_name, fptr_int)
if nargin < 2
  fptr_int = 0;
end
% mwrap marshals array inputs via mxDOUBLE_CLASS; pass ASCII bytes as double.
sym_bytes  = uint8(sym_name);
sym_name_c = double([sym_bytes(:).' 0]);     % null-terminated
sym_len    = double(numel(sym_name_c));
fptr_int   = double(fptr_int);
mex_id_ = 'kernel_lookup(c i char[x], c i int64_t[x], c io int64_t[x])';
[fptr_int] = LineQuaaadrature_mex(mex_id_, sym_name_c, sym_len, fptr_int, sym_len, 1, 1);
end

% --------------------------------------------------------------------------
% lqa_root_initial_guess_mex
%
% [tinit_re, tinit_im] = lqa_root_initial_guess_mex(
%     tgl, x, y, z, n, tx, ty, tz)
% --------------------------------------------------------------------------
