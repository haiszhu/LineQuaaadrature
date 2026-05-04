function [len, lenl, lenr] = lqk_estimate_nearroot_lengths_mex(t_root, root_imag_abs, max_len_each_side, len, lenl, lenr)
t_root = double(t_root);
root_imag_abs = double(root_imag_abs);
max_len_each_side = double(max_len_each_side);
if nargin < 4
  len = 0;
end
if nargin < 5
  lenl = 0;
end
if nargin < 6
  lenr = 0;
end
len = double(len);
lenl = double(lenl);
lenr = double(lenr);
mex_id_ = 'lqk_estimate_nearroot_lengths_mex(c i double[x], c i double[x], c i int64_t[x], c io int64_t[x], c io int64_t[x], c io int64_t[x])';
[len, lenl, lenr] = LineQuaaadrature_mex(mex_id_, t_root, root_imag_abs, max_len_each_side, len, lenl, lenr, 1, 1, 1, 1, 1, 1);
end

% --------------------------------------------------------------------------
