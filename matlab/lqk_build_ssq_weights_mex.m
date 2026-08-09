function [w1, w3, w5, troot, xroot, yroot, zroot, rfc, rfc_ssq] = lqk_build_ssq_weights_mex(m, r0, rho_ssq, nbd, sxbd, sbdnp, nquad, tgl, wgl, Legmat, w1, w3, w5, troot, xroot, yroot, zroot, rfc, rfc_ssq)
m = double(m);
nbd = double(nbd);
sbdnp = double(sbdnp);
nquad = double(nquad);
rho_ssq = double(rho_ssq);
r0 = double(reshape(r0, 3, m));
sxbd = double(reshape(sxbd, 3, nbd));
tgl = double(tgl(:));
wgl = double(wgl(:));
Legmat = double(reshape(Legmat, nquad, nquad));
nsm = sbdnp*m;
if nargin < 11 || isempty(w1), w1 = zeros(nquad, sbdnp, m); end
if nargin < 12 || isempty(w3), w3 = zeros(nquad, sbdnp, m); end
if nargin < 13 || isempty(w5), w5 = zeros(nquad, sbdnp, m); end
if nargin < 14 || isempty(troot), troot = complex(zeros(m, sbdnp)); end
if nargin < 15 || isempty(xroot), xroot = complex(zeros(m, sbdnp)); end
if nargin < 16 || isempty(yroot), yroot = complex(zeros(m, sbdnp)); end
if nargin < 17 || isempty(zroot), zroot = complex(zeros(m, sbdnp)); end
if nargin < 18 || isempty(rfc), rfc = zeros(m, sbdnp); end
if nargin < 19 || isempty(rfc_ssq), rfc_ssq = zeros(m, 1); end
w1 = reshape(w1, nquad, nsm);
w3 = reshape(w3, nquad, nsm);
w5 = reshape(w5, nquad, nsm);
mex_id_ = 'lqk_build_ssq_weights_mex(c i int64_t[x], c i double[xx], c i double[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[xx], c io double[xx], c io double[xx], c io double[xx], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx], c io int64_t[xx], c io int64_t[x])';
[w1, w3, w5, troot, xroot, yroot, zroot, rfc, rfc_ssq] = LineQuaaadrature_mex(mex_id_, m, r0, rho_ssq, nbd, sxbd, sbdnp, nquad, tgl, wgl, Legmat, w1, w3, w5, troot, xroot, yroot, zroot, rfc, rfc_ssq, 1, 3, m, 1, 1, 3, nbd, 1, 1, nquad, nquad, nquad, nquad, nquad, nsm, nquad, nsm, nquad, nsm, m, sbdnp, m, sbdnp, m, sbdnp, m, sbdnp, m, sbdnp, m);
w1 = reshape(w1, nquad, sbdnp, m);
w3 = reshape(w3, nquad, sbdnp, m);
w5 = reshape(w5, nquad, sbdnp, m);
end

% LineQuaaadrature.mw
% mwrap interface for linequaaadrature_mod + solidangle_mod
%
% Build flags: mwrap -c99complex -i8 -mex
% After generation: perl -pi -e 's/_{2,}/_/g' LineQuaaadrature_mex.c
%
% ============================================================
% mwrap syntax:
%
%   @function [inout1, inout2] = SubroutineName_mex(in1, in2, inout1, inout2)
%   % optional MATLAB pre-processing (reshape etc.)
%   # FORTRAN SubroutineName(type[dims] in1, type[dims] in2, inout type[dims] inout1, inout type[dims] inout2);
%   end
%
% Convention: inout variables appear on BOTH sides of the MATLAB call.
% This keeps Fortran and MATLAB calling sequences exactly the same.
% type keywords: int, double, dcomplex
% dims: [1] scalar, [n] vector, [3,n] array (col-major)
%
% ============================================================
% New entries go below.
% --------------------------------------------------------------------------
