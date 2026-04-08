function [funvals] = line_kernel_eval_mex(m, r0, nbd, sbdnp, nquad, sxbd, sxpbd, stangbd, fptr_int, kdata, funvals)
m        = double(m);
nbd      = double(nbd);
sbdnp    = double(sbdnp);
nquad    = double(nquad);
fptr_int = double(fptr_int);
nqs     = nquad * sbdnp;
funvals = zeros(nqs, m);
mex_id_ = 'line_kernel_eval_mex(c i int64_t[x], c i double[xx], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[xx], c i int64_t[x], c i double[xx], c io double[xx])';
[funvals] = LineQuaaadrature_mex(mex_id_, m, r0, nbd, sbdnp, nquad, sxbd, sxpbd, stangbd, fptr_int, kdata, funvals, 1, 3, m, 1, 1, 1, 3, nbd, 3, nbd, 3, nbd, 1, 3, m, nqs, m);
end

% --------------------------------------------------------------------------
% evaluate_solid_angle_integral_mex: compute Asvestas solid angle integral.
%
% [IalphaAsvestas] = evaluate_solid_angle_integral_mex(
%     m, tx, n, sx, snx, sw, r_vert, nbd, sxbd_in, IalphaAsvestas)
%
% tx(3,m)        : target points
% sx(3,n)        : source positions (VR nodes)
% snx(3,n)       : source normals (for qhat)
% sw(n)          : source weights
% r_vert(3,3)    : triangle vertices for circumcircle transform
% nbd            : number of boundary quad nodes (= 3*nquad)
% sxbd_in(3,nbd) : analytic boundary positions at GL nodes
% IalphaAsvestas(m) : output solid angle values
% --------------------------------------------------------------------------
