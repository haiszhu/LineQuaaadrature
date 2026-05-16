function [funvals, sxbdw] = line_quad_compress_mex(m, r0, nbd, sbdnp, nquad, sxbd, sxpbd, stangbd, sspbd, tgl, wgl, Dgl, w_bclag, Legmat, bclagmatlr, fptr_int, kdata, funvals, sxbdw)
m        = double(m);
nbd      = double(nbd);
sbdnp    = double(sbdnp);
nquad    = double(nquad);
fptr_int = double(fptr_int);
nqs     = nquad * sbdnp;
if nargin < 18 || isempty(funvals)
  funvals = zeros(nqs, m);
else
  funvals = reshape(funvals, nqs, m);
end
if nargin < 19 || isempty(sxbdw)
  sxbdw = zeros(nqs, m);
else
  sxbdw = reshape(sxbdw, nqs, m);
end
mex_id_ = 'line_quad_compress_mex(c i int64_t[x], c i double[xx], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[xx], c i double[x], c i double[x], c i double[x], c i double[xx], c i double[x], c i double[xx], c i double[xx], c i int64_t[x], c i double[xx], c io double[xx], c io double[xx])';
[funvals, sxbdw] = LineQuaaadrature_mex(mex_id_, m, r0, nbd, sbdnp, nquad, sxbd, sxpbd, stangbd, sspbd, tgl, wgl, Dgl, w_bclag, Legmat, bclagmatlr, fptr_int, kdata, funvals, sxbdw, 1, 3, m, 1, 1, 1, 3, nbd, 3, nbd, 3, nbd, nbd, nquad, nquad, nquad, nquad, nquad, nquad, nquad, nquad, 2, 1, 3, m, nqs, m, nqs, m);
end


