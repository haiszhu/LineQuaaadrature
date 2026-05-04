function [funvals] = lqk_eval_mex(m, r0, nbd, sbdnp, nquad, sxbd, sxpbd, stangbd, fptr_int, kdata, funvals)
m        = double(m);
nbd      = double(nbd);
sbdnp    = double(sbdnp);
nquad    = double(nquad);
fptr_int = double(fptr_int);
nqs     = nquad * sbdnp;
funvals = zeros(nqs, m);
mex_id_ = 'lqk_eval_mex(c i int64_t[x], c i double[xx], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[xx], c i int64_t[x], c i double[xx], c io double[xx])';
[funvals] = LineQuaaadrature_mex(mex_id_, m, r0, nbd, sbdnp, nquad, sxbd, sxpbd, stangbd, fptr_int, kdata, funvals, 1, 3, m, 1, 1, 1, 3, nbd, 3, nbd, 3, nbd, 1, 3, m, nqs, m);
end

