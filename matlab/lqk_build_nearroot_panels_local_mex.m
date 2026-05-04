function [tpan, upan, wpan, xpan, ypan, zpan, xdisp, ydisp, zdisp, stangpan, sppan, dswpan] = lqk_build_nearroot_panels_local_mex(t0, nquad, npan, lenl, lenr, tgl, wgl, xjhat, yjhat, zjhat, rbase, tpan, upan, wpan, xpan, ypan, zpan, xdisp, ydisp, zdisp, stangpan, sppan, dswpan)
t0 = double(t0);
nquad = double(nquad);
npan = double(npan);
lenl = double(lenl);
lenr = double(lenr);
nquad3 = 3*nquad;
tgl = double(tgl(:));
wgl = double(wgl(:));
xjhat = double(xjhat(:));
yjhat = double(yjhat(:));
zjhat = double(zjhat(:));
rbase = double(rbase(:));
tpan = double(reshape(tpan, nquad, npan));
upan = double(reshape(upan, nquad, npan));
wpan = double(reshape(wpan, nquad, npan));
xpan = double(reshape(xpan, nquad, npan));
ypan = double(reshape(ypan, nquad, npan));
zpan = double(reshape(zpan, nquad, npan));
xdisp = double(reshape(xdisp, nquad, npan));
ydisp = double(reshape(ydisp, nquad, npan));
zdisp = double(reshape(zdisp, nquad, npan));
stangpan = double(reshape(stangpan, nquad3, npan));
sppan = double(reshape(sppan, nquad, npan));
dswpan = double(reshape(dswpan, nquad, npan));
mex_id_ = 'lqk_build_nearroot_panels_local_mex(c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx])';
[tpan, upan, wpan, xpan, ypan, zpan, xdisp, ydisp, zdisp, stangpan, sppan, dswpan] = LineQuaaadrature_mex(mex_id_, t0, nquad, npan, lenl, lenr, tgl, wgl, xjhat, yjhat, zjhat, rbase, tpan, upan, wpan, xpan, ypan, zpan, xdisp, ydisp, zdisp, stangpan, sppan, dswpan, 1, 1, 1, 1, 1, nquad, nquad, nquad, nquad, nquad, 3, nquad, npan, nquad, npan, nquad, npan, nquad, npan, nquad, npan, nquad, npan, nquad, npan, nquad, npan, nquad, npan, nquad3, npan, nquad, npan, nquad, npan);
stangpan = double(reshape(stangpan, 3, nquad, npan));
end

% --------------------------------------------------------------------------
