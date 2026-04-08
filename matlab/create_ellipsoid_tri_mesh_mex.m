function [x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr] = create_ellipsoid_tri_mesh_mex(mp, np, p, nq, ratio, nquad_bdry, nvr, ntri, x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr)
% nq = 0  : Vioreanu-Rokhlin nodes  (nvr = p*(p+1)/2)
% nq > 0  : areal quadrature nodes  (nvr = 3*nq^2)
mp         = double(mp);
np         = double(np);
p          = double(p);
nq         = double(nq);
ratio      = double(ratio);
nquad_bdry = double(nquad_bdry);
nvr        = double(nvr);
ntri       = double(ntri);
nxb        = 3 * nquad_bdry;
nptr       = ntri + 1;
nx3        = 3 * nvr;
nxbd3      = 3 * nxb;
ntv        = 9;
if nargin < 8  || isempty(x),        x        = zeros(nx3,   ntri); end
if nargin < 9  || isempty(nx),       nx       = zeros(nx3,   ntri); end
if nargin < 10 || isempty(w),        w        = zeros(nvr,   ntri); end
if nargin < 11 || isempty(xbd),      xbd      = zeros(nxbd3, ntri); end
if nargin < 12 || isempty(tri2face), tri2face = zeros(ntri,  1);    end
if nargin < 13 || isempty(tri2cell), tri2cell = zeros(2,     ntri); end
if nargin < 14 || isempty(tri_vert), tri_vert = zeros(ntv,   ntri); end
if nargin < 15 || isempty(ptr),      ptr      = zeros(nptr,  1);    end
x        = reshape(x,        nx3,   ntri);
nx       = reshape(nx,       nx3,   ntri);
xbd      = reshape(xbd,      nxbd3, ntri);
tri_vert = reshape(tri_vert, ntv,   ntri);
mex_id_ = 'create_ellipsoid_tri_mesh(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[x], c io double[xx], c io double[xx], c io double[x])';
[x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr] = LineQuaaadrature_mex(mex_id_, mp, np, p, nq, ratio, nquad_bdry, nvr, ntri, x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr, 1, 1, 1, 1, 1, 1, 1, 1, nx3, ntri, nx3, ntri, nvr, ntri, nxbd3, ntri, ntri, 2, ntri, ntv, ntri, nptr);
x        = reshape(x,        3, nvr, ntri);
nx       = reshape(nx,       3, nvr, ntri);
xbd      = reshape(xbd,      3, 3*nquad_bdry, ntri);
tri_vert = reshape(tri_vert, 3, 3, ntri);
end

% --------------------------------------------------------------------------
