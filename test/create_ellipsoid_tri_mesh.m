function [x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr] = ...
    create_ellipsoid_tri_mesh(mp, np, p, ratio, nquad_bdry, nvr, ntri, ...
                              x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr)
% create_ellipsoid_tri_mesh
% Direct ellipsoid triangular mesh generator (array-first, no structs).
%
% Input:
%   mp, np        : panel counts on two chart directions (per face)
%   p             : Vioreanu order (nodes per triangle: p*(p+1)/2)
%   ratio         : ellipsoid ratio
%   nquad_bdry    : Gauss nodes per edge for boundary sampling
%   x,nx,w,xbd,tri2face,tri2cell,tri_vert,ptr
%                 : preallocated output arrays (in/out style)
%
% Output shapes (recommended):
%   x(3,nvr,ntri), nx(3,nvr,ntri), w(nvr,ntri)
%   xbd(3,3*nquad_bdry,ntri)
%   tri2face(ntri) or tri2face(1,ntri)
%   tri2cell(2,ntri)
%   tri_vert(3,3,ntri)
%   ptr:
%     - vector(ntri+1): flattened-node CSR pointer
%     - or ptr(2,ntri)/(4,ntri): [node_start;node_end;bdry_start;bdry_end]
%
% Hai/Codex 2026-04-07

[uvs, wts] = get_vioreanu_nodes(p-1);
nvr = size(uvs, 2);
ntri = 12 * mp * np;
nbdry = 3 * nquad_bdry;

[xg, ~, ~] = gauss(nquad_bdry);

tpansiz1 = 1/np;
tpansiz2 = 1/mp;
rotation = [0 0 -pi/2 0 -pi -3*pi/2];

% Allocate if empty, otherwise validate dimensions.
if isempty(x),  x  = zeros(3, nvr, ntri); end
if isempty(nx), nx = zeros(3, nvr, ntri); end
if isempty(w),  w  = zeros(nvr, ntri); end
if isempty(xbd), xbd = zeros(3, nbdry, ntri); end
if isempty(tri2face), tri2face = zeros(ntri, 1); end
if isempty(tri2cell), tri2cell = zeros(2, ntri); end
if isempty(tri_vert), tri_vert = zeros(3, 3, ntri); end
if isempty(ptr), ptr = zeros(ntri+1, 1); end

assert(isequal(size(x), [3, nvr, ntri]), 'x must be size [3, nvr, ntri].');
assert(isequal(size(nx), [3, nvr, ntri]), 'nx must be size [3, nvr, ntri].');
assert(isequal(size(w), [nvr, ntri]), 'w must be size [nvr, ntri].');
assert(isequal(size(xbd), [3, nbdry, ntri]), 'xbd must be size [3, 3*nquad_bdry, ntri].');
assert(size(tri2cell,1) == 2 && size(tri2cell,2) == ntri, 'tri2cell must be size [2, ntri].');
assert(isequal(size(tri_vert), [3, 3, ntri]), 'tri_vert must be size [3, 3, ntri].');
assert(numel(tri2face) == ntri, 'tri2face must contain ntri entries.');

itri = 0;
for face = 1:6
  for j = 1:np
    for i = 1:mp
      x1lo = ((-1) - (np-1))*tpansiz1 + 2*(j-1)/np;
      x1hi = ((+1) - (np-1))*tpansiz1 + 2*(j-1)/np;
      x2lo = ((-1) - (mp-1))*tpansiz2 + 2*(i-1)/mp;
      x2hi = ((+1) - (mp-1))*tpansiz2 + 2*(i-1)/mp;

      vll = [x1lo; x2lo];
      vur = [x1hi; x2hi];
      vul = [x1lo; x2hi];
      vlr = [x1hi; x2lo];

      rll = face_map_point(vll(1), vll(2), face, ratio, rotation);
      rur = face_map_point(vur(1), vur(2), face, ratio, rotation);
      rul = face_map_point(vul(1), vul(2), face, ratio, rotation);
      rlr = face_map_point(vlr(1), vlr(2), face, ratio, rotation);

      d1 = sum((rll-rur).^2);
      d2 = sum((rul-rlr).^2);
      if d1 + 1e-13 > d2
        tri12 = { [vll, vlr, vul], [vur, vul, vlr] };
      else
        tri12 = { [vul, vll, vur], [vlr, vur, vll] };
      end

      for kk = 1:2
        itri = itri + 1;
        V = tri12{kk};

        x12 = V(:,1) + (V(:,2)-V(:,1))*uvs(1,:) + (V(:,3)-V(:,1))*uvs(2,:);
        [xt, nxt, spt, tht] = face_map_nodes(x12(1,:), x12(2,:), face, ratio, rotation);
        wwt = cos(tht) ./ (1 + x12(1,:).^2 + x12(2,:).^2);
        jac = abs(det([V(:,2)-V(:,1), V(:,3)-V(:,1)]));
        wt = (spt .* wwt .* wts(:)') * jac;

        e1 = edge_nodes(V(:,1), V(:,2), xg);
        e2 = edge_nodes(V(:,2), V(:,3), xg);
        e3 = edge_nodes(V(:,3), V(:,1), xg);
        xb1 = face_map_nodes(e1(1,:), e1(2,:), face, ratio, rotation);
        xb2 = face_map_nodes(e2(1,:), e2(2,:), face, ratio, rotation);
        xb3 = face_map_nodes(e3(1,:), e3(2,:), face, ratio, rotation);

        x(:,:,itri) = xt;
        nx(:,:,itri) = nxt;
        w(:,itri) = wt(:);
        xbd(:,:,itri) = [xb1, xb2, xb3];
        tri_vert(:,:,itri) = [face_map_point(V(1,1), V(2,1), face, ratio, rotation), ...
                              face_map_point(V(1,2), V(2,2), face, ratio, rotation), ...
                              face_map_point(V(1,3), V(2,3), face, ratio, rotation)];
        tri2face(itri) = face;
        tri2cell(:,itri) = [i; j];
      end
    end
  end
end

% ptr bookkeeping for flattened storage.
if isvector(ptr) && numel(ptr) == ntri + 1
  for k = 1:ntri+1
    ptr(k) = 1 + (k-1) * nvr;
  end
elseif size(ptr,2) == ntri && size(ptr,1) >= 2
  for k = 1:ntri
    ptr(1,k) = 1 + (k-1) * nvr;
    ptr(2,k) = k * nvr;
    if size(ptr,1) >= 4
      ptr(3,k) = 1 + (k-1) * nbdry;
      ptr(4,k) = k * nbdry;
    end
  end
end

end

function e = edge_nodes(a, b, s)
e = 0.5*(1-s(:)') .* a + 0.5*(1+s(:)') .* b;
end

function r = face_map_point(x1, x2, face, ratio, rotation)
[x, ~, ~] = face_map_nodes(x1, x2, face, ratio, rotation);
r = x(:,1);
end

function [x, nx, sp, th] = face_map_nodes(x1, x2, face, ratio, rotation)
one = ones(size(x1));
q = [one; x1; x2] ./ sqrt(one + x1.^2 + x2.^2);
[th, phi] = cart2sph(q(1,:), q(2,:), q(3,:));

if (face == 1) || (face == 4)
  [x, nx, sp] = ellipsoidparam(th + floor(face/3)*pi, phi, ratio);
else
  [x, nx, sp] = ellipsoidparam(th + pi/2, phi, ratio);
  ang = rotation(face);
  Rx = [1 0 0; 0 cos(ang) sin(ang); 0 -sin(ang) cos(ang)];
  x = Rx * x;
  nx = Rx * nx;
end
end

function [x, nx, sp, dp, dt] = ellipsoidparam(p, t, ratio)
if nargin < 3, ratio = 1; end

f  = @(tt,pp) 1./sqrt(cos(pp).^2.*cos(tt).^2 + ratio^2*sin(pp).^2.*cos(tt).^2 + ratio^2*sin(tt).^2);
ft = @(tt,pp) -0.5.*f(tt,pp).^3.*(-2*cos(pp).^2.*cos(tt).*sin(tt) -2*ratio^2*sin(pp).^2.*cos(tt).*sin(tt) +2*ratio^2*sin(tt).*cos(tt));
fp = @(tt,pp) -0.5.*f(tt,pp).^3.*(-2*cos(pp).*sin(pp).*cos(tt).^2 +2*ratio^2*sin(pp).*cos(pp).*cos(tt).^2);

r  = @(tt,pp) [f(tt,pp).*cos(tt).*cos(pp); f(tt,pp).*cos(tt).*sin(pp); f(tt,pp).*sin(tt)];
rt = @(tt,pp) [ft(tt,pp).*cos(tt).*cos(pp)-f(tt,pp).*sin(tt).*cos(pp); ...
               ft(tt,pp).*cos(tt).*sin(pp)-f(tt,pp).*sin(tt).*sin(pp); ...
               ft(tt,pp).*sin(tt)+f(tt,pp).*cos(tt)];
rp = @(tt,pp) [fp(tt,pp).*cos(tt).*cos(pp)-f(tt,pp).*cos(tt).*sin(pp); ...
               fp(tt,pp).*cos(tt).*sin(pp)+f(tt,pp).*cos(tt).*cos(pp); ...
               fp(tt,pp).*sin(tt)];

x = r(t,p);
rts = rt(t,p); rps = rp(t,p);
dt = sqrt(sum(rts.^2,1));
dp = sqrt(sum(rps.^2,1));
nx = cross(rps, rts);
sp = sqrt(sum(nx.^2,1));
nx = nx ./ sp;
end
