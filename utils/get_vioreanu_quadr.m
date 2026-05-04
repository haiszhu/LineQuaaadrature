function [sl,sr] = get_vioreanu_quadr(sx,order,norder)
% get quadrature 
% from quadrilateral patch to two triangular patches
% assume GL nodes as input
%
% add three vertices, with order [0 1 0; 0 0 1] in parametric space
% (r_vert = sbd.trichart([0 1 0; 0 0 1]) in various places )
%
% Hai 01/25/23, add one level of h-refinement

% 'sx' is assumed to be GL node of 'order'
% norder is vioreanu rokhlin node order
% norder4sub is the order for 4 sub triangles

[uvs,wts]=get_vioreanu_nodes(norder-1); 

% interpolate everything from tensor grid
[x0,w0,D]=gauss(order);

% compute derivative numerically
tp_ind = reshape(1:order^2,[order order]);
xts = nan(1,order^2); yts = xts; zts = xts; % partial t (small circle)
for j=1:order
  xk = sx(1,tp_ind(:,j)); yk = sx(2,tp_ind(:,j)); zk = sx(3,tp_ind(:,j));
  xts(tp_ind(:,j)) = D*xk(:); yts(tp_ind(:,j)) = D*yk(:); zts(tp_ind(:,j)) = D*zk(:);
end
rts = [xts(:),yts(:),zts(:)]';
xps = nan(1,order^2); yps = xps; zps = xps; % partial p (big circle)
for j=1:order
  x = sx(1,tp_ind(j,:)); y = sx(2,tp_ind(j,:)); z = sx(3,tp_ind(j,:));
  xps(tp_ind(j,:)) = D*x(:); yps(tp_ind(j,:)) = D*y(:); zps(tp_ind(j,:)) = D*z(:);
end
rps = [xps(:),yps(:),zps(:)]';

% compute 2nd order derivative numerically
xtts = nan(1,order^2); ytts = xtts; ztts = xtts;
for j=1:order
  xtts(tp_ind(:,j)) = D*xts(tp_ind(:,j))'; 
  ytts(tp_ind(:,j)) = D*yts(tp_ind(:,j))'; 
  ztts(tp_ind(:,j)) = D*zts(tp_ind(:,j))';
end
rtts = [xtts(:),ytts(:),ztts(:)]';
xpps = nan(1,order^2); ypps = xpps; zpps = xpps;
for j=1:order
  xpps(tp_ind(j,:)) = D*xps(tp_ind(j,:))'; 
  ypps(tp_ind(j,:)) = D*yps(tp_ind(j,:))'; 
  zpps(tp_ind(j,:)) = D*zps(tp_ind(j,:))';
end
rpps = [xpps(:),ypps(:),zpps(:)]';
xtps = nan(1,order^2); ytps = xtps; ztps = xtps;
for j=1:order
  xtps(tp_ind(j,:)) = D*xts(tp_ind(j,:))'; 
  ytps(tp_ind(j,:)) = D*yts(tp_ind(j,:))'; 
  ztps(tp_ind(j,:)) = D*zts(tp_ind(j,:))';
end
rtps = [xtps(:),ytps(:),ztps(:)]';

% 4 corners, and cut rectangular into two triangles (make a decision for skewed parallelogram)
rll = sx(:,1); rur = sx(:,end); % lower left & upper right
rul = sx(:,order); rlr = sx(:,end-order+1); % upper left & lower right
d1 = sum((rll-rur).^2); d2 = sum((rul - rlr).^2);
if d1+1e-13 > d2 % cut line lower left & upper right
  uvsl = 2*uvs - 1; 
  uvsr = -uvsl;
  rl_vert = [rll rlr rul]; % three vertices
  rr_vert = [rur rul rlr];
else % cut line upper left & lower right
  uvsl = 2*uvs - 1; 
  uvsr = -uvsl;
  uvslc = uvsl(1,:)+1i*uvsl(2,:);
  uvsrc = uvsr(1,:)+1i*uvsr(2,:);
  uvslc = uvslc*exp(-1i*pi/2);
  uvsrc = uvsrc*exp(-1i*pi/2);
  uvsl = [real(uvslc);imag(uvslc)];
  uvsr = [real(uvsrc);imag(uvsrc)];
  rl_vert = [rul rll rur]; % three vertices
  rr_vert = [rlr rur rll];
%   uvsl(2,:) = -uvsl(2,:); % flip y 
%   uvsr(2,:) = -uvsr(2,:); % this changes orientation, normal direction won't be consistent
end

% interpolate from tensor grid to vioreanu grid
tmpl = interpval(uvsl,[sx;rts;rps;rtts;rpps;rtps],x0);
slx = tmpl(1:3,:);
rlts = tmpl(4:6,:);  % partial t
rlps = tmpl(7:9,:);  % partial p
rltts = tmpl(10:12,:);  % partial tt
rlpps = tmpl(13:15,:);  % partial pp
rltps = tmpl(16:18,:);  % partial tp 
slnx =  cross(rlps, rlts); % outward normal
slsp = sqrt(sum(slnx.^2,1)); % speeds
slnx = slnx./slsp;  % surface normal
slw = slsp.*wts*4;  % quadrature weights
sl.x = slx; sl.nx = slnx; sl.w = slw;
sl.r_vert = rl_vert; % three vertices ~ s.xlo & s.xhi in 2d
sl.p = order;
sl.uvs = uvs; % 
sl.uvsl = uvsl;
sl.xpt = rlts; sl.xps = rlps; % partial derivative
sl.xptt = rltts; sl.xpss = rlpps; sl.xpts = rltps; % 2nd order derivative

tmpr = interpval(uvsr,[sx;rts;rps;rtts;rpps;rtps],x0);
srx = tmpr(1:3,:);
rrts = tmpr(4:6,:);  % partial t
rrps = tmpr(7:9,:);  % partial p
rrtts = tmpr(10:12,:);  % partial tt
rrpps = tmpr(13:15,:);  % partial pp
rrtps = tmpr(16:18,:);  % partial tp 
srnx =  cross(rrps, rrts); % outward normal
srsp = sqrt(sum(srnx.^2,1)); % speeds
srnx = srnx./srsp;  % surface normal
srw = srsp.*wts*4;  % quadrature weights
sr.x = srx; sr.nx = srnx; sr.w = srw;
sr.r_vert = rr_vert;
sr.p = order;
sr.uvs = uvs; % 
sr.uvsr = uvsr; 
sr.xpt = rrts; sr.xps = rrps; % partial derivative
sr.xptt = rrtts; sr.xpss = rrpps; sr.xpts = rrtps; % 2nd order derivative

% curvature
fl1E = sum(sl.xpt.*sl.xpt); fl1F = sum(sl.xpt.*sl.xps); fl1G = sum(sl.xps.*sl.xps);
fl2L = sum(sl.xptt.*sl.nx); fl2M = sum(sl.xpts.*sl.nx); fl2N = sum(sl.xpss.*sl.nx);
sl.cur = (fl1E.*fl2N - 2*fl1F.*fl2M + fl1G.*fl2L)./(2*(fl1E.*fl1G - fl1F.^2));
sl.kcur = (fl2L.*fl2N - fl2M.^2)./(fl1E.*fl1G - fl1F.^2);
fr1E = sum(sl.xpt.*sl.xpt); fr1F = sum(sl.xpt.*sl.xps); fr1G = sum(sl.xps.*sl.xps);
fr2L = sum(sl.xptt.*sl.nx); fr2M = sum(sl.xpts.*sl.nx); fr2N = sum(sl.xpss.*sl.nx);
sr.cur = (fr1E.*fr2N - 2*fr1F.*fr2M + fr1G.*fr2L)./(2*(fr1E.*fr1G - fr1F.^2));
sr.kcur = (fr2L.*fr2N - fr2M.^2)./(fr1E.*fr1G - fr1F.^2);

end

function t_val = interpval(t,s_val,x0)
% using 1d interpolation, assume GL nodes.
% interpolate 3d coordinates s_val = (x,y,z) to subpanel parametrized by t
%
% easier to debug if interpolation is accurate for tiny patches
% assume RCIP works with 1d interpolation, then this should also work?
%
% interpolate from 3xN source to 3*Nt targets
% s_val is source value, assume x otimes x grid, where x is 1d GL
% t_val is target value, assume [t(1,:);t(2,:)]...
%
% 12/22/22 Hai

% (x,y,z) of s_val
p = length(x0);

% intptmp1 & intptmp2
intptmp = interpmat_1d([t(1,:)';t(2,:)'],x0); % 1d GL nodes along t(2) direction

% 1st interpolate to (t(1),x0), where x0 is gauss legendre nodes
% intptmp1 = interpmat_1d(t(1,:)',x0); % 1d GL nodes along t(2) direction
intptmp1 = intptmp(1:end/2,:);

% then interpolate to (t(1),t(2))
% intptmp2 = interpmat_1d(t(2,:)',x0); 
intptmp2 = intptmp(end/2+1:end,:);

% more general
t_val = zeros(size(s_val,1),size(t,2));
for k=1:size(s_val,1) % assume each row of s_val needs to be interpilated to target t
  s_valk = reshape(s_val(k,:),p,p);
  s_valtmp = intptmp1*s_valk';
  t_val(k,:) = sum((s_valtmp.*intptmp2)');
end

end

function L = interpmat_1d(t,s)
% INTERPMAT_1D   interpolation matrix from nodes in 1D to target nodes
%
% L = interpmat_1d(t,s) returns interpolation matrix taking values on nodes s
%  to target nodes t. Computed in Helsing style.

% bits taken from qplp/interpmatrix.m from 2013
% Barnett 7/17/16

if nargin==0, test_interpmat_1d; return; end

p = numel(s); q = numel(t); s = s(:); t = t(:);       % all col vecs
n = p; % set the polynomial order we go up to (todo: check why bad if not p)
V = ones(p,n); for j=2:n, V(:,j) = V(:,j-1).*s; end   % polyval matrix on nodes
R = ones(q,n); for j=2:n, R(:,j) = R(:,j-1).*t; end   % polyval matrix on targs
L = (V'\R')'; % backwards-stable way to do it (Helsing) See corners/interpdemo.m
end 