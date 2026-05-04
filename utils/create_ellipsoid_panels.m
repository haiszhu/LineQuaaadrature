function [pan N] = create_ellipsoid_panels(so,o)
% create panels for part of ellipsoid according to id
%
% Hai 04/06/20

if nargin<2, o=[]; end
if isfield(o,'p'), p = o.p; else p = 8; end  % # nodes per side of panel (order)

if isfield(so,'mp'), mp = so.mp; else mp = 2; end  % # pans on small circle
if isfield(so,'np'), np = so.np; else np = 2; end  % # pans on big circle
if isfield(so,'ratio'), ratio = so.ratio; else ratio = 1; end
if isfield(so,'tau'), func = @(x) so.tau(x(1,:),x(2,:),x(3,:)); else func = @(x) ellipsoidcurvature(x(1,:),x(2,:),x(3,:)); end
tpansiz1 = 1/np; tpansiz2 = 1/mp;     % panel half-sizes in parameter
npan = mp*np; N = npan*p^2;
panind = reshape(1:npan,[mp np]);
pan = cell(6*npan,1);
for ell=1:6
  pan_ell = cell(npan,1);
  for j=1:np, for i=1:mp, k=panind(i,j);     % loop over panels (in any order)
    pan_ell{k}.psiz = [tpansiz1,tpansiz2];     % param sizes
    pan_ell{k}.inds = p^2*(k-1)+(1:p^2);     % set up indices of unknowns
    pan_ell{k}.chart = @(t) ellipsoidparam(t(1,:),t(2,:),ratio); 
    pan_ell{k}.ltrichart = @(uvs) ellipsoidparam( 2*(uvs(1,:)-1/2), 2*(uvs(2,:)-1/2),ratio); 
    pan_ell{k}.utrichart = @(uvs) ellipsoidparam(-2*(uvs(1,:)-1/2),-2*(uvs(2,:)-1/2),ratio); 
    pan_ell{k}.hcurvature = @(x) func(x); 
    pan_ell{k}.spfac = tpansiz1*tpansiz2;  % factor to convert speed to chart z
    if mp>2, if i>1 && i<mp, idxi = [i-1 i i+1]; elseif i==1, idxi = [i i+1]; else,  idxi = [i-1 i]; end, else, idxi=1:mp; end
    if np>2, if j>1 && j<np, idxj = [j-1 j j+1]; elseif j==1, idxj = [j j+1]; else,  idxj = [j-1 j]; end, else, idxj=1:np; end
    nei = panind(idxi,idxj);
    nei = nei(nei~=k);    % remove self from list, to leave 8
    pan_ell{k}.nei = nei(:);
    pan_ell{k}.mp=mp; pan_ell{k}.np=np; pan_ell{k}.topo = 'ellipsoid'; pan_ell{k}.p = p;
    pan_ell{k}.treal=[];
  end,end
pan_ell = panel_smooth_quad(pan_ell,p,ell);  % does all panels
pan(npan*(ell-1)+(1:npan)) = pan_ell;
end

%%%

function pan = panel_smooth_quad(pan,p,id)   % setup G-L tensor native quadr
[x w D] = gauss(p);   % Gauss-Legendre nodes & weights on [-1,1]
% [x w] = cheby(p);
[x10 x20] = meshgrid(x); xx = [x10(:)';x20(:)'];  % 2*p^2 parameter col vecs in R2
ww0 = w(:)*w(:)'; ww0 = ww0(:)';                     % 1*p^2 weights
One = ones(numel(x10(:)),1);

npan = pan{1}.mp*pan{1}.np; panind = reshape(1:npan,[pan{1}.mp pan{1}.np]);
for j=1:pan{1}.np, for i=1:pan{1}.mp, k=panind(i,j); 
  x1(:)=(x10(:)-(pan{1}.np-1))*pan{k}.psiz(1)+2*(j-1)/pan{1}.np;
  x2(:)=(x20(:)-(pan{1}.mp-1))*pan{k}.psiz(1)+2*(i-1)/pan{1}.mp;
  pan{k}.t = [x1(:)';x2(:)'];  
  
  x = [One,x1(:),x2(:)];
  x = x./repmat(sqrt(One+x1(:).^2+x2(:).^2),1,3);
  [pan{k}.th,pan{k}.phi] = cart2sph(x(:,1),x(:,2),x(:,3));
  pan{k}.ww = cos(pan{k}.th')./(One'+x1.^2+x2.^2);
end,end

% Unified face map: each id is represented by a latitude shift and a
% rotation around x-axis. This is algebraically equivalent to the previous
% special-case branch for id==1/4.
th_shift = [0, pi/2, pi/2, pi, pi/2, pi/2];
rot_x    = [0, 0, -pi/2, 0, -pi, -3*pi/2];

for m=1:numel(pan)
  [pan{m}.x, pan{m}.nx, pan{m}.sp] = pan{m}.chart([pan{m}.th'+th_shift(id); pan{m}.phi']);
  Rx = [1 0 0; 0 cos(rot_x(id)) sin(rot_x(id)); 0 -sin(rot_x(id)) cos(rot_x(id))];
  pan{m}.x  = Rx * pan{m}.x;
  pan{m}.nx = Rx * pan{m}.nx;
  pan{m}.hx = pan{m}.hcurvature(pan{m}.x);
  pan{m}.w = pan{m}.sp.*pan{m}.ww.*ww0*pan{m}.spfac;
  pan{m}.N = p*p;
  %

end

function [hx] = ellipsoidcurvature(x,y,z)
% this needs to be changed, just try to see if axisymmetric density can be
% assigned this way
%
% Hai 04/08/20

fH   = @(x,yz) ((x-1/2).^2).*sin(yz);
hx = fH(x,sqrt(y.^2+z.^2));


function [ x, nx, sp, dp, dt] = ellipsoidparam(p,t,ratio)
% be aware that we switched p and t in input... need to figure out how
% 
% Hai 04/05/20.
if nargin<3 , ratio=1; end


f = @(t,p) 1./sqrt( cos(p).^2.*cos(t).^2 + ratio^2*sin(p).^2.*cos(t).^2 + ratio^2*sin(t).^2);
ft = @(t,p) -1/2.*f(t,p).^3.*(-2*cos(p).^2.*cos(t).*sin(t)-2*ratio^2*sin(p).^2.*cos(t).*sin(t)+2*ratio^2*sin(t).*cos(t));
fp = @(t,p) -1/2.*f(t,p).^3.*(-2*cos(p).*sin(p).*cos(t).^2+2*ratio^2*sin(p).*cos(p).*cos(t).^2);
  
r  = @(t,p) [ f(t,p).*cos(t).*cos(p); f(t,p).*cos(t).*sin(p); f(t,p).*sin(t)];
rt = @(t,p) [ ft(t,p).*cos(t).*cos(p)-f(t,p).*sin(t).*cos(p); ft(t,p).*cos(t).*sin(p)-f(t,p).*sin(t).*sin(p); ft(t,p).*sin(t)+f(t,p).*cos(t)];
rp = @(t,p) [ fp(t,p).*cos(t).*cos(p)-f(t,p).*cos(t).*sin(p); fp(t,p).*cos(t).*sin(p)+f(t,p).*cos(t).*cos(p); fp(t,p).*sin(t)];

x  = r(t,p);
rts = rt(t,p); rps = rp(t,p);
dt = sqrt(sum(rts.^2,1));
dp = sqrt(sum(rps.^2,1));
nx =  cross(rps, rts);       % outward normal
sp = sqrt(sum(nx.^2,1)); % speeds
nx = nx.*repmat(1./sp, [3 1]);
