% hello test
%
% 04/07/26 Hai

rootdir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootdir, '../matlab'))
addpath('~/git/FMM3D/matlab/')

side = 'e';
type = 'ellipsoid'; 
nplotpts = 52; 

order = 14;
mp = 8;
np = 8;
ratio = 1.0;
nquad_bdry = order + 8;

nterms = order;
nvr  = order*(order+1)/2;
ntri = 12*mp*np;

% % preallocate mesh arrays (Fortran-style)
% x        = zeros(3, nvr, ntri);
% nx       = zeros(3, nvr, ntri);
% w        = zeros(nvr, ntri);
% xbd      = zeros(3, 3*nquad_bdry, ntri);
% tri2face = zeros(ntri,1);
% tri2cell = zeros(2,ntri);
% tri_vert = zeros(3,3,ntri);
% ptr      = zeros(ntri+1,1);

% ellipsoid triangular mesh (direct, no struct internals)
% [x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr] = ...
%   create_ellipsoid_tri_mesh(mp, np, order, ratio, nquad_bdry, nvr, ntri, ...
%                             x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr);
% [x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr] = ...                                                                               
%       create_ellipsoid_tri_mesh_mex(mp, np, order, ratio, nquad_bdry, nvr, ntri, ...                                                         
%                               x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr);  

% VR (existing):                                                                                                                       
nq = 0;  nvr = order*(order+1)/2;    
x        = zeros(3, nvr, ntri);
nx       = zeros(3, nvr, ntri);
w        = zeros(nvr, ntri);
xbd      = zeros(3, 3*nquad_bdry, ntri);
tri2face = zeros(ntri,1);
tri2cell = zeros(2,ntri);
tri_vert = zeros(3,3,ntri);
ptr      = zeros(ntri+1,1);
[x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr] = ...                                                                               
      create_ellipsoid_tri_mesh_mex(mp, np, order, nq, ratio, nquad_bdry, nvr, ntri, ...                                                         
                              x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr);                                            
                                                                                                                                       
% % Areal:                                                                                                                               
% nq = 12;  nvr = 3*nq^2;   
% x        = zeros(3, nvr, ntri);
% nx       = zeros(3, nvr, ntri);
% w        = zeros(nvr, ntri);
% xbd      = zeros(3, 3*nquad_bdry, ntri);
% tri2face = zeros(ntri,1);
% tri2cell = zeros(2,ntri);
% tri_vert = zeros(3,3,ntri);
% ptr      = zeros(ntri+1,1);
% [x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr] = ...                                                                               
%       create_ellipsoid_tri_mesh_mex(mp, np, order, nq, ratio, nquad_bdry, nvr, ntri, ...                                                         
%                               x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr);   

N = ntri * nvr;
sx = reshape(x, 3, []);
snx = reshape(nx, 3, []);
sw = reshape(w, 1, []);
s.x = sx; s.nx = snx; s.w = sw;
so.ratio = ratio;

% targets
rng(10)
domain = [-0.5 2 -0.5 2 -1.25 1.25];
[xx, yy, zz] = meshgrid(linspace(domain(1), domain(2), nplotpts), ...
                        linspace(domain(3), domain(4), nplotpts), ...
                        linspace(domain(5), domain(6), nplotpts));
t.x = [xx(:) yy(:) zz(:)]';
so.sc = [0;0;0]; 
phi = inoutfun(so,side); 
iout = phi(t.x(1,:),t.x(2,:),t.x(3,:));
t.x = t.x(:,logical(iout));

% quadrature
tic
count = 0;
tmpidx = 1:numel(t.x(1,:));
len = nvr;
K_Atcxvec = zeros(numel(t.x(1,:)),1);
% use_nearroot = false;
use_nearroot = true;
for k=1:ntri
% for k=13
  sjx = x(:,:,k);
  sjn = nx(:,:,k);
  sjw = w(:,k).';
  sjxbd = xbd(:,:,k);
  sj = struct('x', sjx, 'nx', sjn, 'w', sjw);

  % targets close to this patch
  % qradii = 2.10*sqrt(sum(sj.w)); 
  qradii = 1.75*sqrt(sum(sjw));
  qpoint = mean(sjx,2)';
  idxc = sum((t.x-qpoint').^2)<sum(qradii.^2);
  idxc = tmpidx(idxc);
  tcj.x = t.x(:,idxc);
  ntc = length(tcj.x(1,:));

  if numel(tcj.x)
    % k
    count = count + 1;
    % replace with accurate evaluation
    js = (k-1)*len+(1:len); 
    
    %
    IalphaAsvestas = zeros(ntc,1);
    IalphaAsvestas = evaluate_solid_angle_integral_mex(ntc, tcj.x, len, sjx, sjn, sjw, tri_vert(:,:,k), 3*nquad_bdry, sjxbd, use_nearroot, IalphaAsvestas); 
    K_ij_naive = Lap3dDLPmat(tcj,sj);
    K_Atcxvec(idxc) = K_Atcxvec(idxc) - 1/(4*pi)*IalphaAsvestas - K_ij_naive*ones(len,1);
    % keyboard
  end

end

loop_time = toc; % Stop timer
fprintf('Serial loop time: %f seconds\n', loop_time);
fmm_eps = 1e-15;
u = Lap3dDLPfmm(t,s,-ones(size(sx(1,:)))',fmm_eps) + K_Atcxvec*(-1);

% visualization
err1 = abs(u);
v = NaN(nplotpts,nplotpts,nplotpts);
v(logical(iout)) = log10(err1);

figure(3),clf,
[h,handles,xx,yy,zz,holdState] = myplot_init(domain,nplotpts);
nSteps = 15;
[h,handles] = myplot(domain,nplotpts,h,handles,xx,yy,zz,holdState,v,nSteps);
axis off

keyboard

function phi = inoutfun(so,side)
my_eps = 1e-8;
if side == 'e'
  phi = @(x,y,z) (x-so.sc(1)).^2 + ((y-so.sc(2))*so.ratio).^2 + ((z-so.sc(3))*so.ratio).^2 > 1+my_eps; 
else
  phi = @(x,y,z) (x-so.sc(1)).^2 + ((y-so.sc(2))*so.ratio).^2 + ((z-so.sc(3))*so.ratio).^2 < 1-my_eps;
end
end
