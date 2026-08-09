
rootdir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootdir, '../../matlab'))
addpath(fullfile(rootdir, '../../utils'))
addpath('~/git/FMM3D/matlab/')

do_plot = 1;

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

N = ntri * nvr;
sx = reshape(x, 3, []);
snx = reshape(nx, 3, []);
sw = reshape(w, 1, []);
s.x = sx; s.nx = snx; s.w = sw;
so.ratio = ratio;

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

count = 0;
tmpidx = 1:numel(t.x(1,:));
len = nvr;
K_base = zeros(numel(t.x(1,:)),1);
K_fast = zeros(numel(t.x(1,:)),1);
use_nearroot = true;
t_base = 0;
t_fast = 0;
worst_pair = 0;
worst_abs = 0;
npair = 0;
for k=1:ntri
  sjx = x(:,:,k);
  sjn = nx(:,:,k);
  sjw = w(:,k).';
  sjxbd = xbd(:,:,k);
  sj = struct('x', sjx, 'nx', sjn, 'w', sjw);

  qradii = 1.75*sqrt(sum(sjw));
  qpoint = mean(sjx,2)';
  idxc = sum((t.x-qpoint').^2)<sum(qradii.^2);
  idxc = tmpidx(idxc);
  tcj.x = t.x(:,idxc);
  ntc = length(tcj.x(1,:));

  if numel(tcj.x)
    count = count + 1;
    js = (k-1)*len+(1:len);

    tic
    Ib = zeros(ntc,1);
    Ib = evaluate_solid_angle_integral_mex(ntc, tcj.x, len, sjx, sjn, sjw, tri_vert(:,:,k), 3*nquad_bdry, sjxbd, use_nearroot, Ib);
    t_base = t_base + toc;

    tic
    Ia = zeros(ntc,1);
    Ia = lqs_evaluate_solid_angle_integral_fast_driver_mex(ntc, tcj.x, len, sjx, sjn, sjw, tri_vert(:,:,k), 3*nquad_bdry, sjxbd, Ia);
    t_fast = t_fast + toc;

    worst_pair = max(worst_pair, max(abs(Ia - Ib)./max(abs(Ib), 1e-8)));
    worst_abs  = max(worst_abs,  max(abs(Ia - Ib)));
    npair = npair + ntc;

    K_ij_naive = Lap3dDLPmat(tcj,sj);
    K_base(idxc) = K_base(idxc) - 1/(4*pi)*Ib - K_ij_naive*ones(len,1);
    K_fast(idxc) = K_fast(idxc) - 1/(4*pi)*Ia - K_ij_naive*ones(len,1);
  end

end

fmm_eps = 1e-15;
u_far  = Lap3dDLPfmm(t,s,-ones(size(sx(1,:)))',fmm_eps);
u_base = u_far + K_base*(-1);
u_fast = u_far + K_fast*(-1);

fprintf('=== test_solid_angle_fast ===\n');
fprintf('  ntri = %d   patches with near targets = %d   (triangle,target) pairs = %d\n', ...
        ntri, count, npair);
fprintf('  max |fast - base| on IalphaAsvestas : %.4e abs   %.4e rel\n', worst_abs, worst_pair);
fprintf('\n  DLP error over %d exterior targets\n', numel(u_base));
fprintf('    base   max |u| = %.4e    l2 |u| = %.4e\n', max(abs(u_base)), norm(u_base)/sqrt(numel(u_base)));
fprintf('    fast   max |u| = %.4e    l2 |u| = %.4e\n', max(abs(u_fast)), norm(u_fast)/sqrt(numel(u_fast)));
fprintf('    max |u_fast - u_base| = %.4e\n', max(abs(u_fast - u_base)));
fprintf('\n  solid-angle loop time (serial)\n');
fprintf('    base = %.3f s    fast = %.3f s\n', t_base, t_fast);

if do_plot
  err1 = abs(u_fast);
  v = NaN(nplotpts,nplotpts,nplotpts);
  v(logical(iout)) = log10(err1);
  figure(3),clf,
  [h,handles,xx,yy,zz,holdState] = myplot_init(domain,nplotpts);
  nSteps = 15;
  [h,handles] = myplot(domain,nplotpts,h,handles,xx,yy,zz,holdState,v,nSteps);
  axis off
end

function phi = inoutfun(so,side)
my_eps = 1e-8;
if side == 'e'
  phi = @(x,y,z) (x-so.sc(1)).^2 + ((y-so.sc(2))*so.ratio).^2 + ((z-so.sc(3))*so.ratio).^2 > 1+my_eps;
else
  phi = @(x,y,z) (x-so.sc(1)).^2 + ((y-so.sc(2))*so.ratio).^2 + ((z-so.sc(3))*so.ratio).^2 < 1-my_eps;
end
end
