function [u,un,unx,uny,unz] = Lap3dSLPfmm(t,s,tau,fmm_eps,Backend,Context)

srcinfo.sources = s.x;
srcinfo.charges = tau(:)'.*s.w(:)';
targ = t.x;
if nargin < 4 || isempty(fmm_eps), fmm_eps = 1e-15; end
if nargin < 5 || isempty(Backend), Backend = 'fmm3d'; end
ifppreg = 0;

if nargout < 2
  ifppregtarg = 1;
  if strcmpi(Backend, 'pvfmm')
    if nargin < 6 || isempty(Context) || ~isfield(Context, 'MPI_HANDLE')
      error('MPI_Context with a valid handle is required for pvfmm backend.');
    end
    if exist('laplace_slp_fmm_mex', 'file') ~= 2
      error('laplace_slp_fmm_mex not found. Build pvfmm MEX first.');
    end
    ns = size(s.x, 2);
    nt = size(t.x, 2);
    u = laplace_slp_fmm_mex(Context.MPI_HANDLE, s.x(:), ...
        tau(:).*s.w(:), t.x(:), ns, nt, Context.MULT_ORDER, ...
        Context.MEM_SIZE, Context.MAX_PTS, zeros(nt,1));
  else
    U = lfmm3d(fmm_eps,srcinfo,ifppreg,targ,ifppregtarg);
    u = U.pottarg';
  end
elseif nargout < 3
  if strcmpi(Backend, 'pvfmm')
    error('pvfmm backend currently supports only nargout < 2.');
  end
  ifppregtarg = 2;
  U = lfmm3d(fmm_eps,srcinfo,ifppreg,targ,ifppregtarg);
  u = U.pottarg';
  un = sum(t.nx.*U.gradtarg)';
else
  if strcmpi(Backend, 'pvfmm')
    error('pvfmm backend currently supports only nargout < 2.');
  end
  ifppregtarg = 2;
  U = lfmm3d(fmm_eps,srcinfo,ifppreg,targ,ifppregtarg);
  u = U.pottarg';
  un = sum(t.nx.*U.gradtarg)';
  unx = U.gradtarg(1,:)';
  uny = U.gradtarg(2,:)';
  unz = U.gradtarg(3,:)';
end
end
