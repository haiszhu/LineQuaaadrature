function [u,un,unx,uny,unz] = Lap3dDLPfmm(t,s,tau,fmm_eps)
% LAP3DDLPFMM. Laplace DLP fmm from sources to targets 
% similar interface as Lap3dDLPmat & Lap3dDLP_closepanel

% 
srcinfo.sources = s.x;
srcinfo.dipoles = bsxfun(@times,s.nx,tau(:)'.*s.w(:)'); % srcinfo = rmfield(srcinfo,'charges'); if necessary
targ = t.x;
if nargin < 4, fmm_eps = 1e-15; end
ifppreg = 0;

if nargout < 2
  ifppregtarg = 1; % potential at target
  U = lfmm3d(fmm_eps,srcinfo,ifppreg,targ,ifppregtarg);
  u = (U.pottarg)'; 
elseif nargout < 3 % need t.nx...
  ifppregtarg = 2; % potential and its gradient at target
  U = lfmm3d(fmm_eps,srcinfo,ifppreg,targ,ifppregtarg);
  u = (U.pottarg)'; 
  un = (sum(t.nx.*U.gradtarg))'; 
else
  ifppregtarg = 2; % potential and its gradient at target
  U = lfmm3d(fmm_eps,srcinfo,ifppreg,targ,ifppregtarg);
  u = (U.pottarg)'; 
  un = (sum(t.nx.*U.gradtarg))'; 
  unx = (U.gradtarg(1,:))';
  uny = (U.gradtarg(2,:))';
  unz = (U.gradtarg(3,:))';
end

end