function [A,An,Anx,Any,Anz] = Lap3dSLPmat(t,s)

d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:));
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
rr = d1.^2+d2.^2+d3.^2;
weights = s.w(:)';
A = bsxfun(@times,1./sqrt(rr),weights*(1/(4*pi)));

if nargout > 1
  ddottn = bsxfun(@times,d1,t.nx(1,:)') + ...
      bsxfun(@times,d2,t.nx(2,:)') + ...
      bsxfun(@times,d3,t.nx(3,:)');
  An = bsxfun(@times,ddottn./(sqrt(rr).*rr),weights*(-1/(4*pi)));
  if nargout > 2
    Anx = bsxfun(@times,d1./(sqrt(rr).*rr),weights*(-1/(4*pi)));
    Any = bsxfun(@times,d2./(sqrt(rr).*rr),weights*(-1/(4*pi)));
    Anz = bsxfun(@times,d3./(sqrt(rr).*rr),weights*(-1/(4*pi)));
  end
end
end
