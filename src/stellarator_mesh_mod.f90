module stellarator_mesh_mod
  use linequaaadrature_mod, only: gauss_r64, r64
  implicit none
  private

  real(r64), parameter :: ASPECT = 1.65_r64
  integer(8), parameter, public :: W7X_NFP = 5_8
  integer(8), parameter, public :: W7X_NMODE = 288_8

  type :: chart_t
    real(r64) :: a1, b1, a2, b2
    real(r64) :: jj, ii
  end type chart_t

  public :: stellarator_mesh_init_r64
  public :: create_stellarator_tri_mesh_r64
  public :: stellarator_tri_uv2x_r64
  public :: get_w7x_modes_r64

contains

  subroutine get_w7x_modes_r64(nfp, nmode, mn, rc, zs)
    integer(8), intent(out) :: nfp, nmode
    integer(8), intent(out) :: mn(2*W7X_NMODE)
    real(r64), intent(out) :: rc(W7X_NMODE), zs(W7X_NMODE)

    nfp = W7X_NFP
    nmode = W7X_NMODE
    INCLUDE 'w7x-modes-dat.txt'
  end subroutine get_w7x_modes_r64

  subroutine stellarator_mesh_init_r64(mp, np, p, nfp, nmode, mn, rc, zs, &
                                       restol, cap, charts, nchart, ntri, ier)
    integer(8), intent(in) :: mp, np, p, nfp, nmode, cap
    integer(8), intent(in) :: mn(2,nmode)
    real(r64), intent(in) :: rc(nmode), zs(nmode), restol
    real(r64), intent(out) :: charts(6,cap)
    integer(8), intent(out) :: nchart, ntri, ier

    integer(8) :: i
    real(r64) :: pi, ts1, ts2
    real(r64), allocatable :: x0(:), w0(:), D(:,:)
    type(chart_t), allocatable :: work(:)

    charts = 0.0_r64
    nchart = 0_8
    ntri = 0_8
    ier = 1_8
    if (mp < 1_8 .or. np < 1_8 .or. p < 1_8 .or. p > 20_8) return
    if (nfp < 1_8 .or. nmode < 0_8 .or. cap < 1_8) return
    pi = acos(-1.0_r64)
    ts1 = pi/real(np, r64)
    ts2 = pi/real(mp, r64)
    allocate(x0(p), w0(p), D(p,p))
    allocate(work(max(cap,2_8*mp*np)))
    call gauss_r64(p, x0, w0, D)
    if (nmode == 0_8) then
      call build_charts_r64(mp, np, p, x0, ts1, ts2, work, nchart)
    else
      call build_charts_adaptive_r64(mp,np,p,x0,D,ts1,ts2,nfp,nmode, &
           mn,rc,zs,restol,cap,work,nchart,ier)
      if (ier /= 0_8) return
    end if
    if (nchart > cap) then
      nchart = 0_8
      ier = 2_8
      return
    end if
    do i = 1_8, nchart
      call pack_chart_r64(work(i), charts(:,i))
    end do
    ntri = 2_8*nchart
    ier = 0_8
  end subroutine stellarator_mesh_init_r64

  subroutine create_stellarator_tri_mesh_r64(mp, np, p, nfp, nmode, mn, rc, zs, &
      nchart, charts, ntri, sx, snx, sw, rts, rps, ier)
    integer(8), intent(in) :: mp, np, p, nfp, nmode, nchart, ntri
    integer(8), intent(in) :: mn(2,nmode)
    real(r64), intent(in) :: rc(nmode), zs(nmode), charts(6,nchart)
    real(r64), intent(out) :: sx(3,ntri*p*(p+1_8)/2_8)
    real(r64), intent(out) :: snx(3,ntri*p*(p+1_8)/2_8)
    real(r64), intent(out) :: sw(ntri*p*(p+1_8)/2_8)
    real(r64), intent(out) :: rts(3,ntri*p*(p+1_8)/2_8)
    real(r64), intent(out) :: rps(3,ntri*p*(p+1_8)/2_8)
    integer(8), intent(out) :: ier

    integer(8) :: hdim, i
    real(r64), allocatable :: x0(:), w0(:), D(:,:)
    real(r64), allocatable :: uvs(:,:), wts(:)
    type(chart_t), allocatable :: work(:)

    ier = 1_8
    if (mp < 1_8 .or. np < 1_8 .or. p < 1_8 .or. p > 20_8) return
    if (nfp < 1_8 .or. nmode < 0_8 .or. nchart < 1_8) return
    if (ntri /= 2_8*nchart) return

    hdim = p*(p+1_8)/2_8
    allocate(x0(p), w0(p), D(p,p), uvs(2,hdim), wts(hdim), work(nchart))
    call gauss_r64(p, x0, w0, D)
    call get_vr_nodes_wts_r64(p-1_8, hdim, uvs, wts)
    do i = 1_8, nchart
      call unpack_chart_r64(charts(:,i), work(i))
    end do
    call build_mesh_r64(work,nchart,mp,np,p,x0,D,uvs,wts,nfp,nmode,mn,rc,zs, &
                        sx, snx, sw, rts, rps)
    ier = 0_8
  end subroutine create_stellarator_tri_mesh_r64

  subroutine stellarator_tri_uv2x_r64(mp, np, p, nfp, nmode, mn, rc, zs, &
      nchart, charts, itri, nuv, uv, x, ier)
    integer(8), intent(in) :: mp, np, p, nfp, nmode, nchart, itri, nuv
    integer(8), intent(in) :: mn(2,nmode)
    real(r64), intent(in) :: rc(nmode), zs(nmode), charts(6,nchart)
    real(r64), intent(in) :: uv(2,nuv)
    real(r64), intent(out) :: x(3,nuv)
    integer(8), intent(out) :: ier

    integer(8) :: ipan, side, i, idim
    real(r64) :: pi, ts1, ts2, d1, d2, t1, t2
    real(r64) :: xll(3), xul(3), xlr(3), xur(3), delta(3), vertices(2,3)
    real(r64), allocatable :: x0(:), w0(:), D(:,:)
    type(chart_t) :: chart

    x = 0.0_r64
    ier = 1_8
    if (mp < 1_8 .or. np < 1_8 .or. p < 1_8 .or. p > 20_8) return
    if (nfp < 1_8 .or. nmode < 0_8 .or. nchart < 1_8 .or. nuv < 0_8) return
    if (itri < 1_8 .or. itri > 2_8*nchart) return

    ipan = (itri-1_8)/2_8+1_8
    side = modulo(itri-1_8,2_8)
    call unpack_chart_r64(charts(:,ipan), chart)
    allocate(x0(p), w0(p), D(p,p))
    call gauss_r64(p, x0, w0, D)
    pi = acos(-1.0_r64)
    ts1 = pi/real(np,r64)
    ts2 = pi/real(mp,r64)

    call chart_eval_surface_r64(chart,ts1,ts2,x0(1),x0(1), &
         nfp,nmode,mn,rc,zs,xll)
    call chart_eval_surface_r64(chart,ts1,ts2,x0(p),x0(p), &
         nfp,nmode,mn,rc,zs,xur)
    call chart_eval_surface_r64(chart,ts1,ts2,x0(1),x0(p), &
         nfp,nmode,mn,rc,zs,xul)
    call chart_eval_surface_r64(chart,ts1,ts2,x0(p),x0(1), &
         nfp,nmode,mn,rc,zs,xlr)
    d1 = 0.0_r64
    d2 = 0.0_r64
    do idim = 1_8, 3_8
      delta(idim) = xll(idim)-xur(idim)
      d1 = d1+delta(idim)*delta(idim)
      delta(idim) = xul(idim)-xlr(idim)
      d2 = d2+delta(idim)*delta(idim)
    end do
    if (d1+1.0e-13_r64 > d2) then
      if (side == 0_8) then
        vertices = reshape([-1.0_r64,-1.0_r64, 1.0_r64,-1.0_r64, &
                            -1.0_r64,1.0_r64],[2,3])
      else
        vertices = reshape([1.0_r64,1.0_r64, -1.0_r64,1.0_r64, &
                            1.0_r64,-1.0_r64],[2,3])
      end if
    else
      if (side == 0_8) then
        vertices = reshape([-1.0_r64,1.0_r64, -1.0_r64,-1.0_r64, &
                            1.0_r64,1.0_r64],[2,3])
      else
        vertices = reshape([1.0_r64,-1.0_r64, 1.0_r64,1.0_r64, &
                            -1.0_r64,-1.0_r64],[2,3])
      end if
    end if
    do i = 1_8, nuv
      t1 = vertices(1,1)+(vertices(1,2)-vertices(1,1))*uv(1,i) &
           +(vertices(1,3)-vertices(1,1))*uv(2,i)
      t2 = vertices(2,1)+(vertices(2,2)-vertices(2,1))*uv(1,i) &
           +(vertices(2,3)-vertices(2,1))*uv(2,i)
      call chart_eval_surface_r64(chart,ts1,ts2,t1,t2,nfp,nmode,mn,rc,zs,x(:,i))
    end do
    ier = 0_8
  end subroutine stellarator_tri_uv2x_r64

  subroutine pack_chart_r64(chart, packed)
    type(chart_t), intent(in) :: chart
    real(r64), intent(out) :: packed(6)
    packed = [chart%a1,chart%b1,chart%a2,chart%b2,chart%jj,chart%ii]
  end subroutine pack_chart_r64

  subroutine unpack_chart_r64(packed, chart)
    real(r64), intent(in) :: packed(6)
    type(chart_t), intent(out) :: chart
    chart%a1 = packed(1)
    chart%b1 = packed(2)
    chart%a2 = packed(3)
    chart%b2 = packed(4)
    chart%jj = packed(5)
    chart%ii = packed(6)
  end subroutine unpack_chart_r64

  subroutine get_vr_nodes_wts_r64(norder, npols, uvs, wts)
    integer(8), intent(in) :: norder, npols
    real(r64), intent(out) :: uvs(2,npols), wts(npols)

    INCLUDE 'koorn-uvs-dat.txt'
    INCLUDE 'koorn-wts-dat.txt'
  end subroutine get_vr_nodes_wts_r64

  subroutine stellarator_param_r64(t, p, x)
    real(r64), intent(in) :: t, p
    real(r64), intent(out) :: x(3)

    real(r64) :: ct, st, c2pt, s2pt, c2p, s2p
    real(r64) :: cp, sp, cpt, spt, cmpt, smpt

    ct = cos(t)
    st = sin(t)
    c2pt = cos(2.0_r64*p-t)
    s2pt = sin(2.0_r64*p-t)
    c2p = cos(2.0_r64*p)
    s2p = sin(2.0_r64*p)
    cp = cos(p)
    sp = sin(p)
    cpt = cos(p+t)
    spt = sin(p+t)
    cmpt = cos(-p+t)
    smpt = sin(-p+t)

    x(1) = 0.17_r64*ct*c2pt + 0.11_r64*ct*c2p + 1.0_r64*ct*cp + 4.5_r64*ct &
         - 0.25_r64*ct*cp + 0.07_r64*ct*cpt - 0.45_r64*ct*cmpt
    x(2) = 0.17_r64*st*c2pt + 0.11_r64*st*c2p + 1.0_r64*st*cp + 4.5_r64*st &
         - 0.25_r64*st*cp + 0.07_r64*st*cpt - 0.45_r64*st*cmpt
    x(3) = 0.17_r64*s2pt + 0.11_r64*s2p + 1.0_r64*sp + 0.0_r64 &
         + 0.25_r64*sp + 0.07_r64*spt - 0.45_r64*smpt
  end subroutine stellarator_param_r64

  subroutine surface_eval_r64(nfp, nmode, mn, rc, zs, t, p, x)
    integer(8), intent(in) :: nfp, nmode, mn(2,nmode)
    real(r64), intent(in) :: rc(nmode), zs(nmode), t, p
    real(r64), intent(out) :: x(3)

    integer(8) :: k
    real(r64) :: radius, height, angle

    if (nmode == 0_8) then
      call stellarator_param_r64(t,p,x)
      return
    end if
    radius = 0.0_r64
    height = 0.0_r64
    do k = 1_8, nmode
      angle = real(mn(1,k),r64)*p-real(nfp*mn(2,k),r64)*t
      radius = radius+rc(k)*cos(angle)
      height = height+zs(k)*sin(angle)
    end do
    x(1) = radius*cos(t)
    x(2) = radius*sin(t)
    x(3) = height
  end subroutine surface_eval_r64

  subroutine chart_eval_surface_r64(chart,ts1,ts2,t1,t2,nfp,nmode,mn,rc,zs,x)
    type(chart_t), intent(in) :: chart
    real(r64), intent(in) :: ts1,ts2,t1,t2
    integer(8), intent(in) :: nfp,nmode,mn(2,nmode)
    real(r64), intent(in) :: rc(nmode),zs(nmode)
    real(r64), intent(out) :: x(3)

    real(r64) :: u,v

    u=((chart%a1*t1+chart%b1)+2.0_r64*chart%jj)*ts1
    v=((chart%a2*t2+chart%b2)+2.0_r64*chart%ii)*ts2
    call surface_eval_r64(nfp,nmode,mn,rc,zs,u,v,x)
  end subroutine chart_eval_surface_r64

  subroutine chart_eval_r64(chart, ts1, ts2, t1, t2, x)
    type(chart_t), intent(in) :: chart
    real(r64), intent(in) :: ts1, ts2, t1, t2
    real(r64), intent(out) :: x(3)

    real(r64) :: u, v

    u = ((chart%a1*t1+chart%b1)+2.0_r64*chart%jj)*ts1
    v = ((chart%a2*t2+chart%b2)+2.0_r64*chart%ii)*ts2
    call stellarator_param_r64(u, v, x)
  end subroutine chart_eval_r64

  subroutine build_charts_r64(mp, np, p, x0, ts1, ts2, charts, nchart)
    integer(8), intent(in) :: mp, np, p
    real(r64), intent(in) :: x0(p), ts1, ts2
    type(chart_t), intent(out) :: charts(2_8*mp*np)
    integer(8), intent(out) :: nchart

    integer(8) :: i, j
    real(r64) :: xll(3), xul(3), xlr(3), dt, dp
    type(chart_t) :: base

    base%a1 = 1.0_r64
    base%b1 = 0.0_r64
    base%a2 = 1.0_r64
    base%b2 = 0.0_r64
    nchart = 0_8

    do j = 1_8, np
      do i = 1_8, mp
        base%jj = real(j, r64)
        base%ii = real(i, r64)
        call chart_eval_r64(base, ts1, ts2, x0(1), x0(1), xll)
        call chart_eval_r64(base, ts1, ts2, x0(1), x0(p), xul)
        call chart_eval_r64(base, ts1, ts2, x0(p), x0(1), xlr)
        dt = sqrt(sum((xll-xlr)**2))
        dp = sqrt(sum((xll-xul)**2))

        if (dt/dp > ASPECT) then
          nchart = nchart+1_8
          charts(nchart) = base
          charts(nchart)%a1 = 0.5_r64
          charts(nchart)%b1 = -0.5_r64
          nchart = nchart+1_8
          charts(nchart) = base
          charts(nchart)%a1 = 0.5_r64
          charts(nchart)%b1 = 0.5_r64
        else if (dp/dt > ASPECT) then
          nchart = nchart+1_8
          charts(nchart) = base
          charts(nchart)%a2 = 0.5_r64
          charts(nchart)%b2 = -0.5_r64
          nchart = nchart+1_8
          charts(nchart) = base
          charts(nchart)%a2 = 0.5_r64
          charts(nchart)%b2 = 0.5_r64
        else
          nchart = nchart+1_8
          charts(nchart) = base
        end if
      end do
    end do
  end subroutine build_charts_r64

  subroutine chart_box_r64(chart,ts1,ts2,box)
    type(chart_t),intent(in)::chart
    real(r64),intent(in)::ts1,ts2
    real(r64),intent(out)::box(4)
    box(1)=((chart%b1-chart%a1)+2.0_r64*chart%jj)*ts1
    box(2)=((chart%b1+chart%a1)+2.0_r64*chart%jj)*ts1
    box(3)=((chart%b2-chart%a2)+2.0_r64*chart%ii)*ts2
    box(4)=((chart%b2+chart%a2)+2.0_r64*chart%ii)*ts2
  end subroutine chart_box_r64

  subroutine chart_grid_r64(chart,ts1,ts2,p,x0,nfp,nmode,mn,rc,zs,sxk)
    type(chart_t),intent(in)::chart
    integer(8),intent(in)::p,nfp,nmode,mn(2,nmode)
    real(r64),intent(in)::ts1,ts2,x0(p),rc(nmode),zs(nmode)
    real(r64),intent(out)::sxk(3,p*p)
    integer(8)::a,b,m
    do b=1_8,p
      do a=1_8,p
        m=a+(b-1_8)*p
        call chart_eval_surface_r64(chart,ts1,ts2,x0(b),x0(a), &
             nfp,nmode,mn,rc,zs,sxk(:,m))
      end do
    end do
  end subroutine chart_grid_r64

  subroutine diff_a_r64(D,p,field,derivative)
    integer(8),intent(in)::p
    real(r64),intent(in)::D(p,p),field(3,p*p)
    real(r64),intent(out)::derivative(3,p*p)
    integer(8)::a,b,k,idim
    real(r64)::value
    do b=1_8,p
      do a=1_8,p
        do idim=1_8,3_8
          value=0.0_r64
          do k=1_8,p
            value=value+D(a,k)*field(idim,k+(b-1_8)*p)
          end do
          derivative(idim,a+(b-1_8)*p)=value
        end do
      end do
    end do
  end subroutine diff_a_r64

  subroutine diff_b_r64(D,p,field,derivative)
    integer(8),intent(in)::p
    real(r64),intent(in)::D(p,p),field(3,p*p)
    real(r64),intent(out)::derivative(3,p*p)
    integer(8)::a,b,k,idim
    real(r64)::value
    do b=1_8,p
      do a=1_8,p
        do idim=1_8,3_8
          value=0.0_r64
          do k=1_8,p
            value=value+D(b,k)*field(idim,a+(k-1_8)*p)
          end do
          derivative(idim,a+(b-1_8)*p)=value
        end do
      end do
    end do
  end subroutine diff_b_r64

  subroutine legendre_vandermonde_r64(p,x0,V)
    integer(8),intent(in)::p
    real(r64),intent(in)::x0(p)
    real(r64),intent(out)::V(p,p)
    integer(8)::i,k
    do i=1_8,p
      V(i,1)=1.0_r64
      if(p>1_8)V(i,2)=x0(i)
      do k=2_8,p-1_8
        V(i,k+1)=((2.0_r64*real(k-1_8,r64)+1.0_r64)*x0(i)*V(i,k) &
             -real(k-1_8,r64)*V(i,k-1))/real(k,r64)
      end do
    end do
  end subroutine legendre_vandermonde_r64

  function panel_resolution_r64(chart,ts1,ts2,p,x0,D,nfp,nmode,mn,rc,zs,V) &
      result(resolution)
    type(chart_t),intent(in)::chart
    integer(8),intent(in)::p,nfp,nmode,mn(2,nmode)
    real(r64),intent(in)::ts1,ts2,x0(p),D(p,p),rc(nmode),zs(nmode),V(p,p)
    real(r64)::resolution
    integer(8)::n2,k,a,b
    real(r64)::normal(3),nn,E,F,G,L,M,Ncoef,den,z,tail,total
    real(r64),allocatable::sxk(:,:),r1(:,:),r2(:,:),r11(:,:),r12(:,:),r22(:,:)
    real(r64),allocatable::H(:),Ccoef(:,:),Tmat(:,:),Awork(:,:)

    n2=p*p
    allocate(sxk(3,n2),r1(3,n2),r2(3,n2),r11(3,n2),r12(3,n2),r22(3,n2))
    allocate(H(n2),Ccoef(p,p),Tmat(p,p),Awork(p,p))
    call chart_grid_r64(chart,ts1,ts2,p,x0,nfp,nmode,mn,rc,zs,sxk)
    call diff_b_r64(D,p,sxk,r1)
    call diff_a_r64(D,p,sxk,r2)
    call diff_b_r64(D,p,r1,r11)
    call diff_b_r64(D,p,r2,r12)
    call diff_a_r64(D,p,r2,r22)

    do k=1_8,n2
      normal(1)=r1(2,k)*r2(3,k)-r1(3,k)*r2(2,k)
      normal(2)=r1(3,k)*r2(1,k)-r1(1,k)*r2(3,k)
      normal(3)=r1(1,k)*r2(2,k)-r1(2,k)*r2(1,k)
      nn=sqrt(sum(normal*normal))
      if(nn<=0.0_r64)nn=1.0_r64
      normal=normal/nn
      E=sum(r1(:,k)*r1(:,k))
      F=sum(r1(:,k)*r2(:,k))
      G=sum(r2(:,k)*r2(:,k))
      L=sum(normal*r11(:,k))
      M=sum(normal*r12(:,k))
      Ncoef=sum(normal*r22(:,k))
      den=2.0_r64*(E*G-F*F)
      if(den==0.0_r64)den=1.0_r64
      H(k)=(E*Ncoef-2.0_r64*F*M+G*L)/den
    end do

    Ccoef=reshape(H,[p,p])
    Awork=V
    call lu_solve_r64(Awork,p,Ccoef,p)
    Tmat=transpose(Ccoef)
    Awork=V
    call lu_solve_r64(Awork,p,Tmat,p)
    Ccoef=transpose(Tmat)
    tail=0.0_r64
    total=0.0_r64
    do a=1_8,p
      do b=1_8,p
        z=Ccoef(a,b)
        total=total+z*z
        if(a>=p-1_8 .or. b>=p-1_8)tail=tail+z*z
      end do
    end do
    resolution=sqrt(tail/merge(total,1.0_r64,total>0.0_r64))
  end function panel_resolution_r64

  subroutine split_chart4_r64(chart,children)
    type(chart_t),intent(in)::chart
    type(chart_t),intent(out)::children(4)
    integer(8)::q
    do q=0_8,3_8
      children(q+1)=chart
      children(q+1)%a1=chart%a1/2.0_r64
      children(q+1)%a2=chart%a2/2.0_r64
      children(q+1)%b1=chart%b1+merge(chart%a1/2.0_r64, &
           -chart%a1/2.0_r64,btest(q,0))
      children(q+1)%b2=chart%b2+merge(chart%a2/2.0_r64, &
           -chart%a2/2.0_r64,btest(q,1))
    end do
  end subroutine split_chart4_r64

  subroutine build_charts_adaptive_r64(mp,np,p,x0,D,ts1,ts2,nfp,nmode,mn, &
      rc,zs,restol,cap,charts,nchart,ier)
    integer(8),intent(in)::mp,np,p,nfp,nmode,mn(2,nmode),cap
    real(r64),intent(in)::x0(p),D(p,p),ts1,ts2,rc(nmode),zs(nmode),restol
    type(chart_t),intent(out)::charts(cap)
    integer(8),intent(out)::nchart,ier
    integer(8)::i,j,k,q,pass,nsplit,n0
    logical,allocatable::done(:)
    real(r64),allocatable::V(:,:)
    real(r64)::xll(3),xul(3),xlr(3),dt,dp,boxk(4),boxq(4)
    real(r64)::du,dv,edgek,edgeq,dj,di
    type(chart_t)::current,children(4)

    nchart=0_8
    ier=2_8
    allocate(done(cap),V(p,p))
    call legendre_vandermonde_r64(p,x0,V)
    do j=1_8,np
      do i=1_8,mp
        if(nchart>=cap)return
        nchart=nchart+1_8
        charts(nchart)%a1=1.0_r64
        charts(nchart)%b1=0.0_r64
        charts(nchart)%a2=1.0_r64
        charts(nchart)%b2=0.0_r64
        charts(nchart)%jj=real(j,r64)
        charts(nchart)%ii=real(i,r64)
      end do
    end do

    do pass=1_8,8_8
      nsplit=0_8
      k=1_8
      do while(k<=nchart)
        current=charts(k)
        call chart_eval_surface_r64(current,ts1,ts2,x0(1),x0(1), &
             nfp,nmode,mn,rc,zs,xll)
        call chart_eval_surface_r64(current,ts1,ts2,x0(1),x0(p), &
             nfp,nmode,mn,rc,zs,xul)
        call chart_eval_surface_r64(current,ts1,ts2,x0(p),x0(1), &
             nfp,nmode,mn,rc,zs,xlr)
        dt=sqrt(sum((xll-xlr)**2))
        dp=sqrt(sum((xll-xul)**2))
        if(dt/dp>ASPECT)then
          if(nchart>=cap)return
          charts(k)%a1=current%a1/2.0_r64
          charts(k)%b1=current%b1-current%a1/2.0_r64
          nchart=nchart+1_8
          charts(nchart)=current
          charts(nchart)%a1=current%a1/2.0_r64
          charts(nchart)%b1=current%b1+current%a1/2.0_r64
          nsplit=nsplit+1_8
        else if(dp/dt>ASPECT)then
          if(nchart>=cap)return
          charts(k)%a2=current%a2/2.0_r64
          charts(k)%b2=current%b2-current%a2/2.0_r64
          nchart=nchart+1_8
          charts(nchart)=current
          charts(nchart)%a2=current%a2/2.0_r64
          charts(nchart)%b2=current%b2+current%a2/2.0_r64
          nsplit=nsplit+1_8
        end if
        k=k+1_8
      end do
      if(nsplit==0_8)exit
    end do

    if(restol>0.0_r64)then
      done(1:nchart)=.false.
      do pass=1_8,12_8
        n0=nchart
        nsplit=0_8
        do k=1_8,n0
          if(done(k))cycle
          if(panel_resolution_r64(charts(k),ts1,ts2,p,x0,D,nfp,nmode, &
               mn,rc,zs,V)<=restol)then
            done(k)=.true.
            cycle
          end if
          if(nchart+3_8>cap)return
          call split_chart4_r64(charts(k),children)
          charts(k)=children(1)
          done(k)=.false.
          do q=2_8,4_8
            nchart=nchart+1_8
            charts(nchart)=children(q)
            done(nchart)=.false.
          end do
          nsplit=nsplit+1_8
        end do
        if(nsplit==0_8)exit
      end do

      do pass=1_8,20_8
        n0=nchart
        nsplit=0_8
        done(1:n0)=.false.
        do k=1_8,n0
          call chart_box_r64(charts(k),ts1,ts2,boxk)
          do q=1_8,n0
            if(q==k)cycle
            dj=abs(charts(q)%jj-charts(k)%jj)
            if(dj>1.5_r64 .and. abs(dj-real(np,r64))>1.5_r64)cycle
            di=abs(charts(q)%ii-charts(k)%ii)
            if(di>1.5_r64 .and. abs(di-real(mp,r64))>1.5_r64)cycle
            call chart_box_r64(charts(q),ts1,ts2,boxq)
            du=1.0e-9_r64*ts1
            dv=1.0e-9_r64*ts2
            if((abs(boxk(2)-boxq(1))<du .or. abs(boxk(1)-boxq(2))<du) &
                .and. min(boxk(4),boxq(4))-max(boxk(3),boxq(3))>dv)then
              edgek=boxk(4)-boxk(3)
              edgeq=boxq(4)-boxq(3)
              if(edgek>2.0_r64*edgeq+dv)done(k)=.true.
            else if((abs(boxk(4)-boxq(3))<dv .or. abs(boxk(3)-boxq(4))<dv) &
                .and. min(boxk(2),boxq(2))-max(boxk(1),boxq(1))>du)then
              edgek=boxk(2)-boxk(1)
              edgeq=boxq(2)-boxq(1)
              if(edgek>2.0_r64*edgeq+du)done(k)=.true.
            end if
          end do
        end do
        do k=1_8,n0
          if(.not.done(k))cycle
          if(nchart+3_8>cap)return
          call split_chart4_r64(charts(k),children)
          charts(k)=children(1)
          do q=2_8,4_8
            nchart=nchart+1_8
            charts(nchart)=children(q)
          end do
          nsplit=nsplit+1_8
        end do
        if(nsplit==0_8)exit
      end do
    end if
    ier=0_8
  end subroutine build_charts_adaptive_r64

  subroutine lu_solve_r64(A, p, B, q)
    integer(8), intent(in) :: p, q
    real(r64), intent(inout) :: A(p,p), B(p,q)

    integer(8) :: k, i, j, ipv
    real(r64) :: pmax, tmp, multiplier, sm

    do k = 1_8, p
      ipv = k
      pmax = abs(A(k,k))
      do i = k+1_8, p
        if (abs(A(i,k)) > pmax) then
          pmax = abs(A(i,k))
          ipv = i
        end if
      end do
      if (ipv /= k) then
        do j = 1_8, p
          tmp = A(k,j)
          A(k,j) = A(ipv,j)
          A(ipv,j) = tmp
        end do
        do j = 1_8, q
          tmp = B(k,j)
          B(k,j) = B(ipv,j)
          B(ipv,j) = tmp
        end do
      end if
      do i = k+1_8, p
        multiplier = A(i,k)/A(k,k)
        A(i,k) = multiplier
        do j = k+1_8, p
          A(i,j) = A(i,j)-multiplier*A(k,j)
        end do
        do j = 1_8, q
          B(i,j) = B(i,j)-multiplier*B(k,j)
        end do
      end do
    end do

    do j = 1_8, q
      do i = p, 1_8, -1_8
        sm = B(i,j)
        do k = i+1_8, p
          sm = sm-A(i,k)*B(k,j)
        end do
        B(i,j) = sm/A(i,i)
      end do
    end do
  end subroutine lu_solve_r64

  subroutine interpolation_matrix_r64(t, q, s, p, L)
    integer(8), intent(in) :: q, p
    real(r64), intent(in) :: t(q), s(p)
    real(r64), intent(out) :: L(q,p)

    integer(8) :: i, j
    real(r64), allocatable :: A(:,:), B(:,:)

    allocate(A(p,p), B(p,q))
    do i = 1_8, p
      A(1,i) = 1.0_r64
      do j = 2_8, p
        A(j,i) = A(j-1_8,i)*s(i)
      end do
    end do
    do i = 1_8, q
      B(1,i) = 1.0_r64
      do j = 2_8, p
        B(j,i) = B(j-1_8,i)*t(i)
      end do
    end do
    call lu_solve_r64(A, p, B, q)
    do i = 1_8, q
      do j = 1_8, p
        L(i,j) = B(j,i)
      end do
    end do
  end subroutine interpolation_matrix_r64

  subroutine build_tensor_interpolation_r64(uv, hdim, x0, p, Mt)
    integer(8), intent(in) :: hdim, p
    real(r64), intent(in) :: uv(2,hdim), x0(p)
    real(r64), intent(out) :: Mt(p*p,hdim)

    integer(8) :: i, a, b, m
    real(r64), allocatable :: t1(:), t2(:), L1(:,:), L2(:,:)

    allocate(t1(hdim), t2(hdim), L1(hdim,p), L2(hdim,p))
    t1 = uv(1,:)
    t2 = uv(2,:)
    call interpolation_matrix_r64(t1, hdim, x0, p, L1)
    call interpolation_matrix_r64(t2, hdim, x0, p, L2)
    do i = 1_8, hdim
      do a = 1_8, p
        do b = 1_8, p
          m = b+(a-1_8)*p
          Mt(m,i) = L1(i,a)*L2(i,b)
        end do
      end do
    end do
  end subroutine build_tensor_interpolation_r64

  subroutine build_mesh_r64(charts,nchart,mp,np,p,x0,D,uvs,wts, &
                            nfp,nmode,mn,rc,zs,sx,snx,sw,rts_out,rps_out)
    integer(8), intent(in) :: nchart,mp,np,p,nfp,nmode,mn(2,nmode)
    type(chart_t), intent(in) :: charts(nchart)
    real(r64), intent(in) :: x0(p), D(p,p)
    real(r64), intent(in) :: uvs(2,p*(p+1_8)/2_8), wts(p*(p+1_8)/2_8)
    real(r64), intent(in) :: rc(nmode),zs(nmode)
    real(r64), intent(out) :: sx(3,nchart*p*(p+1_8))
    real(r64), intent(out) :: snx(3,nchart*p*(p+1_8))
    real(r64), intent(out) :: sw(nchart*p*(p+1_8))
    real(r64), intent(out) :: rts_out(3,nchart*p*(p+1_8))
    real(r64), intent(out) :: rps_out(3,nchart*p*(p+1_8))

    integer(8) :: p2, hdim, ipan, i, k, m, a, b, idim, row, col
    integer(8) :: side, off
    logical :: diagonal_a
    real(r64) :: pi, ts1, ts2, sumd, d1, d2, coeff, speed
    real(r64) :: delta(3), acc(9), xt(3), xp(3), normal(3)
    real(r64), allocatable :: uvtmp(:,:)
    real(r64), allocatable, target :: MtlA(:,:), MtrA(:,:), MtlB(:,:), MtrB(:,:)
    real(r64), pointer :: Mt(:,:)
    real(r64), allocatable :: xx1(:), xx2(:), sxk(:,:)
    real(r64), allocatable :: rts(:,:), rps(:,:), S(:,:)

    p2 = p*p
    hdim = p*(p+1_8)/2_8
    pi = acos(-1.0_r64)
    ts1 = pi/real(np, r64)
    ts2 = pi/real(mp, r64)

    allocate(uvtmp(2,hdim))
    allocate(MtlA(p2,hdim), MtrA(p2,hdim))
    allocate(MtlB(p2,hdim), MtrB(p2,hdim))

    uvtmp = 2.0_r64*uvs-1.0_r64
    call build_tensor_interpolation_r64(uvtmp, hdim, x0, p, MtlA)
    uvtmp = -uvtmp
    call build_tensor_interpolation_r64(uvtmp, hdim, x0, p, MtrA)
    do i = 1_8, hdim
      uvtmp(1,i) = 2.0_r64*uvs(2,i)-1.0_r64
      uvtmp(2,i) = -(2.0_r64*uvs(1,i)-1.0_r64)
    end do
    call build_tensor_interpolation_r64(uvtmp, hdim, x0, p, MtlB)
    uvtmp = -uvtmp
    call build_tensor_interpolation_r64(uvtmp, hdim, x0, p, MtrB)

    allocate(xx1(p2), xx2(p2), sxk(3,p2), rts(3,p2), rps(3,p2), S(9,p2))
    do b = 1_8, p
      do a = 1_8, p
        m = a+(b-1_8)*p
        xx1(m) = x0(b)
        xx2(m) = x0(a)
      end do
    end do

    do ipan = 1_8, nchart
      do m = 1_8, p2
        call chart_eval_surface_r64(charts(ipan),ts1,ts2,xx1(m),xx2(m), &
             nfp,nmode,mn,rc,zs,sxk(:,m))
      end do

      do col = 1_8, p
        do row = 1_8, p
          do idim = 1_8, 3_8
            sumd = 0.0_r64
            do k = 1_8, p
              sumd = sumd+D(row,k)*sxk(idim,k+(col-1_8)*p)
            end do
            rts(idim,row+(col-1_8)*p) = sumd
          end do
        end do
      end do
      do row = 1_8, p
        do col = 1_8, p
          do idim = 1_8, 3_8
            sumd = 0.0_r64
            do k = 1_8, p
              sumd = sumd+D(col,k)*sxk(idim,row+(k-1_8)*p)
            end do
            rps(idim,row+(col-1_8)*p) = sumd
          end do
        end do
      end do

      d1 = 0.0_r64
      d2 = 0.0_r64
      do idim = 1_8, 3_8
        delta(idim) = sxk(idim,1)-sxk(idim,p2)
        d1 = d1+delta(idim)*delta(idim)
        delta(idim) = sxk(idim,p)-sxk(idim,p2-p+1_8)
        d2 = d2+delta(idim)*delta(idim)
      end do
      if (d1+1.0e-13_r64 > d2) then
        diagonal_a = .true.
      else
        diagonal_a = .false.
      end if

      do m = 1_8, p2
        S(1:3,m) = sxk(:,m)
        S(4:6,m) = rts(:,m)
        S(7:9,m) = rps(:,m)
      end do

      do side = 1_8, 2_8
        if (diagonal_a .and. side == 1_8) then
          Mt => MtlA
        else if (diagonal_a) then
          Mt => MtrA
        else if (side == 1_8) then
          Mt => MtlB
        else
          Mt => MtrB
        end if
        off = (2_8*(ipan-1_8)+side-1_8)*hdim
        do i = 1_8, hdim
          acc = 0.0_r64
          do m = 1_8, p2
            coeff = Mt(m,i)
            do idim = 1_8, 9_8
              acc(idim) = acc(idim)+S(idim,m)*coeff
            end do
          end do
          xt = acc(4:6)
          xp = acc(7:9)
          normal(1) = xp(2)*xt(3)-xp(3)*xt(2)
          normal(2) = xp(3)*xt(1)-xp(1)*xt(3)
          normal(3) = xp(1)*xt(2)-xp(2)*xt(1)
          speed = sqrt(sum(normal*normal))
          sx(:,off+i) = acc(1:3)
          snx(:,off+i) = normal/speed
          rts_out(:,off+i) = 2.0_r64*xt
          rps_out(:,off+i) = 2.0_r64*xp
          sw(off+i) = speed*wts(i)*4.0_r64
        end do
      end do
    end do
  end subroutine build_mesh_r64

end module stellarator_mesh_mod
