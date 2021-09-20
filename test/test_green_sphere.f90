      implicit real *8 (a-h,o-z)
      
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      integer, allocatable :: iptype(:),norders(:),ixyzs(:)

      real *8, allocatable :: srcvals_nc(:,:),srccoefs_nc(:,:)
      integer, allocatable :: iptype_nc(:),norders_nc(:), &
        ixyzs_nc(:)
      
      complex *16, allocatable :: uval_nc(:),dudnval_nc(:)
      complex *16, allocatable :: pot_s_nc(:),pot_d_nc(:)

      complex *16, allocatable :: uval(:),dudnval(:)
      complex *16, allocatable :: pot_s(:),pot_d(:)

      real *8 xyz0(3),xyz_out(3)

      complex *16 zk,ima,zpars(5)
      integer ppw
      integer, allocatable :: ipatch_id(:),ipatch_id_nc(:)
      real *8, allocatable :: uvs_targ(:,:),uvs_targ_nc(:,:)
      real *8, allocatable :: wts(:),wts_nc(:)
      complex *16 pottmp,zpig
      complex *16 hvals(0:100),hder(0:100)
      complex *16 jvals(0:100),jder(0:100),zfac_s,zfac_d
      data ima/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4


!     
!      npatches = npatches_rhab
 
      zk = 2.0d0*pi
      rlam = 2.0d0*pi/zk
      norder = 5
      npols = (norder+1)*(norder+2)/2
      a = 1.0d0
      b = 1.0d0
      c = 1.0d0
      ppw = 20

      xyz0(1) = 0.0d0
      xyz0(2) = 0.0d0
      xyz0(3) = 0.0d0

      xyz_out(1) = 7.1d0
      xyz_out(2) = 3.1d0
      xyz_out(3) = -2.6d0
      
      call get_polar_uvmem_nc(a,b,c,rlam,ppw,norder,npatches_nc)
      npts_nc = npatches_nc*npols
      allocate(srcvals_nc(12,npts_nc),srccoefs_nc(9,npts_nc))
      allocate(iptype_nc(npatches_nc),ixyzs_nc(npatches_nc+1))
      allocate(norders_nc(npatches_nc))
      allocate(wts_nc(npts_nc))
      allocate(ipatch_id_nc(npts_nc),uvs_targ_nc(2,npts_nc))

      call get_ellipsoid_geom_polar_nc(a,b,c,xyz0,rlam,ppw, &
        norder,npatches_nc,norders_nc,ixyzs_nc,iptype_nc,npts_nc, &
        srcvals_nc,srccoefs_nc)
      call get_qwts(npatches_nc,norders_nc,ixyzs_nc,iptype_nc,npts_nc, &
        srcvals_nc,wts_nc)
      call get_patch_id_uvs(npatches_nc,norders_nc,ixyzs_nc,iptype_nc, &
        npts_nc,ipatch_id_nc,uvs_targ_nc)
      

      call get_rectparapiped_mem(a,b,c,rlam,ppw,norder,npatches)
      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts),iptype(npatches))
      allocate(ixyzs(npatches+1),norders(npatches),wts(npts))
      allocate(ipatch_id(npts),uvs_targ(2,npts))

      call get_ellipsoid_geom(a,b,c,xyz0,rlam,ppw,norder,npatches, &
        norders,ixyzs,iptype,npts,srcvals,srccoefs)
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype, &
        npts,ipatch_id,uvs_targ)

      print *, npatches,npatches_nc

      allocate(uval(npts),dudnval(npts),pot_s(npts),pot_d(npts))
      allocate(uval_nc(npts_nc),dudnval_nc(npts_nc), &
        pot_s_nc(npts_nc),pot_d_nc(npts_nc))

      eps = 1.0d-8
      zpars(1) = zk

      do i=1,npts
        call h3d_slp(xyz_out,3,srcvals(1,i),0,dpars,1,zpars,0,ipars, &
          uval(i))
        call h3d_sprime(xyz_out,12,srcvals(1,i),0,dpars,1,zpars,0, &
          ipars,dudnval(i))
      enddo
      do i=1,npts_nc
        call h3d_slp(xyz_out,3,srcvals_nc(1,i),0,dpars,1,zpars,0,ipars, &
          uval_nc(i))
        call h3d_sprime(xyz_out,12,srcvals_nc(1,i),0,dpars,1,zpars,0, &
          ipars,dudnval_nc(i))
      enddo

!
!
!  evaluate single layer on conforming mesh
! 
      zpars(2) = 1.0d0
      zpars(3) = 0.0d0


      call lpcomp_helm_comb_dir(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,12,npts,srcvals,ipatch_id,uvs_targ,eps,zpars, &
        dudnval,pot_s)

      
!
! evaluate single layer non-conforming mesh
!
!
      call lpcomp_helm_comb_dir(npatches_nc,norders_nc,ixyzs_nc, &
        iptype_nc,npts_nc,srccoefs_nc,srcvals_nc,12,npts_nc,srcvals_nc, &
        ipatch_id_nc,uvs_targ_nc,eps,zpars,dudnval_nc,pot_s_nc)

!
!
!  evaluate double layer on conforming mesh
! 
      zpars(2) = 0.0d0
      zpars(3) = 1.0d0


      call lpcomp_helm_comb_dir(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,12,npts,srcvals,ipatch_id,uvs_targ,eps,zpars, &
        uval,pot_d)

!
! evaluate double layer non-conforming mesh
!
!
      call lpcomp_helm_comb_dir(npatches_nc,norders_nc,ixyzs_nc, &
        iptype_nc,npts_nc,srccoefs_nc,srcvals_nc,12,npts_nc,srcvals_nc, &
        ipatch_id_nc,uvs_targ_nc,eps,zpars,uval_nc,pot_d_nc)

      erra_nc = 0.0d0
      erra_c = 0.0d0
      do i=1,npts_nc
        pottmp = (pot_s_nc(i) - pot_d_nc(i))*2
        erra_nc = erra_nc + abs(pottmp-uval_nc(i))**2*wts_nc(i)
      enddo
      erra_nc = sqrt(erra_nc)


      do i=1,npts
        pottmp = (pot_s(i) - pot_d(i))*2
        erra_c = erra_c + abs(pottmp-uval(i))**2*wts(i)
      enddo
      erra_c = sqrt(erra_c)
      call prin2('error in greens identity conforming=*',erra_c,1)
      call prin2('error in greens identity non-conforming=*',erra_nc,1)





      stop
      end
!
!
!
!
!      
      subroutine get_polar_uvmem_nc(a,b,c,rlam,ppw,norder,npatches)
      implicit real *8 (a-h,o-z)
      integer ppw,norder,npatches

      done = 1.0d0
      pi = atan(done)*4

      nthet = ceiling(2*c*(ppw+0.0d0)/rlam/(norder+0.0d0))
      hthet = (pi+0.0d0)/(nthet+0.0d0)

      npatches = 0
      do ithet=1,nthet
        t0 = (ithet-1)*hthet
        t1 = (ithet)*hthet

        tuse = t0
        if(abs(t0-pi/2).ge.abs(t1-pi/2)) tuse = t1
        if((t0-pi/2)*(t1-pi/2).le.0) tuse = pi/2


        alpha = a*sin(tuse)
        beta = b*sin(tuse)
        hh = (alpha-beta)**2/(alpha+beta)**2 

        ellip_p = pi*(alpha + beta)* &
           (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
        nphi = ceiling(ellip_p*(ppw+0.0d0)/rlam/(norder+0.0d0))
        npatches = npatches + 2*nphi
      enddo
      

      return
      end
!
!
!
!

      subroutine get_ellipsoid_geom_polar_nc(a,b,c,xyz0,rlam,ppw, &
        norder,npatches,norders,ixyzs,iptype,npts,srcvals,srccoefs)
      implicit real *8 (a-h,o-z)

      real *8 v1(3),v2(3),v3(3),v4(3),xyz0(1:3)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      real *8, allocatable :: uvs(:,:),wts(:),umatr(:,:),vmatr(:,:)
      real *8 srcvals(12,npts),srccoefs(9,npts)
      integer iptype(npatches),ixyzs(npatches+1),norders(npatches)
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      integer ppw
      external xtri_ellipsoid_polar_eval

      npols = (norder+1)*(norder+2)/2
      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = (i-1)*npols + 1
        iptype(i) = 1
      enddo
      ixyzs(npatches+1) = npts+1

      allocate(triaskel(3,3,npatches))

      done = 1.0d0
      pi = atan(done)*4

      nthet = ceiling(2*c*(ppw+0.0d0)/rlam/(norder+0.0d0))
      hthet = (pi+0.0d0)/(nthet+0.0d0)

      nthet0 = 1
      istart = 1

      do ithet = 1,nthet
        n0 = 0

        umin = (ithet-1)*hthet
        umax = ithet*hthet
        vmin = 0
        vmax = 2*pi

        tuse = umin
        if(abs(umin-pi/2).ge.abs(umax-pi/2)) tuse = umax
        if((umin-pi/2)*(umax-pi/2).le.0) tuse = pi/2

        nthet0 = 1
        nover = 0
        alpha = a*sin(tuse)
        beta = b*sin(tuse)
        hh = (alpha-beta)**2/(alpha+beta)**2 

        ellip_p = pi*(alpha + beta)* &
           (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
        nphi = ceiling(ellip_p*(ppw+0.0d0)/rlam/(norder+0.0d0))
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,nthet0,nphi,nover, &
          npatches,n0,triaskel(1,1,istart))
        istart = istart + n0
      enddo

      p2(1) = a
      p2(2) = b
      p2(3) = c

      p3(1:3) = xyz0(1:3)
      
      ptr1 => triaskel(1,1,1)
      ptr2 => p2(1)
      ptr3 => p3(1)
      ptr4 => p4(1)

      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),wts(npols),umatr(npols,npols), &
        vmatr(npols,npols))
      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)
      call getgeominfo(npatches,xtri_ellipsoid_polar_eval,ptr1,ptr2,ptr3, &
        ptr4,npols,uvs,umatr,srcvals,srccoefs)


      return
      end
!
!
!
!
!
!


      subroutine xtri_ellipsoid_polar_eval(itri, u, v, xyz, dxyzduv, & 
          triainfo,p2, p3, p4)
      implicit real *8 (a-h,o-z)
      real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*),p2(3),p3(3)


      !
      ! project the triangle itri in triainfo onto the sphere
      !
      !    Input:
      ! itri - triangle number to map
      ! u,v - local uv coordinates on triangle itri
      ! triainfo - flat skeleton triangle info
      ! p2,p3,p4 - dummy parameters
      !
      !    Output:
      ! xyz - point on the sphere
      ! dxyzduv - first derivative information
      !
      !

      x0=triainfo(1,1,itri)
      y0=triainfo(2,1,itri)
      z0=triainfo(3,1,itri)

      x1=triainfo(1,2,itri)
      y1=triainfo(2,2,itri)
      z1=triainfo(3,2,itri)

      x2=triainfo(1,3,itri)
      y2=triainfo(2,3,itri)
      z2=triainfo(3,3,itri)

      !
      ! ... process the geometry, return the point location on the sphere
      ! and the derivatives with respect to u and v
      !
      x=x0+u*(x1-x0)+v*(x2-x0)
      y=y0+u*(y1-y0)+v*(y2-y0)
      z=z0+u*(z1-z0)+v*(z2-z0)

      dxdu = x1-x0
      dydu = y1-y0
      dzdu = z1-z0
    
      dxdv = x2-x0
      dydv = y2-y0
      dzdv = z2-z0

      !
      ! second derivatives are zero...
      !

      !
      ! project onto the sphere
      !
      xyz(1)=p2(1)*sin(x)*cos(y) + p3(1)
      xyz(2)=p2(2)*sin(x)*sin(y) + p3(2)
      xyz(3)=p2(3)*cos(x) + p3(3)


      ! du
      dxyzduv(1,1) = p2(1)*(cos(x)*cos(y)*dxdu - sin(x)*sin(y)*dydu)
      dxyzduv(2,1) = p2(2)*(cos(x)*sin(y)*dxdu + sin(x)*cos(y)*dydu)
      dxyzduv(3,1) = p2(3)*(-sin(x)*dxdu)

      ! dv
      dxyzduv(1,2) = p2(1)*(cos(x)*cos(y)*dxdv - sin(x)*sin(y)*dydv)
      dxyzduv(2,2) = p2(2)*(cos(x)*sin(y)*dxdv + sin(x)*cos(y)*dydv)
      dxyzduv(3,2) = p2(3)*(-sin(x)*dxdv)

      return
      end subroutine xtri_ellipsoid_polar_eval
!
!
!
!
      subroutine get_ellipsoid_geom(a,b,c,xyz0,rlam,ppw,norder, &
        npatches,norders,ixyzs,iptype,npts,srcvals,srccoefs)
      implicit real *8 (a-h,o-z)

      real *8 v1(3),v2(3),v3(3),v4(3),xyz0(1:3)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      real *8, allocatable :: uvs(:,:),wts(:),umatr(:,:),vmatr(:,:)
      real *8 srcvals(12,*),srccoefs(9,*)
      integer norders(npatches),iptype(npatches),ixyzs(npatches+1)
      integer ppw
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      external xtri_ellipsoid_eval



      allocate(triaskel(3,3,npatches))

      call get_rectparapiped(a,b,c,rlam,ppw,norder, &
        npatches,triaskel)
      
      p2(1) = a
      p2(2) = b
      p2(3) = c

      p3(1:3) = xyz0(1:3)
      
      ptr1 => triaskel(1,1,1)
      ptr2 => p2(1)
      ptr3 => p3(1)
      ptr4 => p4(1)

      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),wts(npols),umatr(npols,npols), &
        vmatr(npols,npols))
      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)
      npts = npols*npatches
      call getgeominfo(npatches,xtri_ellipsoid_eval,ptr1,ptr2,ptr3, &
        ptr4,npols,uvs,umatr,srcvals,srccoefs)
      do i=1,npatches
        ixyzs(i) = (i-1)*npols+1
        iptype(i) = 1
        norders(i) = norder
      enddo
      ixyzs(npatches+1) = npts+1
       

      return
      end
!
!
!
!
!
!

      subroutine get_rectparapiped(aa,bb,cc,rlam,ppw,norder,npatches, &
        triaskel)
      implicit real *8 (a-h,o-z)

      real *8 triaskel(3,3,npatches),vs(3,4)
      real *8 vcube(3,8),xnorm(3)
      integer ppw

      a = aa/sqrt(3.0d0)
      b = bb/sqrt(3.0d0)
      c = cc/sqrt(3.0d0)
      vcube(1,1) = -a
      vcube(2,1) = -b
      vcube(3,1) = -c

      vcube(1,2) = a
      vcube(2,2) = -b
      vcube(3,2) = -c

      vcube(1,3) = a
      vcube(2,3) = b
      vcube(3,3) = -c

      vcube(1,4) = -a
      vcube(2,4) = b
      vcube(3,4) = -c

      vcube(1,5) = -a
      vcube(2,5) = -b
      vcube(3,5) = c

      vcube(1,6) = a
      vcube(2,6) = -b
      vcube(3,6) = c

      vcube(1,7) = a
      vcube(2,7) = b
      vcube(3,7) = c

      vcube(1,8) = -a
      vcube(2,8) = b
      vcube(3,8) = c

      na = ceiling(2*aa*(ppw+0.0d0)/rlam/(norder+0.0d0))
      nb = ceiling(2*bb*(ppw+0.0d0)/rlam/(norder+0.0d0))
      nc = ceiling(2*cc*(ppw+0.0d0)/rlam/(norder+0.0d0))



!       z = -c face      
      vs(1:3,1) = vcube(1:3,1)
      vs(1:3,2) = vcube(1:3,4)
      vs(1:3,3) = vcube(1:3,3)
      vs(1:3,4) = vcube(1:3,2)
      ntri = 0
      call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nb,na, &
        npatches,triaskel(1,1,ntri+1))

      call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
      call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)

      ntri = ntri + 2*na*nb
      

!       z = c face      
      vs(1:3,1:4) = vcube(1:3,5:8)
      call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),na,nb, &
        npatches,triaskel(1,1,ntri+1))

      call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
      call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)

      ntri = ntri + 2*na*nb

!      y = -b face
!
      vs(1:3,1) = vcube(1:3,1)
      vs(1:3,2) = vcube(1:3,2)
      vs(1:3,3) = vcube(1:3,6)
      vs(1:3,4) = vcube(1:3,5)

      call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),na,nc, &
        npatches,triaskel(1,1,ntri+1))
      call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
      call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)


      ntri = ntri + 2*na*nc

!      y = b face
!
      vs(1:3,1) = vcube(1:3,4)
      vs(1:3,2) = vcube(1:3,8)
      vs(1:3,3) = vcube(1:3,7)
      vs(1:3,4) = vcube(1:3,3)

      call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nc,na, &
        npatches,triaskel(1,1,ntri+1))
      call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
      call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)

      ntri = ntri + 2*na*nc



!      x = -a face
!
      vs(1:3,1) = vcube(1:3,1)
      vs(1:3,2) = vcube(1:3,5)
      vs(1:3,3) = vcube(1:3,8)
      vs(1:3,4) = vcube(1:3,4)

      call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nc,nb, &
        npatches,triaskel(1,1,ntri+1))
      call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
      call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)


      ntri = ntri + 2*nb*nc

!      x = 2 face
!
      vs(1:3,1) = vcube(1:3,2)
      vs(1:3,2) = vcube(1:3,3)
      vs(1:3,3) = vcube(1:3,7)
      vs(1:3,4) = vcube(1:3,6)

      call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nb,nc, &
        npatches,triaskel(1,1,ntri+1))
      call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
      call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)

      ntri = ntri + 2*nb*nc


      return
      end


      subroutine get_norm_triaskel(tria,xnorm)
      implicit real *8 (a-h,o-z)
      real *8 tria(3,3),xnorm(3),xu(3),xv(3)
      

      xu(1:3) = tria(1:3,2) - tria(1:3,1)
      xv(1:3) = tria(1:3,3) - tria(1:3,1)

      xnorm(1:3) = 0
      call cross_prod3d(xu,xv,xnorm)


      return
      end
      


      subroutine get_rectparapiped_mem(a,b,c,rlam,ppw,norder,npatches)
      implicit real *8 (a-h,o-z)
      integer na,nb,nc,ppw

      na = ceiling(2*a*(ppw+0.0d0)/rlam/(norder+0.0d0))
      nb = ceiling(2*b*(ppw+0.0d0)/rlam/(norder+0.0d0))
      nc = ceiling(2*c*(ppw+0.0d0)/rlam/(norder+0.0d0))

      npatches = 4*(na*nb + nb*nc + nc*na)

      return
      end subroutine get_rectparapiped_mem



      subroutine xtri_rectmesh_3d(v1,v2,v3,v4,nu,nv,npatches,triaskel)
      implicit real *8 (a-h,o-z)
      real *8 triaskel(3,3,npatches),v1(3),v2(3),v3(3),v4(3)
      real *8 vl(3),vr(3),vb(3),vt(3)
      real *8 uvw1(3),uvw2(3),uvw3(3),uvw4(3)

      vl(1:3) = v4(1:3)-v1(1:3)
      vr(1:3) = v3(1:3)-v2(1:3)



      ntri = 0

      do i=1,nv
        uvw1(1:3) = v1(1:3) + (i-1)*vl(1:3)/(nv+ 0.0d0)
        uvw4(1:3) = uvw1(1:3) + vl(1:3)/(nv+0.0d0)

        vb(1:3) =  v2(1:3) + (i-1)*vr(1:3)/(nv+0.0d0) - uvw1(1:3)
        vt(1:3) = uvw1(1:3) + vb(1:3) + vr(1:3)/(nv+0.0d0) - uvw4(1:3)
        
        uvw2(1:3) = uvw1(1:3) + vb(1:3)/(nu+0.0d0)
        uvw3(1:3) = uvw4(1:3) + vt(1:3)/(nu+0.0d0)

        do j=1,nu
          call xtri_rectmesh0_3d(uvw1,uvw2,uvw3,uvw4, &
            triaskel(1,1,ntri+1))
          ntri = ntri + 2
          uvw1(1:3) = uvw2(1:3)
          uvw2(1:3) = uvw2(1:3) + vb(1:3)/(nu+0.0d0)
          uvw4(1:3) = uvw3(1:3)
          uvw3(1:3) = uvw3(1:3) + vt(1:3)/(nu+0.0d0)
        enddo
      enddo

      end subroutine xtri_rectmesh_3d







      subroutine xtri_rectmesh0_3d(v1,v2,v3,v4,triaskel)
      implicit real *8 (a-h,o-z)
      real *8 v1(3),v2(3),v3(3),v4(3),triaskel(3,3,2)

      do i=1,3
        triaskel(i,1,1) = v1(i)
        triaskel(i,2,1) = v2(i)
        triaskel(i,3,1) = v4(i)
        triaskel(i,1,2) = v3(i)
        triaskel(i,2,2) = v4(i)
        triaskel(i,3,2) = v2(i)
      enddo

      return
      end subroutine xtri_rectmesh0_3d



      subroutine xtri_vtk_flat(iw, ntri, xtri1s, title)
      implicit real *8 (a-h,o-z)
      real *8 :: xtri1s(3,3,ntri)
      character(len=*) :: title

      character(len=1024) :: filename, dataname, valsname, imgname
      character(len=1024) :: trisname, vecsname, centname
      character(len=12) :: fmt, fmt3, fmt4
      character(len=25) :: fmt2

!
!  This routine plots a sequence of FLAT TRIANGLES with surface
! color vals.
!
! Input:
!   iw - plot number, controls the filenames
!   ntri - number of flat triangles
!   xtri1s - full triangle information
!
! Output:
!   files which can be executed in matlab to plot the surface
!
!

      if (iw .lt. 10) then
        fmt = "(A4,I1,A4)"
        fmt3 = "(A8,I1,A4)"
        fmt4 = "(A5,I1,A4)"
      elseif (iw .lt. 100) then
        fmt = "(A4,I2,A4)"
        fmt3 = "(A8,I2,A4)"
        fmt4 = "(A5,I2,A4)"
      elseif (iw .lt. 1000) then
        fmt = "(A4,I3,A4)"
        fmt3 = "(A8,I3,A4)"
        fmt4 = "(A5,I3,A4)"
      elseif (iw .lt. 10000) then
        fmt = "(A4,I4,A4)"
        fmt3 = "(A8,I4,A4)"
        fmt4 = "(A5,I4,A4)"
      end if

      write(filename, fmt) 'plot', iw, '.vtk'

      !
      ! write the vtk plotting script
      !
      iunit1 = 877
      open(unit = iunit1, file=trim(filename), status='replace')
    
      write(iunit1,'(a)') "# vtk DataFile Version 3.0"
      write(iunit1,'(a)') "vtk output"
      write(iunit1,'(a)') "ASCII"
      !write(iunit1,'(a)') "DATASET POLYDATA"
      write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
      write(iunit1,'(a,i8,a)') "POINTS ", ntri*3, " float"

      fmt2 = "(E11.5,2X,E11.5,2X,E11.5)"
      do i = 1,ntri
        write(iunit1,fmt2) xtri1s(1,1,i), xtri1s(2,1,i), xtri1s(3,1,i)
        write(iunit1,fmt2) xtri1s(1,2,i), xtri1s(2,2,i), xtri1s(3,2,i)
        write(iunit1,fmt2) xtri1s(1,3,i), xtri1s(2,3,i), xtri1s(3,3,i)
      end do


      write(iunit1,'(a,i8,i8)') "CELLS ", ntri, ntri*4

      do i = 1,ntri
        i1 = 3*(i-1) + 1
        write(iunit1,'(a,i8,i8,i8)') "3 ", i1-1, i1, i1+1
      end do

      write(iunit1,'(a,i8)') "CELL_TYPES ", ntri
      do i = 1,ntri
        write(iunit1,'(a)') "5"
      end do

      write(iunit1,'(a,i8)') "POINT_DATA ", ntri*3
      write(iunit1,'(a)') "SCALARS scalars float 1"
      write(iunit1,'(a)') "LOOKUP_TABLE default"
      do i = 1,ntri
        do j = 1,3
          write(iunit1,'(E11.5)') xtri1s(3,j,i)
        end do
      end do



      write(iunit1,'(a)') ""
      write(iunit1,'(a)') ""
      write(iunit1,'(a,i8)') "CELL_DATA ", ntri
      write(iunit1,'(a)') "SCALARS scalars float 1"
      write(iunit1,'(a)') "LOOKUP_TABLE default"
      do i = 1,ntri
        write(iunit1,'(E13.5)') (xtri1s(3,1,i) + &
        xtri1s(3,2,i) + xtri1s(3,3,i))/3
      end do
    
      close(iunit1)

      return
      end subroutine xtri_vtk_flat





      subroutine xtri_ellipsoid_eval(itri, u, v, xyz, dxyzduv, & 
          triainfo,p2, p3, p4)
      implicit real *8 (a-h,o-z)
      real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*),p2(3),p3(3)


      !
      ! project the triangle itri in triainfo onto the sphere
      !
      !    Input:
      ! itri - triangle number to map
      ! u,v - local uv coordinates on triangle itri
      ! triainfo - flat skeleton triangle info
      ! p2,p3,p4 - dummy parameters
      !
      !    Output:
      ! xyz - point on the sphere
      ! dxyzduv - first derivative information
      !
      !

      x0=triainfo(1,1,itri)
      y0=triainfo(2,1,itri)
      z0=triainfo(3,1,itri)

      x1=triainfo(1,2,itri)
      y1=triainfo(2,2,itri)
      z1=triainfo(3,2,itri)

      x2=triainfo(1,3,itri)
      y2=triainfo(2,3,itri)
      z2=triainfo(3,3,itri)

      !
      ! ... process the geometry, return the point location on the sphere
      ! and the derivatives with respect to u and v
      !
      x=x0+u*(x1-x0)+v*(x2-x0)
      y=y0+u*(y1-y0)+v*(y2-y0)
      z=z0+u*(z1-z0)+v*(z2-z0)

      dxdu = x1-x0
      dydu = y1-y0
      dzdu = z1-z0
    
      dxdv = x2-x0
      dydv = y2-y0
      dzdv = z2-z0

      !
      ! second derivatives are zero...
      !

      !
      ! project onto the sphere
      !
      r=sqrt(x**2+y**2+z**2)
      xyz(1)=p2(1)*x/r + p3(1)
      xyz(2)=p2(2)*y/r + p3(2)
      xyz(3)=p2(3)*z/r + p3(3)

      a = x0*(x1-x0) + y0*(y1-y0) + z0*(z1-z0)
      b = (x1-x0)*(x2-x0) + (y1-y0)*(y2-y0) + (z1-z0)*(z2-z0)
      c = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0)

      drdu = (a + v*b + u*c)/r
      drdu2 = (r*c - r*drdu*drdu)/r/r

      e = x0*(x2-x0) + y0*(y2-y0) + z0*(z2-z0)
      f = b
      g = (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0)

      drdv = (e + u*f + v*g)/r
      drdv2 = (r*g - r*drdv*drdv)/r/r

      drduv = (r*b - r*drdu*drdv)/r/r

      ! du
      dxyzduv(1,1) = p2(1)*(r*dxdu-x*drdu)/r/r
      dxyzduv(2,1) = p2(2)*(r*dydu-y*drdu)/r/r
      dxyzduv(3,1) = p2(3)*(r*dzdu-z*drdu)/r/r

      ! dv
      dxyzduv(1,2) = p2(1)*(r*dxdv-x*drdv)/r/r
      dxyzduv(2,2) = p2(2)*(r*dydv-y*drdv)/r/r
      dxyzduv(3,2) = p2(3)*(r*dzdv-z*drdv)/r/r

      return
      end subroutine xtri_ellipsoid_eval








   



      subroutine l3getsph(nmax,mm,nn,ndx,xyzs,ynms,npts,ynm)
      implicit real *8 (a-h,o-z)
      real *8 :: xyzs(ndx,npts)
      complex *16 ynms(npts),ima
      real *8 rat1(10000),rat2(10000)
      real *8 ynm(0:nmax,0:nmax)
      data ima/(0.0d0,1.0d0)/
  
      call ylgndrini(nmax, rat1, rat2)
  
      do i=1,npts
        x=xyzs(1,i)
        y=xyzs(2,i)
        z=xyzs(3,i)
        r=sqrt(x**2+y**2+z**2)
        call cart2polar(xyzs(1,i),r,theta,phi)
        ctheta = cos(theta)
        call ylgndrf(nmax, ctheta, ynm, rat1, rat2)
        ynms(i) = ynm(nn,abs(mm))*exp(ima*mm*phi)        
      enddo
       
      return
      end



