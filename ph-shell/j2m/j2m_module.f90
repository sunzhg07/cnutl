MODULE constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
  ! min and max isospin projection
  INTEGER, PUBLIC :: itzmin, itzmax
  ! min and max total two-body angular momentum in lab frame
  INTEGER, PUBLIC :: j_lab_min, j_lab_max
  ! min and max total two-body angular momentum in Rel-CoM frame
  INTEGER, PUBLIC :: jmin, jmax, n_nucleus, z_nucleus
  INTEGER, PUBLIC :: occ_protons, occ_neutrons
  real*8, public :: com_beta,com_switch, j2_beta, eom_energy
  REAL*8, public :: hcom_value,occnum3_cut,e0,o0,total_energy, eccsd
  INTEGER, PUBLIC :: mass_nucleus, tot_orbs, below_ef, above_ef,nlmax_cut, & 
       e3max_cut, l3max_cut, j3_cut, j2_max
  REAL(DP), PUBLIC ::  oscl, hbar_omega, efermi
  REAL(DP) , PARAMETER, PUBLIC :: p_mass = 938.926_dp
  REAL(DP) , PARAMETER, PUBLIC :: g_A = 1.29d0
  REAL(DP) , PARAMETER, PUBLIC :: fpi = 92.4d0
  REAL(DP), PARAMETER, PUBLIC :: theta_rot = 0.0_dp! 0.125_dp
  REAL(DP), PARAMETER, PUBLIC :: hbarc = 197.326968_dp
  REAL(DP), PARAMETER, PUBLIC :: hb2ip = hbarc*hbarc/p_mass
  REAL(DP), PUBLIC, PARAMETER :: pi = 3.141592741012573_dp
  REAL(DP), PUBLIC, PARAMETER :: pi_2 = 1.570796370506287_dp
  REAL(DP), PUBLIC, PARAMETER :: pi_4 = 0.7853981852531433_dp
  LOGICAL :: switch_density, switch_triples
  CHARACTER(LEN=100) :: operator_type, gs_type
  CHARACTER(LEN=40) :: clem_j
  integer, public :: stime=0 ! whether show the time of running
  
END MODULE constants



!
!            
!     This module contains the angular momentun functions
!     and transformation coefficients when going from 
!     lab system  <--> cm system
!

MODULE ang_mom_functions
  integer, private, parameter:: maxjj=100
  REAL*8, PRIVATE :: f_mb(maxjj),g_mb(maxjj),w_mb(maxjj)
  INTEGER, PRIVATE :: kh(4*maxjj)
  REAL*8, PARAMETER, PRIVATE :: pi=3.141592654
  REAL*8, PRIVATE :: q(maxjj,maxjj), cn(0:maxjj+1,0:maxjj+1)
 
CONTAINS
  !
  !     factorials for 3j,6j and 9j symbols            
  !     for moshinsky trans brackets and for           
  !     vector brackets                                
  !
  SUBROUTINE commons_to_angmom
    IMPLICIT NONE
    INTEGER :: l, k, i, j
    REAL*8 :: a , sq_pi, fj, tfj, fk
    !    3j, 6j and 9j symbols
    kh=1
    kh(2*maxjj) =0
    DO l=1,maxjj
       q(l,1)=1.0d0
       q(l,l)=1.0d0
       kh(l+l+2*maxjj)=0
    ENDDO
    DO l=2,maxjj-1
       DO k=2,l
          q(l+1,k)=q(l,k-1)+q(l,k)
       ENDDO
    ENDDO
    !    Moshinsky brackets
    f_mb(1)=0.
    g_mb(1)=LOG(0.5D0)
    w_mb(1)=0.
    DO i=2,maxjj
       a=i-1
       f_mb(i)=f_mb(i-1)+LOG(a)
       g_mb(i)=g_mb(i-1)+LOG(a+0.5D0)
       w_mb(i)=LOG(a+a+1.)
    ENDDO
    !    spherical harmonics
    cn=0.
    sq_pi=1./SQRT(2.*pi)
    DO j=0,maxjj+1
       cn(0,j)=SQRT(0.5*(2.*j+1.))
    ENDDO
    DO j=1,maxjj+1
       tfj=2.*j
       cn(j,j)=cn(j-1,j-1)*SQRT((tfj+1.)/tfj)
    ENDDO
    DO j=0,maxjj+1
       fj=FLOAT(j)
       DO k=1,j-1
          fk=FLOAT(k)
          cn(k,j)=cn(k-1,j)*SQRT((fj+fk)*(fj-fk+1.))*0.5/fk
       ENDDO
    ENDDO
    cn=cn*sq_pi

  END SUBROUTINE commons_to_angmom
  
  !
  !     calculates CG-symbols            
  !
  REAL*8 FUNCTION cgc(j_a,j_b,j_c,m_a,m_b,m_c)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j_a,j_b,j_c,m_a,m_b,m_c
    integer :: iph 

    cgc = tjs(j_a,j_b,j_c,m_a,m_b,-m_c) * iph( (j_a-j_b+m_c)/2) * sqrt(j_c+1.) 
  end FUNCTION cgc

  !
  !     calculates 3j-symbols           
  !
  REAL*8 FUNCTION tjs(j_a,j_b,j_c,m_a,m_b,m_c)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j_a,j_b,j_c,m_a,m_b,m_c
    INTEGER :: ja, jb, jc, mb, ma, mc, la, lb, lc, lt, ld, ja2, jb2, &
         jc2, i, k0, k1, k, ip
    REAL*8 :: x, fn, p
    logical :: triag 
    
    
    tjs=0.
    IF ( triag( j_a, j_b, j_c) ) return 
    IF ( m_a + m_b + m_c /= 0 ) return

    ja=(j_a+m_a)/2+1
    ma=(j_a-m_a)/2+1
    jb=(j_b+m_b)/2+1
    mb=(j_b-m_b)/2+1
    jc=(j_c+m_c)/2+1
    mc=(j_c-m_c)/2+1
    la=(j_b+j_c-j_a)/2+1
    lb=(j_c+j_a-j_b)/2+1
    lc=(j_a+j_b-j_c)/2+1
    lt=(j_a+j_b+j_c)/2+1
    ld=MIN(ja,jb,jc,ma,mb,mc,la,lb,lc)
    IF(((m_a+m_b+m_c) <= 0).AND.(ld > 0)) THEN
       ja2=j_a+m_a
       jb2=j_b+m_b
       jc2=j_c+m_c
       i=ja2+jb2+jc2-ja2/2*2-jb2/2*2-jc2/2*2
       IF(i == 0) then
          fn=q(ja+ma-1,lc)*q(jb+mb-1,lc)/(q(lt,jc+mc-1)*q(lt+1,2) &
               *q(ja+ma-1,ja)*q(jb+mb-1,jb)*q(jc+mc-1,jc))
          k0=MAX(0,lc-ja,lc-mb)+1
          k1=MIN(lc,ma,jb)
          x=0.
          DO k=k0,k1
             x=-x-q(lc,k)*q(lb,ma-k+1)*q(la,jb-k+1)
          ENDDO
          ip=k1+lb+jc
          p=1-2*(ip-ip/2*2)
          tjs=p*x*SQRT(fn)
       ENDIF
    ENDIF
  END FUNCTION tjs

   !
  !     calculates 6j-symbols           
  !
  REAL*8 FUNCTION sjs(j_a,j_b,j_c,l_a,l_b,l_c)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j_a,j_b,j_c,l_a,l_b,l_c
    INTEGER :: ja,jb,jc,la,lb,lc,i,mt,ma,mb,mc,na,nb,nc,ka,&
         kb,kc,l,l0,l1, ihlp, idum
    REAL*8 :: x, fs, fss
    LOGICAL :: triag
    
    ihlp=2*maxjj-1

    sjs=0.0d0
    
    IF ( triag( j_a, j_b, j_c) ) return 
    IF ( triag( j_c, l_a, l_b) ) return 
    IF ( triag( j_a, l_b, l_c) ) return 
    IF ( triag( j_b, l_a, l_c) ) return 
    
    ja=j_a + 1
    jb=j_b + 1
    jc=j_c + 1
    la=l_a + 1
    lb=l_b + 1
    lc=l_c + 1
    i=kh(ja+jb-jc+ihlp)+kh(jb+jc-ja+ihlp)+kh(jc+ja-jb+ihlp)+kh(ja+lb-lc+ihlp) &
         +kh(lb+lc-ja+ihlp)+kh(lc+ja-lb+ihlp)+kh(la+jb-lc+ihlp)+kh(jb+lc-la+ihlp) &
         +kh(lc+la-jb+ihlp)+kh(la+lb-jc+ihlp)+kh(lb+jc-la+ihlp)+kh(jc+la-lb+ihlp)
    IF(i <= 0) THEN
       mt=(j_a+j_b+j_c)/2 + 2
       ma=(j_a+l_b+l_c)/2+ 2
       mb=(l_a+j_b+l_c)/2+ 2
       mc=(l_a+l_b+j_c)/2+ 2
       na=mt-ja
       nb=mt-jb
       nc=mt-jc
       ka=ma-lc
       kb=mb-lc
       kc=mc-jc

       idum=max(mt,ja+1,nc,ma,ka,mb,la+1,kb,mc,kc)

       if(idum.gt.maxjj) then
          write(6,*) 'increase maxjj in MODULE ang_mom_functions from', maxjj, 'to', idum
          stop
       end if

       fss=q(mt,ja+1)*q(ja,nc)/(q(ma,ja+1)*q(ja,ka)*q(mb,la+1)* &
            q(la,kb)*q(mc,la+1)*q(la,kc))
       fs=SQRT(fss)/(l_a + 1.)
       l0=MAX(mt,ma,mb,mc)+1
       l1=MIN(ma+na,mb+nb,mc+nc)
       x=0.
       DO l=l0,l1
          x=-x+q(l-1,mt)*q(na,l-ma)*q(nb,l-mb)*q(nc,l-mc)
       ENDDO
       sjs=-(1+2*(l1/2*2-l1))*fs*x
    ENDIF

  END FUNCTION sjs
  !
  !     calculates ninej-symbols
  !      
  REAL*8 FUNCTION snj (ia,ib,ie,ic,id,if,ig,ih,it)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ia,ib,ie,ic,id,if,ig,ih,it
    INTEGER :: ja,jb,je,jc,jd,jf,jg,jh,jt,i,la,ld,ma,mc,na,nb,le,lf,&
         lg,me,mf,mg,ne,nf,ng,lx,mx,nx,jsi,jsf, js,is,lb, lc, &
         mb, ly, my,ny,l,l0,m0,n0,l1,m1,n1,m,n,ihx, ihlp
    REAL*8 :: x, fn, fd, ps, fs, u, y, z, ud, p

    ihlp=2*maxjj-1

    snj=0.
    ja=ia+1
    jb=ib+1
    jc=ic+1
    jd=id+1
    je=ie+1
    jf=IF+1
    jg=ig+1
    jh=ih+1
    jt=it+1
    i=kh(ja+jb-je+ihlp)+kh(jb+je-ja+ihlp)+kh(je+ja-jb+ihlp)+kh(jc+jd-jf+ihlp) &
         +kh(jd+jf-jc+ihlp)+kh(jf+jc-jd+ihlp)+kh(jg+jh-jt+ihlp)+kh(jh+jt-jg+ihlp) &
         +kh(jt+jg-jh+ihlp)+kh(ja+jc-jg+ihlp)+kh(jc+jg-ja+ihlp)+kh(jg+ja-jc+ihlp) &
         +kh(jb+jd-jh+ihlp)+kh(jd+jh-jb+ihlp)+kh(jh+jb-jd+ihlp)+kh(je+jf-jt+ihlp) &
         +kh(jf+jt-je+ihlp)+kh(jt+je-jf+ihlp)
    IF(i <= 0) THEN
       la=(ie+IF+it)/2+2
       ld=(ig+ih+it)/2+2
       ma=(ia+ic+ig)/2+2
       mc=(IF+ic+id)/2+2
       na=(ib+id+ih)/2+2
       nb=(ib+ie+ia)/2+2
       le=(ie+IF-it)/2+1
       lf=(IF+it-ie)/2+1
       lg=(it+ie-IF)/2+1
       me=(ia+ic-ig)/2+1
       mf=(ic+ig-ia)/2+1
       mg=(ig+ia-ic)/2+1
       ne=(ib+id-ih)/2+1
       nf=(id+ih-ib)/2+1
       ng=(ih+ib-id)/2+1
       lx=(it+ig-ih)/2+1
       mx=(ic+id-IF)/2+1
       nx=(ib+ie-ia)/2+1
       fn=q(la,jt+1)*q(jt,lg)*q(ma,jc+1)*q(jc,mf)*q(na,jb+1)*q(jb,ne)
       fd=q(ld,jt+1)*q(jt,lx)*q(mc,jc+1)*q(jc,mx)*q(nb,jb+1)*q(jb,nx)
       jsi=MAX(ABS(je-jh),ABS(jg-jf),ABS(ja-jd))+1
       jsf=MIN(je+jh,jg+jf,ja+jd)-1
       ps=-1-2*(jsi/2*2-jsi)
       fs=ps*SQRT(fn/fd)/FLOAT((ig+1)*(ie+1))
       u=0.
       DO js=jsi,jsf,2
          is=js-1
          lb=(ie+ih+is)/2+2
          lc=(ig+IF+is)/2+2
          mb=(ia+id+is)/2+2
          ly=(ie+ih-is)/2+1
          my=(ig+IF-is)/2+1
          ny=(ia-id+is)/2+1
          ud=q(lb,je+1)*q(je,ly)*q(lc,jg+1)*q(jg,my)*q(mb,js+1)*q(js,ny)
          l0=MAX(la,lb,lc,ld)+1
          m0=MAX(ma,mb,mc,lc)+1
          n0=MAX(na,nb,mb,lb)+1
          l1=MIN(le+ld,lf+lb,lg+lc)
          m1=MIN(me+lc,mf+mb,mg+mc)
          n1=MIN(ne+lb,nf+nb,ng+mb)
          x=0.
          DO l=l0,l1
             x=-x-q(l-1,la)*q(le,l-ld)*q(lf,l-lb)*q(lg,l-lc)
          ENDDO
          y=0.
          DO m=m0,m1
             y=-y-q(m-1,ma)*q(me,m-lc)*q(mf,m-mb)*q(mg,m-mc)
          ENDDO
          z=0.
          DO n=n0,n1
             z=-z-q(n-1,na)*q(ne,n-lb)*q(nf,n-nb)*q(ng,n-mb)
          ENDDO
          ihx=l1+m1+n1
          p=1+2*(ihx/2*2-ihx)
          u=u+p*x*y*z/ud
       ENDDO
       snj=u*fs
    ENDIF

  END FUNCTION snj

  !
  !     This routine calculates the moshinsky vector bracket      
  !     Note that D=mass1/mass2                                   
  !     Ref  m.sotona and m.gmitro  comp.phys.comm 3(1972)53      
  !
  DOUBLE PRECISION FUNCTION gmosh &
       (n,l,nc,lc,n1,l1,n2,l2,lr,d)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n,l,nc,lc,n1,l1,n2,l2,lr
    DOUBLE PRECISION, INTENT(IN) :: d
    INTEGER :: ip,ixf,ix, iyi, iyf, j1f,j2,k1i,k1f,m1f,iy,m2f,k2, &
         m2,m2i,m1,j1,k2f,k2i,k1
    DOUBLE PRECISION :: dl, d1l, bb, ba, anorm, y, p, bc, cfac, bm , &
         sm, s, sxy, bxy

    gmosh=0.
    IF(n+n+nc+nc+l+lc-n1-n1-n2-n2-l1-l2 /= 0 ) RETURN
    IF(l+lc-lr < 0 ) RETURN
    IF(l1+l2-lr < 0 ) RETURN
    IF(ABS(l-lc)-lr > 0 ) RETURN
    IF(ABS(l1-l2)-lr > 0 ) RETURN
    DL=LOG(D)
    D1L=LOG(D+1.)
    bb=f_mb(n1+1)+f_mb(n2+1)+f_mb(n+1)-f_mb(nc+1)+ &
         g_mb(n1+l1+1)+g_mb(n2+l2+1) &
         -g_mb(n+l+1)-g_mb(nc+lc+1)
    ba=w_mb(l1+1)+w_mb(l2+1)+w_mb(lc+1)+w_mb(l+1)+ &
         f_mb(l1+l2-lr+1)+f_mb(l+lc+lr+2) &
         +f_mb(l+lc-lr+1)+f_mb(lc+lr-l+1)+ &
         f_mb(lr+l-lc+1)-f_mb(l1+l2+lr+2) &
         -f_mb(l1+lr-l2+1)-f_mb(l2+lr-l1+1)-DBLE(l)*d1l
    ip=lr+n+n1+n2
    p=1+2*(ip/2*2-ip)
    anorm=p*EXP(0.5D0*(bb+ba))
    y=0.
    j1f=l+1
    DO j1=1,j1f
       j2=l+2-j1
       k1i=ABS(l1-j1+1)+1
       k1f=l1+j1
       DO k1=k1i,k1f,2
          m1f=n1-(j1+k1-l1)/2+2
          IF(m1f-1 < 0 )  CYCLE
          k2i=MAX(ABS(l2-j2+1),ABS(lc-k1+1))+1
          k2f=MIN(l2+j2,lc+k1)
          IF(k2i-k2f > 0 ) CYCLE
          DO k2=k2i,k2f,2
             m2f=n2-(j2+k2-l2)/2+2
             IF(m2f-1 < 0 )  CYCLE
             ip=j2-1+(l1+k1+j1+l2+k2+j2)/2
             p=1+2*(ip/2*2-ip)
             bc=0.5D0*(DBLE(k1+j2-2)*dl-DBLE(k1+k2-2)*d1l) &
                  +f_mb(k1+l1-j1+1)+f_mb(k1+k2-lc-1)+ &
                  f_mb(k2+l2-j2+1)-f_mb(k1+l1+j1)-f_mb(k1+k2+lc)- &
                  f_mb(k2+l2+j2)+w_mb(k1)+w_mb(k2)+f_mb((k1+l1+j1)/2)+ &
                  f_mb((k1+k2+lc)/2)+f_mb((k2+l2+j2)/2)- &
                  f_mb((k1+l1-j1)/2+1)-f_mb((l1+j1-k1)/2+1)- &
                  f_mb((j1+k1-l1)/2)-f_mb((k1+k2-lc)/2)- &
                  f_mb((k2+lc-k1)/2+1)-f_mb((lc+k1-k2)/2+1) &
                  -f_mb((k2+l2-j2)/2+1)-f_mb((l2+j2-k2)/2+1)- &
                  f_mb((j2+k2-l2)/2)
             cfac=p*EXP(bc)
             sxy=0.
             ixf=MIN(k1+k1,k1+k2-lc)-1
             DO ix=1,ixf
                iyi=MAX(1,ix+j1+l2-k1-lr)
                iyf=MIN(l2+l2+1,l1+l2-lr+1,l2+lc+ix-k1-j2+2)
                IF(iyi-iyf > 0 ) CYCLE
                DO iy=iyi,iyf
                   ip=ix+iy
                   p=1+2*(ip/2*2-ip)
                   bxy=f_mb(k1+k1-ix)+f_mb(l2+l2-iy+2)+ &
                        f_mb(k2+lc-k1+ix)+f_mb(l1+lr-l2+iy) &
                        -f_mb(ix)-f_mb(iy)-f_mb(k1+k2-lc-ix)- &
                        f_mb(l1+l2-lr-iy+2)-f_mb(k1-l2+lr-j1+iy-ix+1)- &
                        f_mb(l2-k1+lc-j2+ix-iy+3)
                   sxy=sxy+p*EXP(bxy)
                ENDDO
             ENDDO
             s=cfac*sxy
             sm=0.
             DO m1=1,m1f
                m2i=MAX(1,nc-m1-(k1+k2-lc)/2+3)
                IF(m2i-m2f > 0 ) CYCLE
                DO m2=m2i,m2f
                   ip=m1+m2
                   p=1+2*(ip/2*2-ip)
                   bm=DBLE(m1-1)*DL-DBLE(m1+m2-2)*d1l+g_mb(1) &
                        +g_mb(m1+m2+(k1+k2+lc)/2-2)-g_mb(k1+m1-1)- &
                        g_mb(k2+m2-1)+f_mb(m1+m2+(k1+k2-lc)/2-2)- &
                        f_mb(m1)-f_mb(m2)-f_mb(n1-m1-(j1+k1-l1)/2+3)- &
                        f_mb(n2-m2-(j2+k2-l2)/2+3) &
                        -f_mb(m1+m2-nc+(k1+k2-lc)/2-2)
                   sm=sm+p*EXP(bm)
                ENDDO
             ENDDO
             y=y+s*sm
          ENDDO
       ENDDO
    ENDDO
    gmosh=anorm*y

  END FUNCTION  gmosh

  !     calculates the vector bracket                           
  !     <ak,l,akk,ll,lam | ak1,l1,ak2,l2,lam>                   
  !     Ref. Kung et al, PRC 19 (1979) 1063                     
  !     Momenta in units of MeV, vector bracket in Units of MeV^-4

  DOUBLE PRECISION FUNCTION vector_trcoefficients &
       (p_mass,ak,akk,ak1,ak2,l,ll,lam,l1,it1,l2,it2)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)  :: p_mass(0:3)
    DOUBLE PRECISION, INTENT(IN) :: ak,akk,ak1,ak2
    INTEGER , INTENT(IN) :: l,ll,lam,l1,it1,l2,it2
    INTEGER :: ic1, ic2, m, mm, mu, mmu, mi, mp, ix, ixx, md
    DOUBLE PRECISION :: sak, dak, tmass, xm, xm1, xm2, xm3, x, xx, &
         xa, xb, xc, aiii, xd, psa, sl2, sll, sggn,&
         sumin, sign, sl, sl1, dr, dr1, dr2

    vector_trcoefficients=0.
    sak=ak1+ak2
    dak=DABS(ak1-ak2)
    IF (akk > sak) RETURN
    IF (akk < dak) RETURN
    tmass=p_mass(it1)+p_mass(it2)
    xm=p_mass(it1)/tmass
    xm1=(p_mass(it2)/tmass)**2
    xm2=xm**2
    xm3=p_mass(it1)*p_mass(it2)/(tmass**2)
    x=(ak1*ak1-ak*ak-akk*akk*xm2)/(2.D0*xm*ak*akk)
    xx=x*x
    IF (xx > 1.) RETURN
    ic1=l1+l2+l+ll
    ic2=ic1-2*(ic1/2)
    IF (ic2 /= 0) RETURN
    aiii=0.
    xa=(ak1*ak1+ak*ak-akk*akk*xm2)/(2.d0*ak*ak1)
    xb=(ak1*ak1+akk*akk*xm2-ak*ak)/(2.D0*xm*ak1*akk)
    xc=(ak1*ak1*xm1+ak2*ak2*xm2-ak*ak)/(2.d0*ak1*ak2*xm3)
    md=0
    xd=1.
    sl1=spherical_harmonics(md,l1,xd)
    IF (sl1 == 0.) RETURN
    DO m=-l,l
       sl=spherical_harmonics(m,l,xa)
       IF (sl == 0.)  CYCLE
       DO mm=-ll,ll
          mu=m+mm
          IF (IABS(mu) > l2) CYCLE
          mmu=-mu
          dr1=tjs(l+l,ll+ll,lam+lam,m+m,mm+mm,mmu+mmu)
          dr2=tjs(l1+l1,l2+l2,lam+lam,md+md,mu+mu,mmu+mmu)
          dr=dr1*dr2
          sign=1.
          mi=m-l-l1+l2+ll
          mp=mi-2*(mi/2)
          IF(mp /= 0) sign=-1.0d0
          sll= spherical_harmonics(mm,ll,xb)
          IF (sll == 0) CYCLE
          sl2= spherical_harmonics(mu,l2,xc)
          IF (sl2 == 0.) CYCLE
          psa=sl*sll*sl1*sl2
          sumin=0.
          sumin=dr*psa*sign
          aiii=aiii+dr*psa*sign
       ENDDO
    ENDDO
    vector_trcoefficients=16.0d0*pi*pi*aiii
    sggn=1.
    ix=(l1+l2-l-ll)/2
    ixx=ix-2*(ix/2)
    IF (ixx /= 0) sggn=-1.
    vector_trcoefficients=vector_trcoefficients/(ak*akk*ak1*ak2)*sggn

  END FUNCTION  vector_trcoefficients

  !  Spherical harmonics from Num. Recipes 

  DOUBLE PRECISION FUNCTION spherical_harmonics(m1,l,x)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m1, l
    DOUBLE PRECISION, INTENT(IN) ::  x
    DOUBLE PRECISION, DIMENSION(0:51) :: y
    INTEGER :: iphase, m, j
    DOUBLE PRECISION :: fj, z, fac, div, sum, a, b, c
    spherical_harmonics=0.
    m=IABS(m1)
    IF(m.LT.0) m=-m1
    y(0)=1.
    IF(l.EQ.0) THEN
       sum=y(0)
    ELSE
       a=m-l
       b=l+m+1
       c=m+1
       z=0.5-x*0.5
       DO j=1,l-m+1
          fj=j-1
          y(j)=y(j-1)*(a+fj)*(b+fj)*z
          div=(c+fj)*(fj+1.)
          y(j)=y(j)/div
       ENDDO
       IF(m > 0) then
          fac=(1.-x*x)**m
          fac=SQRT(fac)
       ELSE
          fac=1.
       ENDIF
       sum=0.
       DO j=0,l-m
          sum=sum+y(j)
       ENDDO
       iphase=m
       IF(m1.LT.0) then
          iphase=0
       ENDIF
       sum=sum*fac*((-1)**iphase)
    ENDIF
    spherical_harmonics=cn(m,l)*sum

  END FUNCTION spherical_harmonics



END MODULE ang_mom_functions




MODULE single_particle_orbits
  TYPE, PUBLIC :: single_particle_descript
     INTEGER :: total_orbits
     INTEGER, DIMENSION(:), POINTER :: nn, ll, jj,jz, itzp,jorder,obst
     REAL*8, DIMENSION(:), POINTER :: e
  END TYPE single_particle_descript

  TYPE (single_particle_descript), PUBLIC :: msp,jjsp
  
CONTAINS
  SUBROUTINE allocate_sp_array(this_array,n)
    TYPE (single_particle_descript), INTENT(INOUT) :: this_array
    INTEGER , INTENT(IN) :: n
    integer :: i 

    IF (ASSOCIATED (this_array%nn) ) DEALLOCATE(this_array%nn)
    ALLOCATE(this_array%nn(n))
    IF (ASSOCIATED (this_array%ll) ) DEALLOCATE(this_array%ll)
    ALLOCATE(this_array%ll(n))
    IF (ASSOCIATED (this_array%jj) ) DEALLOCATE(this_array%jj)
    ALLOCATE(this_array%jj(n))
    IF (ASSOCIATED (this_array%itzp) ) DEALLOCATE(this_array%itzp)
    ALLOCATE(this_array%itzp(n))

    IF (ASSOCIATED (this_array%jz) ) DEALLOCATE(this_array%jz)
    ALLOCATE(this_array%jz(n))

    IF (ASSOCIATED (this_array%jorder) ) DEALLOCATE(this_array%jorder)
    ALLOCATE(this_array%jorder(n))

    IF (ASSOCIATED (this_array%e) ) DEALLOCATE(this_array%e)
    ALLOCATE(this_array%e(n))
    IF (ASSOCIATED (this_array%obst) ) DEALLOCATE(this_array%obst)
    ALLOCATE(this_array%obst(n))
    !           blank all characters and zero all other values

    DO i= 1, n
       this_array%e(i)=0.
       this_array%nn(i)=0
       this_array%ll(i)=0
       this_array%jj(i)=0
       this_array%jz(i)=0
       this_array%jorder(i)=0
       this_array%itzp(i)=0
       this_array%obst(i)=0
    ENDDO

  END SUBROUTINE allocate_sp_array
  
  SUBROUTINE deallocate_sp_array(this_array)
    TYPE (single_particle_descript), INTENT(INOUT) :: this_array
    DEALLOCATE(this_array%nn) ; DEALLOCATE(this_array%ll)
    DEALLOCATE(this_array%jj) ;DEALLOCATE(this_array%itzp)
    DEALLOCATE(this_array%e) ;DEALLOCATE(this_array%jz)
    DEALLOCATE(this_array%jorder) 
    DEALLOCATE(this_array%obst) 

  END SUBROUTINE deallocate_sp_array
END MODULE single_particle_orbits

 module c_funcs
   use iso_c_binding
   implicit none
 interface
     subroutine qsort(array,elem_count,elem_size,compare) bind(C,name="qsort")
     use iso_c_binding
     implicit none
 
       type(c_ptr),value       :: array
       integer(c_size_t),value :: elem_count
       integer(c_size_t),value :: elem_size
       type(c_funptr),value    :: compare !int(*compare)(const void *, const void *)
     end subroutine qsort !standard C library qsort
   end interface
 end module c_funcs



module m_int
  use c_funcs
  implicit none
  type, bind(c),public ::  mat
    integer(C_INT128_T) :: idx
    real(C_DOUBLE) :: val
  end type 
     interface operator(<=)
        module procedure LessOrEqual
     end interface operator(<=)


  integer:: nj1b,nj2b,nm1b,nm2b
  type(mat),allocatable,target  :: jj1b(:),jj2b(:)
  type(mat),allocatable,target  :: m1b(:),m2b(:)



  contains

        function compare(v1, v2) bind(C)
           integer(C_INT) compare
           type(C_PTR), value :: v1, v2
           type(mat), pointer :: p1, p2
  
           call C_F_POINTER(v1,p1)
           call C_F_POINTER(v2,p2)
           if(.NOT.(p1 <= p2)) then
              compare = 1
           else if(.NOT.(p2 <= p1)) then
              compare = -1
           else
              compare = 0
           end if

        end function compare

    function LessOrEqual(v1, v2)
       logical LessOrEqual
       type(mat), intent(in) :: v1, v2
    
       LessOrEqual = v1%idx <= v2%idx
    end function LessOrEqual



    subroutine fetch_mat(m,k,val)
      integer*16, intent(in) :: m
      integer,intent(in) :: k
      real*8, intent(inout):: val
      type(mat), pointer :: pt(:)
      integer :: low,high,mid,flag
      val=0.d0
      low=1
      flag=0
      select case(k)
      case(1)
        pt=>jj1b
        high=nj1b
      case(2)
        pt=>jj2b
        high=nj2b
      end select

      if(pt(1)%idx==m)then
        val=pt(1)%val
        return
      endif

      if(pt(high)%idx==m)then
        val=pt(high)%val
        return
      endif


      do
      mid=(low+high)/2
      if(mid==low)exit

      if(pt(mid)%idx==m)then
        val=pt(mid)%val
        exit
      elseif(pt(mid)%idx>m)then
        high=mid
      elseif(pt(mid)%idx<m)then
        low=mid
      endif
      
      enddo
      return
      



    end subroutine fetch_mat


    subroutine sort_mat(k)

           integer(c_size_t) :: n
           integer,intent(in):: k



           if(k==1)then
             n=nj1b
             call qsort(c_loc(jj1b(1)),n,sizeof(jj1b(1)),c_funloc(compare))
           else if(k==2)then
             n=nj2b
             call qsort(c_loc(jj2b(1)),n,sizeof(jj2b(1)),c_funloc(compare))
           else if(k==3)then
             call qsort(c_loc(m1b(1)),n,sizeof(m1b(1)),c_funloc(compare))
           elseif(k==4)then
             call qsort(c_loc(m2b(1)),n,sizeof(m2b(1)),c_funloc(compare))
           endif
   
    end subroutine sort_mat

  end module m_int



