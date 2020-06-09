program main
  use m_int
  use single_particle_orbits
  use constants
  use ang_mom_functions
  implicit none
  integer:: njjsp,nmsp,i,j,ctr
  integer*16:: a,b,c,d,vidx,aa,bb,cc,dd,ee,ff
  integer:: ja,jb,jc,jd,ma,mb,mc,md,bra,ket
  integer:: jt,ipar,iso_z,jtot
  integer,allocatable:: mconfig(:)
  integer:: nmconfig
  real*8:: angmon_fact,val,vtmp
  real*8:: dij
  logical:: hermit
  integer:: iph
  integer:: ph_type
  integer:: phase_ab,phase_cd
  character(len=20):: ctmp
  character(len=100):: jjsp_file ,&
              jj1b_file ,&
              jj2b_file ,&
              jj3b_file ,&
              msp_file ,&
              m1b_file ,&
              m2b_file ,&
              m3b_file 



  open(unit=5,file='input.dat')

  read(5,*)
  read(5,*)hermit
  read(5,*)
  read(5,*)jjsp_file

   call commons_to_angmom
  open(unit=10,file=jjsp_file)
  read(10,*)njjsp
  call allocate_sp_array(jjsp,njjsp)
  nmsp=0

  do i=1,njjsp
    read(10,*)j,jjsp%nn(i),&
      jjsp%ll(i),&
      jjsp%jj(i),&
      jjsp%itzp(i),ctmp
     if(ctmp=='hole')then
             jjsp%obst(i)=-1
     else
             jjsp%obst(i)=1
     endif
    nmsp=nmsp+jjsp%jj(i)+1
  enddo

  call allocate_sp_array(msp,nmsp)

  ctr=0

  do i=1,njjsp
    if(jjsp%itzp(i)/=1)cycle
    do j=-jjsp%jj(i),jjsp%jj(i),2
      ctr=ctr+1
      msp%nn(ctr)=jjsp%nn(i)
      msp%ll(ctr)=jjsp%ll(i)
      msp%jj(ctr)=jjsp%jj(i)
      msp%itzp(ctr)=jjsp%itzp(i)
      msp%jorder(ctr)=i
      msp%obst(ctr)=jjsp%obst(i)
      msp%jz(ctr)=j
      if(msp%obst(ctr)==-1)then
      msp%jz(ctr)=-j
      endif

    enddo
  enddo


  do i=1,njjsp
    if(jjsp%itzp(i)/=-1)cycle
    do j=-jjsp%jj(i),jjsp%jj(i),2
      ctr=ctr+1
      msp%nn(ctr)=jjsp%nn(i)
      msp%ll(ctr)=jjsp%ll(i)
      msp%jj(ctr)=jjsp%jj(i)
      msp%itzp(ctr)=jjsp%itzp(i)
      msp%jz(ctr)=j
      msp%jorder(ctr)=i
      msp%jz(ctr)=j
      msp%obst(ctr)=jjsp%obst(i)
      if(msp%obst(ctr)==-1)then
      msp%jz(ctr)=-j
      endif
    enddo
  enddo

  close(10)

  read(5,*)
  read(5,*)msp_file
  open(unit=15,file=msp_file)
  do i=1,ctr
     ctmp='particle'
    if(msp%obst(i)==-1)ctmp='hole'
    write(15,'(6I4,2X,A20)')i,msp%nn(i),&
      msp%ll(i),&
      msp%jj(i),&
      msp%jz(i),&
      msp%itzp(i),&
      ctmp
  enddo

close(15)

  if(ctr/= nmsp) stop 'error number of single particle m orbit'

read(5,*)
read(5,*)jj1b_file
  open(unit=11,file=jj1b_file)
  read(11,*)nj1b
  allocate(jj1b(nj1b))
  do i=1,nj1b
    read(11,*)a,b,jj1b(i)%val
    jj1b(i)%idx=ishft(a,7)+b
  enddo
  call sort_mat(1)
  close(11)



read(5,*)
read(5,*)m1b_file
  open(unit=13,file=m1b_file)
  ctr=0
  do a=1,nmsp
    aa=msp%jorder(a)
    do b=1,nmsp
      bb=msp%jorder(b)
       if(msp%ll(a)/=msp%ll(b))cycle
       if(msp%jj(a)/=msp%jj(b))cycle
       if(msp%jz(a)/=msp%jz(b))cycle
       ctr=ctr+1

    enddo
  enddo
  write(13,'(I4)')ctr

  do a=1,nmsp
    aa=msp%jorder(a)
    do b=1,nmsp
      bb=msp%jorder(b)
       if(msp%ll(a)/=msp%ll(b))cycle
       if(msp%jj(a)/=msp%jj(b))cycle
       if(msp%jz(a)/=msp%jz(b))cycle

      vidx=ishft(aa,7)+bb
      call fetch_mat(vidx,1,val)
      write(13,'(2I4,F12.6)')a,b,val
    enddo
  enddo
  close(13)





read(5,*)
read(5,*)jj2b_file

  open(unit=12,file=jj2b_file)
  read(12,*)nj2b
  allocate(jj2b(nj2b))

  jmin=100
  jmax=0
  do i=1,nj2b
    read(12,*)a,b,c,d,jt,jj2b(i)%val
    jj2b(i)%idx=ishft(a,35)+ishft(b,28)+ishft(c,14)+ishft(d,7)+jt
    if(jmin>jt)jmin=jt
    if(jmax<jt)jmax=jt
  enddo

  call sort_mat(2)
  close(12)

read(5,*)
read(5,*)m2b_file

  open(unit=14,file=m2b_file)

  !! pp

    ctr=0
  do jt=-jmax,jmax
    do ipar=0,1
      do iso_z=-1,1

        nmconfig=0
        do a=1,nmsp
          if(msp%obst(a)==-1)cycle
          do b=a,nmsp
          if(msp%obst(b)==-1)cycle
            if(a==b)cycle
            if(msp%jz(a)+msp%jz(b) /=jt*2 )cycle
            if(mod(msp%ll(a)+msp%ll(b),2)/=ipar )cycle
            if(msp%itzp(a)+msp%itzp(b) /=2*iso_z )cycle
            nmconfig=nmconfig+1
          enddo
        enddo
        ctr=nmconfig*nmconfig+ctr

        enddo
      enddo
    enddo
    write(14,'(I8)')ctr


  do jt=-jmax,jmax
    do ipar=0,1
      do iso_z=-1,1

        nmconfig=0
        do a=1,nmsp
          if(msp%obst(a)==-1)cycle
          do b=a,nmsp
          if(msp%obst(b)==-1)cycle
            if(a==b)cycle
            if(msp%jz(a)+msp%jz(b) /=jt*2 )cycle
            if(mod(msp%ll(a)+msp%ll(b),2)/=ipar )cycle
            if(msp%itzp(a)+msp%itzp(b) /=2*iso_z )cycle
            nmconfig=nmconfig+1
          enddo
        enddo
        allocate(mconfig(nmconfig*2))

        nmconfig=0
        do a=1,nmsp
          if(msp%obst(a)==-1)cycle
          do b=a,nmsp
          if(msp%obst(b)==-1)cycle
            if(a==b)cycle
            if(msp%jz(a)+msp%jz(b) /=jt*2 )cycle
            if(mod(msp%ll(a)+msp%ll(b),2)/=ipar )cycle
            if(msp%itzp(a)+msp%itzp(b) /=2*iso_z )cycle
            nmconfig=nmconfig+1
            mconfig(2*nmconfig-1)= a
            mconfig(2*nmconfig)  = b
          enddo
        enddo

        do bra=1,nmconfig
          a=mconfig(2*bra-1)
          b=mconfig(2*bra)
          ma=msp%jz(a)
          mb=msp%jz(b)
          aa=msp%jorder(a)
          bb=msp%jorder(b)
          ja=jjsp%jj(aa)
          jb=jjsp%jj(bb)
          do ket=1,nmconfig
            c=mconfig(2*ket-1)
            d=mconfig(2*ket)
            mc=msp%jz(c)
            md=msp%jz(d)
            cc=msp%jorder(c)
            dd=msp%jorder(d)
            jc=jjsp%jj(cc)
            jd=jjsp%jj(dd)


            jmin=max(abs(ja-jb)/2,abs(jc-jd)/2)
            jmax=min(abs(ja+jb)/2,abs(jc+jd)/2)
            vtmp=0.d0
            if(jmin<=jmax)then
            do jtot=jmin,jmax
              angmon_fact=cgc(ja,jb,2*jtot,ma,mb,2*jt)*&
                cgc(jc,jd,2*jtot,mc,md,2*jt)
                phase_ab=1
                phase_cd=1


                if(iso_z/=0)then
                if(aa>bb)then
                  ee=aa
                  aa=bb
                  bb=ee
                  phase_ab=iph((ja+jb)/2+jtot+1)
                endif
                if(cc>dd)then
                  ee=cc
                  cc=dd
                  dd=ee
                  phase_cd=iph((jc+jd)/2+jtot+1)
                endif
                 endif

              vidx=ishft(aa,35)+ishft(bb,28)+ishft(cc,14)+ishft(dd,7)+jtot

              call fetch_mat(vidx,2,val)
              vtmp=vtmp+val*angmon_fact*phase_ab*phase_cd
            enddo
          endif
            write(14,'(4I4,F12.8)')a,b,c,d,vtmp
        enddo
      enddo

      deallocate(mconfig)


    enddo
  enddo
enddo



  !! hh

    ctr=0
  do jt=-jmax,jmax
    do ipar=0,1
      do iso_z=-1,1

        nmconfig=0
        do a=1,nmsp
          if(msp%obst(a)==1)cycle
          do b=a,nmsp
          if(msp%obst(b)==1)cycle
            if(a==b)cycle
            if(msp%jz(a)+msp%jz(b) /=jt*2 )cycle
            if(mod(msp%ll(a)+msp%ll(b),2)/=ipar )cycle
            if(msp%itzp(a)+msp%itzp(b) /=2*iso_z )cycle
            nmconfig=nmconfig+1
          enddo
        enddo
        ctr=nmconfig*nmconfig+ctr

        enddo
      enddo
    enddo

    write(14,'(I8)')ctr


  do jt=-jmax,jmax
    do ipar=0,1
      do iso_z=-1,1

        nmconfig=0
        do a=1,nmsp
          if(msp%obst(a)==1)cycle
          do b=a,nmsp
          if(msp%obst(b)==1)cycle
            if(a==b)cycle
            if(msp%jz(a)+msp%jz(b) /=jt*2 )cycle
            if(mod(msp%ll(a)+msp%ll(b),2)/=ipar )cycle
            if(msp%itzp(a)+msp%itzp(b) /=2*iso_z )cycle
            nmconfig=nmconfig+1
          enddo
        enddo
        allocate(mconfig(nmconfig*2))

        nmconfig=0
        do a=1,nmsp
          if(msp%obst(a)==1)cycle
          do b=a,nmsp
          if(msp%obst(b)==1)cycle
            if(a==b)cycle
            if(msp%jz(a)+msp%jz(b) /=jt*2 )cycle
            if(mod(msp%ll(a)+msp%ll(b),2)/=ipar )cycle
            if(msp%itzp(a)+msp%itzp(b) /=2*iso_z )cycle
            nmconfig=nmconfig+1
            mconfig(2*nmconfig-1)= a
            mconfig(2*nmconfig)  = b
          enddo
        enddo

        do bra=1,nmconfig
          a=mconfig(2*bra-1)
          b=mconfig(2*bra)
          ma=msp%jz(a)
          mb=msp%jz(b)
          aa=msp%jorder(a)
          bb=msp%jorder(b)
          ja=jjsp%jj(aa)
          jb=jjsp%jj(bb)
          do ket=1,nmconfig
            c=mconfig(2*ket-1)
            d=mconfig(2*ket)
            mc=msp%jz(c)
            md=msp%jz(d)
            cc=msp%jorder(c)
            dd=msp%jorder(d)
            jc=jjsp%jj(cc)
            jd=jjsp%jj(dd)


            jmin=max(abs(ja-jb)/2,abs(jc-jd)/2)
            jmax=min(abs(ja+jb)/2,abs(jc+jd)/2)
            vtmp=0.d0
            if(jmin<=jmax)then
            do jtot=jmin,jmax
              angmon_fact=cgc(ja,jb,2*jtot,ma,mb,2*jt)*&
                cgc(jc,jd,2*jtot,mc,md,2*jt)
                phase_ab=1
                phase_cd=1


                if(iso_z/=0)then
                if(aa>bb)then
                  ee=aa
                  aa=bb
                  bb=ee
                  phase_ab=iph((ja+jb)/2+jtot+1)
                endif
                if(cc>dd)then
                  ee=cc
                  cc=dd
                  dd=ee
                  phase_cd=iph((jc+jd)/2+jtot+1)
                endif
                 endif

              vidx=ishft(aa,35)+ishft(bb,28)+ishft(cc,14)+ishft(dd,7)+jtot

              call fetch_mat(vidx,2,val)
              vtmp=vtmp+val*angmon_fact*phase_ab*phase_cd
            enddo
          endif
            write(14,'(4I4,F12.8)')a,b,c,d,vtmp
        enddo
      enddo

      deallocate(mconfig)


    enddo
  enddo
enddo


!ph


    ctr=0
  do jt=-jmax,jmax
    do ipar=0,1
      do iso_z=-1,1

        nmconfig=0
        do a=1,nmsp
          do b=a,nmsp
          if(msp%obst(b)*msp%obst(a)/=-1)cycle
            if(a==b)cycle
            if(msp%jz(a)+msp%jz(b) /=jt*2 )cycle
            if(mod(msp%ll(a)+msp%ll(b),2)/=ipar )cycle
            if(msp%itzp(a)+msp%itzp(b) /=2*iso_z )cycle
            nmconfig=nmconfig+1
          enddo
        enddo
        ctr=nmconfig*nmconfig+ctr

        enddo
      enddo
    enddo

    write(14,'(I8)')ctr


  do jt=-jmax,jmax
    do ipar=0,1
      do iso_z=-1,1

        nmconfig=0
        do a=1,nmsp
          do b=a,nmsp
          if(msp%obst(b)*msp%obst(a)/=-1)cycle
            if(a==b)cycle
            if(msp%jz(a)+msp%jz(b) /=jt*2 )cycle
            if(mod(msp%ll(a)+msp%ll(b),2)/=ipar )cycle
            if(msp%itzp(a)+msp%itzp(b) /=2*iso_z )cycle
            nmconfig=nmconfig+1
          enddo
        enddo
        allocate(mconfig(nmconfig*2))

        nmconfig=0
        do a=1,nmsp
          do b=a,nmsp
          if(msp%obst(b)*msp%obst(a)/=-1)cycle
            if(a==b)cycle
            if(msp%jz(a)+msp%jz(b) /=jt*2 )cycle
            if(mod(msp%ll(a)+msp%ll(b),2)/=ipar )cycle
            if(msp%itzp(a)+msp%itzp(b) /=2*iso_z )cycle
            nmconfig=nmconfig+1
            mconfig(2*nmconfig-1)= a
            mconfig(2*nmconfig)  = b
          enddo
        enddo

        do bra=1,nmconfig
          a=mconfig(2*bra-1)
          b=mconfig(2*bra)
          ma=msp%jz(a)
          mb=msp%jz(b)
          aa=msp%jorder(a)
          bb=msp%jorder(b)
          ja=jjsp%jj(aa)
          jb=jjsp%jj(bb)
          if(msp%obst(a)==-1)then 
                  ph_type=1
                  write(*,*)'hp'
          else
                  ph_type=-1
                  write(*,*)'ph'
          endif
          

          do ket=1,nmconfig
            c=mconfig(2*ket-1)
            d=mconfig(2*ket)
            mc=msp%jz(c)
            md=msp%jz(d)
            cc=msp%jorder(c)
            dd=msp%jorder(d)
            jc=jjsp%jj(cc)
            jd=jjsp%jj(dd)


            jmin=max(abs(ja-jb)/2,abs(jc-jd)/2)
            jmax=min(abs(ja+jb)/2,abs(jc+jd)/2)
            vtmp=0.d0
            if(jmin<=jmax)then
            do jtot=jmin,jmax
              if(ph_type==1)then !! hp
              angmon_fact=cgc(ja,jb,2*jtot,ma,mb,2*jt)*&
                cgc(jc,jd,2*jtot,mc,md,2*jt)

              vidx=ishft(bb,35)+ishft(aa,28)+ishft(dd,14)+ishft(cc,7)+jtot
              call fetch_mat(vidx,2,val)
              vtmp=vtmp+val*angmon_fact
             else
                     !! ph
              angmon_fact=cgc(ja,jb,2*jtot,ma,mb,2*jt)*&
                cgc(jc,jd,2*jtot,mc,md,2*jt)

              vidx=ishft(aa,35)+ishft(bb,28)+ishft(cc,14)+ishft(dd,7)+jtot

              call fetch_mat(vidx,2,val)
              vtmp=vtmp+val*angmon_fact
             endif
            enddo
          endif
            write(14,'(4I4,F12.8)')a,b,c,d,vtmp
        enddo
      enddo

      deallocate(mconfig)


    enddo
  enddo
enddo

close(14)



end program main
