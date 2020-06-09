subroutine sm_matrix
  use mpi_mapping
  use lanczos
  use m_sp

  implicit none
  integer:: bra_min,bra_max,ket_min,ket_max
  integer:: bra,ket
  
  interface

    function moverlap(a,b)
      real*8:: moverlap
     integer,intent(in)::a,b
    end function moverlap

  end interface


   bra_min=mapping(sub_row,1)
   bra_max=mapping(sub_row,2)
   ket_min=mapping(sub_col,1)
   ket_max=mapping(sub_col,2)

   do bra=bra_min,bra_max
     do ket=ket_min,ket_max
       haml(bra,ket)=moverlap(bra,ket)
     enddo
   enddo

end subroutine sm_matrix
 


 function moverlap(a,b)
   use m_sp
   use constants
   use lanczos
   use tiny_functions
   use storage

   implicit none
   real*8 :: moverlap
   integer,intent(in)::a,b
   integer::channel,channel1,channel2,ii,ij
   integer:: i,j,k,l,il,id,ie,ik
   integer*16::ia,ib,ic
   real*8:: tmp,tmp1,val,val1
   integer*16:: bra_p,ket_p,bra_n,ket_n,c1,c2,c3,c4,nucbit,ucbit,bit_back
   integer:: bra,ket
   integer*16:: nxor,one
  integer:: bt_maskp(20),bt_maskn(20),bt_mask(40),bit_change(6),ncbitn,ncbitp,bit_uc(40)
  integer :: sig
  moverlap=0.d0

   bra_p=mbsp(mbs_row(a,1))
   ket_p=mbsp(mbs_col(b,1))
   ket_n=mbsn(mbs_col(b,2))
   bra_n=mbsn(mbs_row(a,2))




      c1=xor(bra_p,ket_p)
      c2=xor(bra_n,ket_n)
      nxor=popcnt(c1)+popcnt(c2)
      if(nxor>6)return

      one=01

      bt_mask=0
      bt_maskn=0
      bt_maskp=0

      tmp=0.d0
      tmp1=0

      select case(nxor)
      case (0)
        ! diagonal
        ! add 1b 2b 3b
        c1=bra_p
        c2=bra_n
        bt_maskn=0
        bt_maskp=0
        bt_mask=0
        ib=0;ic=0
        !    write(*,*)c1,c2
        ib=0;ic=0
        do ia=1,max(nn_orb,np_orb)
          if(and(c1, one)== one)then
            ib=ib+1
            bt_mask(ib+nn)=ia+nn_orb
          endif
          if(and(c2, one)== one)then
            ic=ic+1
            bt_mask(ic)=ia
            !       write(*,*)'hi',ic,c2
          endif
          c1=ishft(c1,-1)
          c2=ishft(c2,-1)
        enddo

        if(ib/=np .or. ic/= nn)then
          write(*,*)bra_n,bra_p
          write(*,*)ib,ic,np,nn
        write(*,*)'fatal error in case 0'
        stop
         endif
        ic=ib+ic


        !write(*,*)'ibic',ib,ic,bra_p,bra_n


        ! add 1b interaction
        do id=1,ic
          i=bt_mask(id)
          val=vt1(i,i)
          tmp=tmp+val
        enddo


        ! add 2b interaction <ij|V|ij>
        do ib=1,ic-1
          do id=ib+1,ic
            i=bt_mask(ib)
            j=bt_mask(id)
            channel=ip2(i,j,1)
            bra=ip2(i,j,2)
            val=vt2(channel)%val(bra,bra)
            tmp=tmp+val
          enddo
        enddo

        ! add 3b interaction <ijk|V|ijk>
        if(number_interactions==3)then
          do ib=1,ic-2
            do id=ib+1,ic-1
              do ie=id+1,ic

                i=bt_mask(ib)
                j=bt_mask(id)
                k=bt_mask(ie)
                channel=ip3(i,j,k,1)
                bra=ip3(i,j,k,2)
                if(channel*bra==0)cycle
                val=vt3(channel)%val(bra,bra)
                tmp=tmp+val
              enddo
            enddo
          enddo
        endif

      case (2)
        bt_mask=0
        ! onebody+ twobody+ 3b
        c1=xor(bra_p , ket_p)
        c2=xor(bra_n , ket_n)

        if(c1/=0)then
          if(c2/=0)stop 'error in case 2'
          ucbit=and(not(c1) ,bra_p)
          bit_back=ucbit
          c2=and(c1,ket_p)
          c1=and(c1,bra_p)
          i=0; j=0
          do ia=1,np_orb
            if(and(c1,one)==one)then
              if(i/=0)write(*,*)'fatal error'
              i=ia
            endif
            if(and(c2,one)==one)then
              if(j/=0)write(*,*)'fatal error'
              j=ia
            endif
            c1=ishft(c1,-1)
            c2=ishft(c2,-1)
          enddo

          i=i+nn_orb
          j=j+nn_orb
          ! i,j is the chaged bit in bra and ket

          ib=0
          bt_mask=0
          do ia=1,np_orb
            if(and(ucbit, one)== one)then
              ib=ib+1
              bt_mask(ib+nn)=ia+nn_orb
            endif
            ucbit=ishft(ucbit,-1)
          enddo
          nucbit=ib

          ucbit=bra_n
          ib=0
          do ia=1,nn_orb
            if(and(ucbit, one)== one)then
              ib=ib+1
              bt_mask(ib)=ia
            endif
            ucbit=ishft(ucbit,-1)
          enddo
          nucbit=nucbit+ib





        else
          if(c2==0)stop 'case 2'
          ucbit=and(not(c2), bra_n)
          bit_back=ucbit
          c1=and(c2,bra_n)
          c2=and(c2,ket_n)
          i=0; j=0
          do ia=1,nn_orb
            if(and(c1,one)==one)then
              if(i/=0)write(*,*)'fatal error'
              i=ia
            endif
            if(and(c2,one)==one)then
              if(j/=0)write(*,*)'fatal error'
              j=ia
            endif
            c1=ishft(c1,-1)
            c2=ishft(c2,-1)
          enddo


          ib=0
          bt_mask=0
          do ia=1,nn_orb
            if(and(ucbit, one)== one)then
              ib=ib+1
              bt_mask(ib)=ia
            endif
            ucbit=ishft(ucbit,-1)
          enddo

          ucbit=bra_p
          do ia=1,np_orb
            if(and(ucbit, one)== one)then
              ib=ib+1
              bt_mask(ib)=ia+nn_orb
            endif
            ucbit=ishft(ucbit,-1)
          enddo
          nucbit=ib


        endif
        


        c1=xor(bra_p , ket_p)
        c2=xor(bra_n , ket_n)
        ! 1b

        sig=0
        if(c1/=0)then
        sig=ij_bits(bra_p,0,i-nn_orb)
        sig=ij_bits(ket_p,0,j-nn_orb)+sig
        else if(c2/=0)then
        sig=ij_bits(bra_n,0,i)
        sig=ij_bits(ket_n,0,j)+sig
        else
        write(*,*)'error in case 2'
        stop
        endif
        

        val=0.0
        val=vt1(i,j)
        tmp=tmp+val*(-1)**sig

        ! 2b

        do ia=1,nucbit
          k=bt_mask(ia)
          ii=min(i,k)
          ij=max(i,k)
          channel=ip2(ii,ij,1)
          if(channel==0)cycle
          bra=ip2(ii,ij,2)
          sig=0
          sig=ch_sign(bra_n,bra_p,ii,ij)


          ii=min(j,k)
          ij=max(j,k)

          sig=sig+ch_sign(ket_n,ket_p,ii,ij)
          if(ip2(ii,ij,1)/=channel)cycle
          ket=ip2(ii,ij,2)
          if(bra*ket==0)cycle
          val=vt2(channel)%val(bra,ket)
          tmp=tmp+val*(-1.0)**sig
        enddo

        ! 3b
        if(number_interactions==3)then
          do ia=1,nucbit-1
            k=bt_mask(ia)
            do ic=ia+1,nucbit
              l=bt_mask(ic)
              if(k>l)stop 'fatal error 3b'
              sig=0
              if(i<k)then
                ii=i
                ij=k
                ik=l
              else if(i>l)then
                ii=k
                ij=l
                ik=i
              else
                ii=k
                ij=i
                ik=l
              endif
              sig=ch_sign(bra_n,bra_p,0,ii)
              sig=ch_sign(bra_n,bra_p,ij,ik)+sig
              channel=ip3(ii,ij,ik,1)
              if(channel==0)cycle
              bra=ip3(ii,ij,ik,2)

              if(j<k)then
                ii=j
                ij=k
                ik=l
              elseif(j>l)then
                ii=k
                ij=l
                ik=j
              else
                ii=k
                ij=j
                ik=l
              endif
              if(channel/=ip3(ii,ij,ik,1))cycle
              ket=ip3(ii,ij,ik,2)
              if(bra*ket==0)cycle

              sig=ch_sign(ket_n,ket_p,0,ii)+sig
              sig=ch_sign(ket_n,ket_p,ij,ik)+sig
              val=vt3(channel)%val(bra,ket)
              tmp=tmp+val*(-1)**sig
            enddo
          enddo
        endif








      case (4)
        ! twobdy+3b
        bit_uc=0


        ! twobody
        c1=xor(bra_p , ket_p)
        c2=xor(bra_n , ket_n)
        ! get unchanged bit
        c3= and(not(c1) , bra_p)
        c4= and(not(c2) , bra_n)
        ncbitn=0
        do ia=1,np_orb
          if(and(c3,one) == one )then
            ncbitn=ncbitn+1
            bit_uc(ncbitn)=ia+nn_orb
          endif
          c3=ishft(c3,-1)
        enddo
        do ia=1,nn_orb
          if(and(c4,one) == one)then
            ncbitn=ncbitn+1
            bit_uc(ncbitn)=ia
          endif
          c4=ishft(c4,-1)
        enddo





        ! add twobody interaction
        c1=xor(bra_p , ket_p)
        c2=xor(bra_n , ket_n)

        c3=and(c1,bra_p)
        c4=and(c2,bra_n)
        ncbitp=0
        bit_change=0
        do ia=1,np_orb
          if(and(c3,one)==one)then
            ncbitp=ncbitp+1
            bit_change(ncbitp)=ia+nn_orb
          endif
          c3=ishft(c3,-1)
        enddo

        do ia=1,nn_orb
          if(and(c4,one)==one)then
            ncbitp=ncbitp+1
            bit_change(ncbitp)=ia
          endif
          c4=ishft(c4,-1)
        enddo

        if(ncbitp/=2)then
          write(*,*)'error in case 4 <bra|ket>, a', ncbitp, and(c2,bra_n),and(c1,bra_p)
          stop
        endif

        c3=and(c1,ket_p)
        c4=and(c2,ket_n)
        do ia=1,np_orb
          if(and(c3,one)==one)then
            ncbitp=ncbitp+1
            bit_change(ncbitp)=ia+nn_orb
          endif
          c3=ishft(c3,-1)
        enddo
        do ia=1,nn_orb
          if(and(c4,one)==one)then
            ncbitp=ncbitp+1
            bit_change(ncbitp)=ia
          endif
          c4=ishft(c4,-1)
        enddo
        if(ncbitp/=4)write(*,*)'error in case 4 <bra|ket> b'



        sig=0


        ii=min(bit_change(1),bit_change(2))
        ij=max(bit_change(1),bit_change(2))
        ik=min(bit_change(3),bit_change(4))
        il=max(bit_change(3),bit_change(4))
        sig=ch_sign(bra_n,bra_p,ii,ij)+ch_sign(ket_n,ket_p,ik,il)

        channel1=ip2(ii,ij,1)
        channel2=ip2(ik,il,1)
        if(channel1==channel2)then

        bra=ip2(ii,ij,2)
        ket=ip2(ik,il,2)

        val=vt2(channel1)%val(bra,ket)
        tmp=tmp+val*(-1.0)**sig
        endif

        sig=0

        ! add 3b
        if(number_interactions==3)then
          do ia=1,ncbitn
            val=0.d0
            sig=0
            j=bit_uc(ia)

            if(j<ii)then
              channel1=ip3(j,ii,ij,1)
              bra=ip3(j,ii,ij,2)
            sig=sig+ch_sign(bra_n,bra_p,0,j)
              sig=ch_sign(bra_n,bra_p,ii,ij)+sig
            elseif(j>ij)then
              channel1=ip3(ii,ij,j,1)
              bra=ip3(ii,ij,j,2)
            sig=sig+ch_sign(bra_n,bra_p,0,ii)
              sig=ch_sign(bra_n,bra_p,ij,j)+sig
            else
              channel1=ip3(ii,j,ij,1)
              bra=ip3(ii,j,ij,2)
            sig=sig+ch_sign(bra_n,bra_p,0,ii)
              sig=ch_sign(bra_n,bra_p,j,ij)+sig
            endif

            if(j<ik)then
              channel2=ip3(j,ik,il,1)
              ket=ip3(j,ik,il,2)

            sig=sig+ch_sign(ket_n,ket_p,0,j)
            sig=sig+ch_sign(ket_n,ket_p,ik,il)
            elseif(j>il)then
              channel2=ip3(ik,il,j,1)
              ket=ip3(ik,il,j,2)
            sig=sig+ch_sign(ket_n,ket_p,0,ik)
            sig=sig+ch_sign(ket_n,ket_p,il,j)
            else
              channel2=ip3(ik,j,il,1)
              ket=ip3(ik,j,il,2)
            sig=sig+ch_sign(ket_n,ket_p,0,ik)
            sig=sig+ch_sign(ket_n,ket_p,j,il)
            endif
            if(channel1/=channel2)cycle
            if(bra*ket==0)cycle

           val=vt3(channel1)%val(bra,ket)
            tmp=tmp+val*(-1)**sig
          enddo
        endif

      case (6)

        if(number_interactions==3)then
          c1=xor(bra_p , ket_p)
          c2=xor(bra_n , ket_n)
          c1=and(c1,bra_p)
          c2=and(c2,bra_n)
          ie=0
          do ia=1,nn_orb
            if(and(c2,one)==one)then
              ie=ie+1
              bit_change(ie)=ia
            endif
            c2=ishft(c2,-1)
          enddo

          do ia=1,np_orb
            if(and(c1,one)==one)then
              ie=ie+1
              bit_change(ie)=ia+nn_orb
            endif
            c1=ishft(c1,-1)
          enddo
          if(ie/=3)stop 'erro in case 6'

          c1=xor(bra_p , ket_p)
          c2=xor(bra_n , ket_n)
          c1=and(c1 ,ket_p)
          c2=and(c2 ,ket_n)
          do ia=1,nn_orb
            if(and(c2,one)==one)then
              ie=ie+1
              bit_change(ie)=ia
            endif
            c2=ishft(c2,-1)
          enddo

          do ia=1,np_orb
            if(and(c1,one)==one)then
              ie=ie+1
              bit_change(ie)=ia+nn_orb
            endif
            c1=ishft(c1,-1)
          enddo


          if(ie/=6)write(*,*)'fatal error in <bra| ket> 3'

          sig=0
            ii=bit_change(1)
            ij=bit_change(2)
            ik=bit_change(3)
            channel1=ip3(ii,ij,ik,1)
            bra=ip3(ii,ij,ik,2)
        
            sig=sig+ch_sign(bra_n,bra_p,0,bit_change(1))
            sig=sig+ch_sign(bra_n,bra_p,bit_change(2),bit_change(3))

            ii=bit_change(4)
            ij=bit_change(5)
            ik=bit_change(6)
            channel2=ip3(ii,ij,ik,1)
            ket=ip3(ii,ij,ik,2)

            sig=sig+ch_sign(ket_n,ket_p,0,bit_change(4))
            sig=sig+ch_sign(ket_n,ket_p,bit_change(5),bit_change(6))
            if((channel1==channel2) .and. (bra*ket/=0))then

              val=vt3(channel1)%val(bra,ket)
              tmp=tmp+val*(-1)**sig
            endif
        endif

        ! pure 3b
      end select
      moverlap=tmp

!@!      write(*,*)a,b,nxor,moverlap

 end function moverlap


