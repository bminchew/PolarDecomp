program PolarDecomp
!****************************************************************


 
!****************************************************************
!**     
!**   FILE NAME: PolarDecomp.f
!**     
!**   DATE WRITTEN: April 2011
!**     
!**   PROGRAMMER: 
!**                        Brent Minchew  
!**               California Institute of Technology 
!**                       bminchew@caltech.edu
!**     
!**   FUNCTIONAL DESCRIPTION: This program decomposes UAVSAR-formatted
!**      quad-pol SAR data and outputs user-selected parameters. 
!**   
!**    Optional outputs:
!**     1)  H/A/alpha parameters         
!**     2)  Pauli decomposition parameters
!**     3)  Averaged intensity (I = p1*lam1 + p2*lam2 + p3*lam3)
!**     4)  Shannon entropy parameters  
!**     5)  Eigenvalues (lam1, lam2, lam3)
!**     6)  Pseudo probabilities (p1, p2, p3)
!**     7)  Additional eigenvalue parameters (beta, gamma, delta)
!**     8)  Invariants of coherency and (equivalently) covariance matrices 
!**   
!**   ROUTINES CALLED:  None
!**     
!**   NOTES: 
!**   
!**
!**   ASSUMPTIONS:  No additional assumptions in the program 
!**      
!**   PERTINENT LITERATURE:
!**
!**
!**
!**   *This program is an ugraded version of a program written 
!**     by the same programmer in April 2010 through funding 
!**     provided by:
!**                         UAVSAR Group  
!**                 NASA Jet Propulsion Laboratory 
!**      
!**   
!**     
!**   UPDATE LOG:
!**
!**   Date Changed        Reason Changed                  CR # and Version #
!**   ------------       ----------------                 -----------------
!**     
!*****************************************************************



   implicit none
!  Initialize variables for cmd file read
   character*250          infile,dataloc,otfile,otfld
   character*250          str,str2,label,labell,label2
   character*250          HHHHfile,VVVVfile,HVHVfile,HHHVfile
   character*250          HHVVfile,HVVVfile
   character*34           basenm
   character*6            trlnm
   character*3            uvsryn,dattyp
   integer                hhhhtrm,vvvvtrm,hvhvtrm,hhhvtrm,geteig,getshan
   integer                hhvvtrm,hvvvtrm,dattrm,otfldtrm,otfiltrm,otstrtrm
   integer                ios,ios2,lne,k,rdtype
   integer                cols,srow,erow,left,right
   real*4                 dbmul,winsize

!  Initialize output variables 
   character*500          otstr
   real*4,allocatable::   paul1(:),paul2(:),paul3(:)
   real*4,allocatable::   t11a(:),t22a(:),t33a(:)
   complex*8,allocatable:: t12a(:),t13a(:),t23a(:)
   complex*8,allocatable:: t12c(:),t13c(:),t23c(:)

!  Initialize working variables 
   real*4,allocatable::   hhhh(:,:),hvhv(:,:),vvvv(:,:)
   complex*8,allocatable:: hhhv(:,:),hhvv(:,:),hvvv(:,:)
   integer,allocatable::  zrmt(:)
   real*4                 hhhv_r,hhhv_i,hhvv_r,hhvv_i,hvvv_r,hvvv_i
   real*4                 hhhh_r,hvhv_r,vvvv_r
   real*8                 t11,t22,t33,winsizinv,t12_r,t12_i,t13_r,t13_i,t23_r,t23_i
   real*8                 paulval1,paulval2,paulval3
   complex*8              t12,t13,t23,hhhv_c,hhvv_c,hvvv_c               
   integer                iii,jjj,eee,mmm,mm,ggg,hhh,zzz,scattnum,winsizside,winsizeeff
   integer                reco,recocount,wlin,wcol,scol,winreadcount,ierr,lines,width
   integer                paulicol,pauliln,zrcnt,nmout
   integer                stv(13),today(3),now(3),colcnt 

!  Defiine parameters
   character(len=11),parameter:: fu='unformatted'
   character(len=6), parameter:: di='direct'
   character(len=7), parameter:: re='replace'

!  Initialize parameter structure
   type :: parstruc
      real*8                 ex,ex2,ex2s,sqrt3,log3i,pi,r2d,exp1,dbmul
      complex*16             nmcp,nmcn
   end type
   
   type(parstruc) parstr

!  Initialize yes/no string (ynstr) structure
   type :: decisval 
      character*3            ent,ani,alp,beta,gama,delta,phi,inten
      character*3            ev1,ev2,ev3,shan,sep,sei,trace,detr,inv2
      character*3            prob1,prob2,prob3,paul1,paul2,paul3
      character*3            alp1,alp2,alp3,beta1,beta2,beta3
      integer                geteigvec,getalp,getbeta,eline
   end type

   type(decisval) ynstr

   print *,' '; print *,' '
   if (iargc().eq.1)goto 987
      print*,'  ~~~   UAVSAR polarimetric decomposition   ~~~'
      print*,' '
      print*,' usage: PolarDecomp_bm name_of_input_command_file '
      print*,' '
      print*,'        input command file template: '
      print*,'            /home/bminchew/Utilities/PolarDecomp_bm.cmd '
      print *,' ';print *,' '
      stop


987 continue
   rdtype    = 1
   nmout     = 0
   geteig    = 0
   getshan   = 0
   ynstr%geteigvec = 0
   ynstr%getalp    = 0
   ynstr%getbeta   = 0
   ynstr%eline     = 0

!  Read command file
   call getarg(1,infile)
   
   open(3,file=infile,status='old')

   ios = 0
   lne = 0 
   do while (ios.eq.0)
   read(3,'(A)',iostat=ios) str


   if (ios.eq.0.and.str(1:3).eq.'   '.and.str(4:6).ne.'***') then
      lne = lne + 1
      if (rdtype.eq.1) then
         k = scan(str,'::')
         label = trim(str(4:k-1))
         str=str(k+2:)
      else
         k = scan(str,'?')
         label = str(4:k-1)
         k = scan(str,'::')
         labell = str(1:k)
         str = str(k+2:)
      endif

   if (k.ne.0) then

   select case (label)
   case('Columns (samples) in image')
      read(str,*,iostat=ios) cols
      if (cols.eq.3300) then
         dattyp = 'mlc'
      else 
         dattyp = 'grd'
      endif
   case('Starting line')
      read(str,*,iostat=ios) srow
   case('Ending line (enter number, end, or all)')
      if (trim(adjustl(str)).eq.'all') then
         erow = 999999
         srow = 1
      elseif (trim(adjustl(str)).eq.'end') then
         erow = 999999
      else
         read(str,*,iostat=ios) erow
         ynstr%eline = 1
      endif
   case('Starting column')
      read(str,*,iostat=ios) left
   case('Ending column (enter number, end, or all)')
      if (trim(adjustl(str)).eq.'all') then
         right = cols 
         left = 1
      elseif (trim(adjustl(str)).eq.'end') then
         right = cols
      else
         read(str,*,iostat=ios) right
      endif
      if (right.gt.cols)right=cols
      width = right - left + 1
   case('Window size (pixels per side)')
      read(str,*,iostat=ios) winsize
      if (winsize.lt.0) then
         print *,'How can the window size be less than 0?!'
         print *,'Stopping...'
         stop
      endif
   case('Input data folder (no / on the end)')
      dataloc = trim(adjustl(str))
   case('Output folder (no / on the end)')
      otfld = trim(adjustl(str))
      otfldtrm = len_trim(otfld)
   case('Output file basename') 
      otfile = trim(adjustl(str))
      otfiltrm = len_trim(otfile)
      otstr = otfld(1:otfldtrm)//'/'//otfile(1:otfiltrm)
      otstrtrm = len_trim(otstr) 
   case('Use UAVSAR file name convention?')
      uvsryn = trim(adjustl(str))
      if (uvsryn.ne.'yes'.and.uvsryn.ne.'no'.and.uvsryn.ne.'y'.and.uvsryn.ne.'n') then
         print*,'Incorrect input for  Use UAVSAR file name convention?'
         stop
      elseif (uvsryn.eq.'yes'.or.uvsryn.eq.'y') then
         k = scan(dataloc,'/',.true.)
         basenm = dataloc(k+1:k+34)
         trlnm  = dataloc(k+35:k+40)
         HHHHfile = basenm//'HHHH'//trlnm//'.'//dattyp
         VVVVfile = basenm//'VVVV'//trlnm//'.'//dattyp
         HVHVfile = basenm//'HVHV'//trlnm//'.'//dattyp
         HHHVfile = basenm//'HHHV'//trlnm//'.'//dattyp
         HHVVfile = basenm//'HHVV'//trlnm//'.'//dattyp
         HVVVfile = basenm//'HVVV'//trlnm//'.'//dattyp
      endif
   case('*** HHHH filename')
      if (uvsryn.eq.'no'.or.uvsryn.eq.'n') HHHHfile = trim(adjustl(str))
   case('*** VVVV filename')
      if (uvsryn.eq.'no'.or.uvsryn.eq.'n') VVVVfile = trim(adjustl(str))
   case('*** HVHV filename')
      if (uvsryn.eq.'no'.or.uvsryn.eq.'n') HVHVfile = trim(adjustl(str))
   case('*** HHHV filename')
      if (uvsryn.eq.'no'.or.uvsryn.eq.'n') HHHVfile = trim(adjustl(str))
   case('*** HHVV filename')
      if (uvsryn.eq.'no'.or.uvsryn.eq.'n') HHVVfile = trim(adjustl(str))
   case('*** HVVV filename')
      if (uvsryn.eq.'no'.or.uvsryn.eq.'n') HVVVfile = trim(adjustl(str))

   case('Output options')
      rdtype = 2
   case('Output entropy')
      ynstr%ent = trim(adjustl(str))
      if (ynstr%ent.ne.'yes'.and.ynstr%ent.ne.'no') then
         print*,'Incorrect input for ',labell
         print*,'   No output file will be created' 
      elseif(ynstr%ent.eq.'yes') then
         geteig = 1
         nmout = nmout + 1
         open(21,file=otstr(1:otstrtrm)//'.ent',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output anisotropy')
      ynstr%ani = trim(adjustl(str))
      if (ynstr%ani.ne.'yes'.and.ynstr%ani.ne.'no') then
         print*,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif(ynstr%ani.eq.'yes') then
         geteig = 1
         nmout = nmout + 1
         open(22,file=otstr(1:otstrtrm)//'.ani',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output alpha (degrees)')
      ynstr%alp = trim(adjustl(str))
      if (ynstr%alp.ne.'yes'.and.ynstr%alp.ne.'no') then
         print*,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif(ynstr%alp.eq.'yes') then
         geteig = 1
         ynstr%geteigvec = 1
         ynstr%getalp    = 1
         nmout = nmout + 1   
         open(23,file=otstr(1:otstrtrm)//'.alp',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output total Shannon entropy')
      ynstr%shan = trim(adjustl(str))
      if (ynstr%shan.ne.'yes'.and.ynstr%shan.ne.'no') then
         print*,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%shan.eq.'yes') then
         getshan = 1
         nmout = nmout + 1
         open(51,file=otstr(1:otstrtrm)//'.se',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output Shannon polarimetric component')
      ynstr%sep = trim(adjustl(str))
      if (ynstr%sep.ne.'yes'.and.ynstr%sep.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%sep.eq.'yes') then
         getshan = 1
         nmout = nmout + 1
         open(52,file=otstr(1:otstrtrm)//'.sep',access=di,form=fu,status=re,recl=4*width)
      endif 
   case('Output Shannon intensity component')
      ynstr%sei = trim(adjustl(str))
      if (ynstr%sei.ne.'yes'.and.ynstr%sei.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%sei.eq.'yes') then
         getshan = 1
         nmout = nmout + 1
         open(53,file=otstr(1:otstrtrm)//'.sei',access=di,form=fu,status=re,recl=4*width)
      endif 
   case('Output largest eigenvalue  (lambda_1)')
      ynstr%ev1 = trim(adjustl(str))
      if (ynstr%ev1.ne.'yes'.and.ynstr%ev1.ne.'no') then
         print*,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif(ynstr%ev1.eq.'yes') then
         geteig = 1
         nmout = nmout + 1
         open(54,file=otstr(1:otstrtrm)//'.ev1',access=di,form=fu,status=re,recl=4*width)
      endif 
   case('Output median eigenvalue   (lambda_2)')
      ynstr%ev2 = trim(adjustl(str))
      if (ynstr%ev2.ne.'yes'.and.ynstr%ev2.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%ev2.eq.'yes') then
         geteig = 1
         nmout = nmout + 1
         open(55,file=otstr(1:otstrtrm)//'.ev2',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output smallest eigenvalue (lambda_3)')
      ynstr%ev3 = trim(adjustl(str))
      if(ynstr%ev3.ne.'yes'.and.ynstr%ev3.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%ev3.eq.'yes') then
         geteig = 1
         nmout = nmout + 1
         open(56,file=otstr(1:otstrtrm)//'.ev3',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output pseudo probability for lambda_1')
      ynstr%prob1 = trim(adjustl(str))
      if (ynstr%prob1.ne.'yes'.and.ynstr%prob1.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%prob1.eq.'yes') then
         geteig = 1
         nmout = nmout + 1
         open(57,file=otstr(1:otstrtrm)//'.p1',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output pseudo probability for lambda_2')
      ynstr%prob2 = trim(adjustl(str))
      if (ynstr%prob2.ne.'yes'.and.ynstr%prob2.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%prob2.eq.'yes') then
         geteig = 1
         nmout = nmout + 1
         open(58,file=otstr(1:otstrtrm)//'.p2',access=di,form=fu,status=re,recl=4*width)
      endif 
   case('Output pseudo probability for lambda_3')
      ynstr%prob3 = trim(adjustl(str))
      if (ynstr%prob3.ne.'yes'.and.ynstr%prob3.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%prob3.eq.'yes') then
         geteig = 1
         nmout = nmout + 1
         open(59,file=otstr(1:otstrtrm)//'.p3',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output averaged intensity')
      ynstr%inten = trim(adjustl(str))
      if (ynstr%inten.ne.'yes'.and.ynstr%inten.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%inten.eq.'yes') then
         geteig = 1
         nmout = nmout + 1
         open(24,file=otstr(1:otstrtrm)//'.lam',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output trace(T_3)')
      ynstr%trace = trim(adjustl(str))
      if (ynstr%trace.ne.'yes'.and.ynstr%trace.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%trace.eq.'yes') then
         nmout = nmout + 1
         open(61,file=otstr(1:otstrtrm)//'.trt',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output det(T_3)')
      ynstr%detr = trim(adjustl(str))
      if (ynstr%detr.ne.'yes'.and.ynstr%detr.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%detr.eq.'yes') then
         nmout = nmout + 1
         open(62,file=otstr(1:otstrtrm)//'.det',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output I_2(T_3)')
      ynstr%inv2 = trim(adjustl(str))
      if (ynstr%inv2.ne.'yes'.and.ynstr%inv2.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%inv2.eq.'yes') then
         nmout = nmout + 1
         open(63,file=otstr(1:otstrtrm)//'.iv2',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output Pauli_1 (|HH+VV|) value')
      ynstr%paul1 = trim(adjustl(str))
      if (ynstr%paul1.ne.'yes'.and.ynstr%paul1.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%paul1.eq.'yes') then
         nmout = nmout + 1
         open(31,file=otstr(1:otstrtrm)//'.t11',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output Pauli_2 (|HH-VV|) value')
      ynstr%paul2 = trim(adjustl(str))
      if (ynstr%paul2.ne.'yes'.and.ynstr%paul2.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%paul2.eq.'yes') then
         nmout = nmout + 1
         open(32,file=otstr(1:otstrtrm)//'.t22',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output Pauli_3 (2|HV|) value')
      ynstr%paul3 = trim(adjustl(str))
      if (ynstr%paul3.ne.'yes'.and.ynstr%paul3.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%paul3.eq.'yes') then
         nmout = nmout + 1
         open(33,file=otstr(1:otstrtrm)//'.t33',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output Beta (degrees)')
      ynstr%beta = trim(adjustl(str))
      if (ynstr%beta.ne.'yes'.and.ynstr%beta.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%beta.eq.'yes') then
         geteig = 1
         ynstr%geteigvec = 1
         ynstr%getbeta   = 1
         nmout = nmout + 1
         open(25,file=otstr(1:otstrtrm)//'.beta',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output gamma (degrees)')
      ynstr%gama = trim(adjustl(str))
      if (ynstr%gama.ne.'yes'.and.ynstr%gama.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%gama.eq.'yes') then
         geteig = 1
         ynstr%geteigvec = 1
         nmout = nmout + 1
         open(26,file=otstr(1:otstrtrm)//'.gamma',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output delta (degrees)')
      ynstr%delta = trim(adjustl(str))
      if (ynstr%delta.ne.'yes'.and.ynstr%delta.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%delta.eq.'yes') then
         geteig = 1
         ynstr%geteigvec = 1
         nmout = nmout + 1
         open(27,file=otstr(1:otstrtrm)//'.delta',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output phi (degrees)')
      ynstr%phi = trim(adjustl(str))
      if (ynstr%phi.ne.'yes'.and.ynstr%phi.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%phi.eq.'yes') then
         geteig = 1
         ynstr%geteigvec = 1
         nmout = nmout + 1
         open(28,file=otstr(1:otstrtrm)//'.phi',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output alpha_1 (degrees)')
      ynstr%alp1 = trim(adjustl(str))
      if (ynstr%alp1.ne.'yes'.and.ynstr%alp1.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%alp1.eq.'yes') then
         geteig = 1
         ynstr%geteigvec = 1
         ynstr%getalp    = 1
         nmout = nmout + 1
         open(41,file=otstr(1:otstrtrm)//'.alp1',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output alpha_2 (degrees)')
      ynstr%alp2 = trim(adjustl(str))
      if (ynstr%alp2.ne.'yes'.and.ynstr%alp2.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%alp2.eq.'yes') then
         geteig = 1
         ynstr%geteigvec = 1
         ynstr%getalp    = 1
         nmout = nmout + 1
         open(42,file=otstr(1:otstrtrm)//'.alp2',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output alpha_3 (degrees)')
      ynstr%alp3 = trim(adjustl(str))
      if (ynstr%alp3.ne.'yes'.and.ynstr%alp3.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%alp3.eq.'yes') then
         geteig = 1
         ynstr%geteigvec = 1
         ynstr%getalp    = 1
         nmout = nmout + 1
         open(43,file=otstr(1:otstrtrm)//'.alp3',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output Beta_1 (degrees)')
      ynstr%beta1 = trim(adjustl(str))
      if (ynstr%beta1.ne.'yes'.and.ynstr%beta1.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%beta1.eq.'yes') then
         geteig = 1
         ynstr%geteigvec = 1
         ynstr%getbeta   = 1
         nmout = nmout + 1
         open(44,file=otstr(1:otstrtrm)//'.beta1',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output Beta_2 (degrees)')
      ynstr%beta2 = trim(adjustl(str))
      if (ynstr%beta2.ne.'yes'.and.ynstr%beta2.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%beta2.eq.'yes') then
         geteig = 1
         ynstr%geteigvec = 1
         ynstr%getbeta   = 1
         nmout = nmout + 1
         open(45,file=otstr(1:otstrtrm)//'.beta2',access=di,form=fu,status=re,recl=4*width)
      endif
   case('Output Beta_3 (degrees)')
      ynstr%beta3 = trim(adjustl(str))
      if (ynstr%beta3.ne.'yes'.and.ynstr%beta3.ne.'no') then
         print *,'Incorrect input for ',labell
         print*,'   No output file will be created'
      elseif (ynstr%beta3.eq.'yes') then
         geteig = 1
         ynstr%geteigvec = 1
         ynstr%getbeta   = 1
         nmout = nmout + 1
         open(46,file=otstr(1:otstrtrm)//'.beta3',access=di,form=fu,status=re,recl=4*width)
      endif
   case('dB multiplier')
      read(str,*,iostat=ios) parstr%dbmul  
   end select
   endif
   endif 
   enddo 

   hhhhtrm  = len_trim(HHHHfile)
   vvvvtrm  = len_trim(VVVVfile)
   hvhvtrm  = len_trim(HVHVfile)
   hhhvtrm  = len_trim(HHHVfile)
   hhvvtrm  = len_trim(HHVVfile)
   hvvvtrm  = len_trim(HVVVfile)
   dattrm   = len_trim(dataloc)

!  Evaluate output array size
   if (ynstr%eline.ne.1) then
      call stat(dataloc(1:dattrm)//'/'//HHHHfile,stv)
      lines = stv(8)/(4*cols)
   else
      lines = erow
   endif 

   if(erow.gt.lines) erow = lines

   if ((erow-srow+1).lt.lines) then
      lines = erow-srow+1
   else
      erow = lines
   endif


   print *,' '; print *,' '
   print *,'~~~~~~~~          Running UAVSAR          ~~~~~~~~'
   print *,'~~~~~~~~    Polarimetric Decomposition    ~~~~~~~~'                 
   print *,' '; print *,' '
   print *,'Decomposing data in:  ';print*,' '
   print *,'         ',dataloc(1:dattrm)
   print *,' ';print *,' '
   print *,'Input files: ';print *,' '
   print *,'         ',HHHHfile(1:hhhhtrm)
   print *,'         ',VVVVfile(1:vvvvtrm)
   print *,'         ',HVHVfile(1:hvhvtrm)
   print *,'         ',HHHVfile(1:hhhvtrm)
   print *,'         ',HHVVfile(1:hhvvtrm)
   print *,'         ',HVVVfile(1:hvvvtrm)
   print *,' ';print *,' '
   print *,'Writing to:  ';print*,' '
   print *,'         ',otfld(1:otfldtrm)
   print *,' ';print *,' '
   print *,'Total width   = ',width
   print *,'Total lines   = ',lines
   print *,' ';print *,' '
   write(*,110) 'Output files = ',lines*width*4*1.0d-6,' (MB each)'
   if (dble(lines*width)*dble(4*nmout)*1.0d-9.gt.1) then
   write(*,113) 'Total output = ',dble(lines*width)*dble(4*nmout)*1.0d-9,' (GB)'
   else
   write(*,113) 'Total output = ',dble(lines*width)*dble(4*nmout)*1.0d-6,' (MB)'
   endif
   print *,' ';print *,' '
   
!  Define parameter structure components
   parstr%ex    = 1.0d0/3.0d0
   parstr%ex2   = 2.0d0**parstr%ex
   parstr%ex2s  = parstr%ex2**2.0d0
   parstr%sqrt3 = sqrt(3.0d0)
   parstr%nmcp  = cmplx(1.0d0,parstr%sqrt3)
   parstr%nmcn  = cmplx(1.0d0,-parstr%sqrt3)
   parstr%pi    = 4.0d0*atan(1.0d0)
   parstr%r2d   = 180.0d0/parstr%pi
   parstr%log3i = 1/log(3.0d0)
   parstr%exp1  = exp(1.0d0)

!  Set averaging window values   
    if (winsize.eq.1.or.winsize.eq.0) then
      winsizside = 0
      winsizeeff = 1
      winsizinv = 1.0d0

    elseif (mod(winsize,2.0d0).eq.0.and.winsize.ne.0) then
      winsizside =  int(winsize/2)
      winsizeeff =  int(winsize + 1)
      winsizinv = 1.0d0/((winsize+1)*(winsize+1))
      recocount = 0
    else
      winsizeeff = int(winsize)
      winsizside = int((winsize-1)/2)
      winsizinv = 1.0d0/((winsize)*(winsize))
      recocount = 1
    endif
   
!  Open input files
   open(12,file=dataloc(1:dattrm)//'/'//HHHHfile(1:hhhhtrm),access='direct',form=fu,status='old',recl=4*cols)
   open(13,file=dataloc(1:dattrm)//'/'//HVHVfile(1:hvhvtrm),access='direct',form=fu,status='old',recl=4*cols)
   open(14,file=dataloc(1:dattrm)//'/'//VVVVfile(1:vvvvtrm),access='direct',form=fu,status='old',recl=4*cols)
   open(15,file=dataloc(1:dattrm)//'/'//HHHVfile(1:hhhvtrm),access='direct',form=fu,status='old',recl=8*cols)
   open(16,file=dataloc(1:dattrm)//'/'//HHVVfile(1:hhvvtrm),access='direct',form=fu,status='old',recl=8*cols)
   open(17,file=dataloc(1:dattrm)//'/'//HVVVfile(1:hvvvtrm),access='direct',form=fu,status='old',recl=8*cols)


!  Allocate coherency matrix element arrays
   allocate(t11a(width),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in t11a'
   allocate(t22a(width),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in t22a'
   allocate(t33a(width),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in t33a'
   allocate(t12a(width),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in t12a'
   allocate(t13a(width),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in t13a'
   allocate(t23a(width),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in t23a'
   allocate(zrmt(width),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in zrmt'
  

!  Calculate the coherency matrix components for the defined window

   imlines:do iii=1,lines

      allocate(hhhh(cols,winsizeeff),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in hhhh'
      allocate(hvhv(cols,winsizeeff),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in hvhv'
      allocate(vvvv(cols,winsizeeff),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in vvvv'
      allocate(hhhv(cols,winsizeeff),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in hhhv'
      allocate(hhvv(cols,winsizeeff),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in hhvv'
      allocate(hvvv(cols,winsizeeff),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in hvvv'

      if (ynstr%paul1.eq.'yes'.or.ynstr%paul2.eq.'yes'.or.ynstr%paul3.eq.'yes') then
         allocate(paul1(width),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in paul1'
         allocate(paul2(width),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in paul2'
         allocate(paul3(width),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in paul3'
      endif


      if (iii.le.winsizside) then
         reco = 1 + srow - 1
         wlin = iii + winsizside
         pauliln = iii
      elseif (iii.gt.(lines-winsizside)) then
         reco = iii-winsizside + srow - 1
         wlin = winsizeeff - (winsizside - (lines - iii))
         pauliln = winsizside + 1
      else
         reco = iii - winsizside + srow - 1
         wlin = winsizeeff
         pauliln = winsizside + 1
      endif

      winreadcount = 0
      do zzz=reco,(reco+wlin-1)
         winreadcount = winreadcount + 1

         read(12,rec=zzz) hhhh(1:cols,winreadcount)
         read(13,rec=zzz) hvhv(1:cols,winreadcount) 
         read(14,rec=zzz) vvvv(1:cols,winreadcount)
         read(15,rec=zzz) hhhv(1:cols,winreadcount)
         read(16,rec=zzz) hhvv(1:cols,winreadcount)
         read(17,rec=zzz) hvvv(1:cols,winreadcount)
      enddo

      call print_status(iii,lines,2)
      
      colcnt = 0
      zrcnt  = 0

      imcols:do jjj=left,right

         colcnt = colcnt + 1

         if (colcnt.le.winsizside) then
            wcol = colcnt + winsizside
            scol = left
            paulicol = colcnt
         elseif (jjj.gt.(right-winsizside)) then
            wcol = winsizeeff - (winsizside - (right - jjj))
            scol = jjj - winsizside
            paulicol = winsizside + 1
         else
            wcol = winsizeeff
            scol = jjj - winsizside
            paulicol = winsizside + 1
         endif

         if (wcol.gt.winsizeeff) print *,'fix the column windowing ass'
         if (wcol.lt.0) print *,'wcol went negative ass clown'

         t11 = 0.0d0
         t22 = 0.0d0
         t33 = 0.0d0
         t12 = (0.0d0,0.0d0)
         t13 = (0.0d0,0.0d0)
         t23 = (0.0d0,0.0d0)

         do ggg=1,wlin
            do hhh=1,wcol
         
            hhhh_r = hhhh((scol+hhh-1),ggg)
            hvhv_r = hvhv((scol+hhh-1),ggg)
            vvvv_r = vvvv((scol+hhh-1),ggg)
            hhhv_c = hhhv((scol+hhh-1),ggg)
            hhvv_c = hhvv((scol+hhh-1),ggg)
            hvvv_c = hvvv((scol+hhh-1),ggg)

            hhhv_r = real(hhhv_c)
            hhhv_i = aimag(hhhv_c)
            hhvv_r = real(hhvv_c)
            hhvv_i = aimag(hhvv_c)
            hvvv_r = real(hvvv_c)
            hvvv_i = aimag(hvvv_c)

            t11 = t11 + 0.5d0*(hhhh_r + 2.0d0*hhvv_r + vvvv_r)
            t22 = t22 + 0.5d0*(hhhh_r - 2.0d0*hhvv_r + vvvv_r)
            t33 = t33 + 2.0d0*hvhv_r    

            t12_r = 0.5d0*(hhhh_r - vvvv_r)
            t12_i = -hhvv_i
            t13_r = hhhv_r + hvvv_r
            t13_i = hhhv_i - hvvv_i
            t23_r = hhhv_r - hvvv_r
            t23_i = hhhv_i + hvvv_i

            t12 = t12 + cmplx(t12_r,t12_i)
            t13 = t13 + cmplx(t13_r,t13_i)
            t23 = t23 + cmplx(t23_r,t23_i)

            

            if (ggg.eq.pauliln.and.hhh.eq.paulicol) then
               if (ynstr%paul1.eq.'yes'.or.ynstr%paul2.eq.'yes'.or.ynstr%paul3.eq.'yes') then
                  
                  paulval1 = sqrt(hhhh_r + 2.0d0*hhvv_r + vvvv_r)
                  paulval2 = sqrt(hhhh_r - 2.0d0*hhvv_r + vvvv_r)
                  paulval3 = sqrt(hvhv_r)

                  paul1(colcnt) = parstr%dbmul*log10(paulval1)
                  paul2(colcnt) = parstr%dbmul*log10(paulval2)
                  paul3(colcnt) = parstr%dbmul*log10(2.0d0*paulval3)

                  if (paulval1.eq.0.0d0) paul1(colcnt)=0.0d0
                  if (paulval2.eq.0.0d0) paul2(colcnt)=0.0d0
                  if (paulval3.eq.0.0d0) paul3(colcnt)=0.0d0

               endif
            endif

            enddo
         enddo

         if (t11.ne.0.0d0.and.t22.ne.0.0d0.and.t33.ne.0.0d0) then

            if(colcnt.le.winsizside.or.colcnt.gt.(right-winsizside)) then
               t11 = t11/(wlin*wcol)
               t22 = t22/(wlin*wcol)
               t33 = t33/(wlin*wcol)
               t12 = t12/(wlin*wcol)
               t13 = t13/(wlin*wcol)
               t23 = t23/(wlin*wcol)
            elseif(iii.le.winsizside.or.iii.gt.(lines-winsizside)) then
               t11 = t11/(wlin*wcol)
               t22 = t22/(wlin*wcol)
               t33 = t33/(wlin*wcol)
               t12 = t12/(wlin*wcol)
               t13 = t13/(wlin*wcol)
               t23 = t23/(wlin*wcol)
            else
               t11 = t11*winsizinv
               t22 = t22*winsizinv
               t33 = t33*winsizinv
               t12 = t12*winsizinv
               t13 = t13*winsizinv
               t23 = t23*winsizinv
            endif           

         else

            t12 = cmplx(0.0d0,0.0d0)
            t13 = t12
            t23 = t12

            zrcnt = zrcnt + 1
            zrmt(zrcnt) = colcnt 
         endif
     
         t11a(colcnt) = t11
         t22a(colcnt) = t22
         t33a(colcnt) = t33
         t12a(colcnt) = t12
         t13a(colcnt) = t13
         t23a(colcnt) = t23
 
      enddo imcols   ! input column loop



      deallocate(hhhh,hvhv,vvvv,hhhv,hhvv,hvvv)

!  Write out desired parameters:
!        Write Pauli elements 
         if (ynstr%paul1.eq.'yes') write(31,rec=iii) (paul1(mm),mm=1,width)
         if (ynstr%paul2.eq.'yes') write(32,rec=iii) (paul2(mm),mm=1,width)
         if (ynstr%paul3.eq.'yes') write(33,rec=iii) (paul3(mm),mm=1,width)
         if (allocated(paul1))deallocate(paul1,paul2,paul3)
      
      if (geteig.eq.1) then
         !  Call eigenvalue and eigvector subroutines as necessary 
         call eigen(iii,width,t11a,t22a,t33a,t12a,t13a,t23a,zrcnt,zrmt,otstrtrm,otstr,ynstr,parstr)
      elseif(ynstr%detr.eq.'yes'.or.ynstr%trace.eq.'yes'.or.ynstr%inv2.eq.'yes') then
         !  Call wtdi2 if trace, det, or i2 are selected but no eigenvalue or eigenvector
         !             parameters are selected in cmd file 
         call wtdi2(iii,width,t11a,t22a,t33a,t12a,t13a,t23a,otstrtrm,otstr,ynstr)
      endif

      if (getshan.eq.1) then
         !   Call shannon entropy subroutine as necessary
         call shan(iii,width,t11a,t22a,t33a,t12a,t13a,t23a,zrcnt,zrmt,otstrtrm,otstr,ynstr,parstr)
      endif

   enddo imlines      ! input image row loop

   open(7,file=otstr(1:otstrtrm)//'.log',status='unknown')

   call idate(today)
   call itime(now)
   write (7,130) 'Rows in output    = ',lines
   write (7,130) 'Samples in output = ',width
   write (7,140) ' ',' ',' '
   write (7,130) 'Starting row      = ',srow
   write (7,130) 'Ending row        = ',erow
   write (7,130) 'Left column       = ',left
   write (7,130) 'Right column      = ',right
   write (7,140) ' ',' ',' '
   write (7,130) 'Width of input = ',cols
   write (7,140) ' ',' ',' '
   write (7,131) 'window size = ',winsizeeff,' x ',winsizeeff
   write (7,141) 'window size is the number of lines x columns used to calculated the averaged'
   write (7,141) '      coherency (t3) matrix such that <t3> = sum(i,j)[T_(ij)]/n'
   write (7,141) '      where i -> lines, j -> columns and n = i*j'  
   write (7,140) ' ',' ',' '
   write (7,151) 'dB multfactor = ',parstr%dbmul
   write (7,140) ' ',' ',' '
   write (7,141) 'raw data location:  '//dataloc(1:dattrm)
   write (7,141) 'output file name prefix:  '//otstr(1:otstrtrm)
   write (7,140) ' ',' ',' '
   write (7,140) ' ',' ',' '
   write (7,166) today(2),today(1),today(3),now

   print*,' ';print*,' ';print*,'Done!';print*,' ';print*,' '


100   format(a8,i6)
110   format(a15,f6.2,a10)
113   format(a15,f6.2,a5)
166   format('Created (mm/dd/yyyy PST):: ',i2.2,'/'i2.2,'/'i4.4,' @ ',i2.2,':',i2.2,':',i2.2)
150   format(a45,es12.4)
151   format(a45,f12.1)
140   format(a9,a30,a)
141   format(a)
130   format(a20,i6)
131   format(a14,i2,a3,i2)

   if (ynstr%paul1.eq.'yes'.or.ynstr%paul2.eq.'yes'.or.ynstr%paul3.eq.'yes') then
      close(31);close(32);close(33)
   endif
999  close(3);close(12);close(13);close(14);close(15);close(16);close(17)   
   end program PolarDecomp







!  ****************************************************************************
!  ****************************************************************************

   subroutine eigen(i,w,t11a,t22a,t33a,t12a,t13a,t23a,zrcnt,zrmt,ott,otstr,yn,parstr)

!****************************************************************



!****************************************************************
!**     
!**   Routine NAME:  eigen
!**     
!**   DATE WRITTEN: April 2011
!**     
!**   PROGRAMMER: 
!**                        Brent Minchew  
!**               California Institute of Technology 
!**                       bminchew@caltech.edu
!**     
!**   FUNCTIONAL DESCRIPTION: Calculates eigenvalues and 
!**            associated eigenvalue parameters of <T3> 
!**         
!**   DEFAULT OUTPUT:  Writes designated parameters to file        
!**   
!**   OPTIONAL OUTPUTS: Eigenvalues, Eigenvalue probabilities, 
!**                  Entropy, Anistropy, Trace(<T3>), Det(<T3>)
!**                  I2(<T3>)
!**   
!**   ROUTINES CALLED:  eigvect (if any eigenvector component is 
!**                  selected in cmd file)   
!**     
!**   NOTES:  
!**
!**   ASSUMPTIONS:  None 
!**      
!**     
!**   UPDATE LOG:
!**
!**   Date Changed        Reason Changed        CR # and Version #
!**   ------------       ----------------        -----------------
!**     
!*****************************************************************

   
   implicit none

   character*500          otstr
   integer                i,j,w,ott,ierr,mm,zrcnt
   integer                zrmt(w)
   real*4                 t11a(w),t22a(w),t33a(w),dummy   
   complex*8              t12a(w),t13a(w),t23a(w)
   real*4,allocatable::   ev1(:),ev2(:),ev3(:)
   real*8,allocatable::   a(:),dett(:),trt(:)
   real*4,allocatable::   p1(:),p2(:),p3(:)
   real*4,allocatable::   ent(:),ani(:),inten(:),alp(:),beta(:),gam(:),delta(:)
   complex*8,allocatable::   t12c(:),t13c(:),t23c(:)
   complex*16,allocatable::  b(:),c(:)

   type :: parstruc
      real*8                 ex,ex2,ex2s,sqrt3,log3i,pi,r2d,exp1,dbmul
      complex*16             nmcp,nmcn
   end type

   type(parstruc) parstr

   type :: struc 
      character*3         ent,ani,alp,beta,gama,delta,phi,inten
      character*3         ev1,ev2,ev3,shan,sep,sei,trace,detr,inv2
      character*3         prob1,prob2,prob3,paul1,paul2,paul3
      character*3         alp1,alp2,alp3,beta1,beta2,beta3
      integer             geteigvec,getalp,getbeta

   end type

   type(struc) yn

   allocate(ev1(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in ev1'
   allocate(ev2(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in ev2'
   allocate(ev3(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in ev3'

   allocate(a(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in a'
   allocate(b(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in b'
   allocate(dett(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in dett'
   allocate(t12c(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in t12c'
   allocate(t13c(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in t13c'
   allocate(t23c(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in t23c'


   t12c = conjg(t12a)
   t13c = conjg(t13a)
   t23c = conjg(t23a)


   dett = t11a*t22a*t33a - t33a*t12a*t12c - t22a*t13a*t13c + t12a*t13c*t23a + &
      & t12c*t13a*t23c - t11a*t23a*t23c

   a = t11a*t22a + t11a*t33a + t22a*t33a - t12a*t12c - t13a*t13c - t23a*t23c
   b = t11a*t11a - t11a*t22a + t22a*t22a - t11a*t33a - t22a*t33a + t33a*t33a + & 
      & 3.0d0*(t12a*t12c + t13a*t13c + t23a*t23c)
   b = cmplx(real(b),0.0d0)

   deallocate(t12c,t13c,t23c)

   allocate(trt(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in trt'
   allocate(c(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in c'

   trt = t11a + t22a + t33a

   c = 27.0d0*dett - 9.0d0*a*trt + 2.0d0*(trt**3.0d0)  
   c = cmplx(real(c),0.0d0)
   c = c + cdsqrt(c**2.0d0 - 4.0d0*b**3.0d0)
   c = c**parstr%ex

   if (yn%detr.eq.'yes') write(62,rec=i) (sngl(dett(mm)),mm=1,w)
   if (yn%trace.eq.'yes')write(61,rec=i) (sngl(trt(mm)),mm=1,w)
   if (yn%inv2.eq.'yes') write(63,rec=i) (sngl(a(mm)),mm=1,w)

   trt = parstr%ex*trt   
   ev1 = trt + parstr%ex2*b/(3.0d0*c) + c/(3.0d0*parstr%ex2) 
   ev2 = trt - parstr%nmcp*b/(3.0d0*parstr%ex2s*c) - parstr%nmcn*c/(6.0d0*parstr%ex2)
   ev3 = trt - parstr%nmcn*b/(3.0d0*parstr%ex2s*c) - parstr%nmcp*c/(6.0d0*parstr%ex2)

!  remove inf values from eigenvalues
   if (zrcnt.ge.1) then
      ev1(zrmt(1:zrcnt)) = 0.0d0
      ev2(zrmt(1:zrcnt)) = 0.0d0
      ev3(zrmt(1:zrcnt)) = 0.0d0
   endif

!  restore trt (changed above for eigenvalue calcs)
   trt = t11a + t22a + t33a 

   deallocate(a,b,c,dett)

   if (yn%ev1.eq.'yes') write(54,rec=i) (ev1(mm),mm=1,w)

!  ensure ev2 >= ev3 for writing 
   if (yn%ev2.eq.'yes'.or.yn%ev3.eq.'yes'.or.yn%prob2.eq.'yes'.or.yn%prob3.eq.'yes'&
        &.or.yn%getalp.eq.1.or.yn%getbeta.eq.1) then
      do j=1,w
         if (ev2(j).lt.ev3(j)) then
            dummy  = ev2(j)
            ev2(j) = ev3(j)
            ev3(j) = dummy
         endif
      enddo
      if (yn%ev2.eq.'yes') write(55,rec=i) (ev2(mm),mm=1,w)
      if (yn%ev3.eq.'yes') write(56,rec=i) (ev3(mm),mm=1,w)
   endif

   if (yn%geteigvec.eq.1) then
      call eigvect(i,w,t11a,t22a,t33a,t12a,t13a,t23a,ev1,ev2,ev3,trt,zrcnt,zrmt,ott,otstr,yn,parstr)
   endif

   if (yn%ani.eq.'yes') then
      allocate(ani(w),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in ani'

      ani = abs(ev2 - ev3)/(ev2 + ev3)
      if (zrcnt.ge.1) ani(zrmt(1:zrcnt)) = 0.0d0

      write(22,rec=i) (ani(mm),mm=1,w)
      deallocate(ani)
   endif
   
   if (yn%ent.eq.'yes'.or.yn%inten.eq.'yes'.or.yn%prob1.eq.'yes'.or.yn%prob2.eq.'yes'.or.yn%prob3.eq.'yes') then
      allocate(p1(w),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in eigen p1'
      allocate(p2(w),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in eigen p2'
      allocate(p3(w),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in eigen p3'
      
      p1 = ev1/trt
      p2 = ev2/trt
      p3 = ev3/trt

      if (zrcnt.ge.1) then
         p1(zrmt(1:zrcnt)) = 0.0d0
         p2(zrmt(1:zrcnt)) = 0.0d0
         p3(zrmt(1:zrcnt)) = 0.0d0
      endif

      if (yn%ent.eq.'yes') then
         allocate(ent(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in ent'
      
         ent = -(p1*log(p1) + p2*log(p2) + p3*log(p3))*parstr%log3i  
         if (zrcnt.ge.1) ent(zrmt(1:zrcnt)) = 0.0d0

         write(21,rec=i) (ent(mm),mm=1,w)
         deallocate(ent) 
      endif

      if (yn%inten.eq.'yes') then
         allocate(inten(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in inten'

         inten = ev1*p1 + ev2*p2 + ev3*p3
         inten = parstr%dbmul*log10(inten)
         if (zrcnt.ge.1) inten(zrmt(1:zrcnt)) = 0.0d0

         write(24,rec=i) (inten(mm),mm=1,w)
         deallocate(inten)
      endif
      
      if (yn%prob1.eq.'yes') write(57,rec=i) (p1(mm),mm=1,w) 
      if (yn%prob2.eq.'yes') write(58,rec=i) (p2(mm),mm=1,w)
      if (yn%prob3.eq.'yes') write(59,rec=i) (p3(mm),mm=1,w) 

      deallocate(p1,p2,p3)

   endif

   deallocate(ev1,ev2,ev3,trt)

   end subroutine 



!  ****************************************************************************
!  ****************************************************************************

   subroutine eigvect(i,w,t11a,t22a,t33a,t12a,t13a,t23a,ev1,ev2,ev3,trt,&
            &zrcnt,zrmt,ott,otstr,yn,parstr) 

!****************************************************************



!****************************************************************
!**     
!**   Routine NAME:  eigvect
!**     
!**   DATE WRITTEN: April 2011
!**     
!**   PROGRAMMER: 
!**                        Brent Minchew  
!**               California Institute of Technology 
!**                       bminchew@caltech.edu
!**     
!**   FUNCTIONAL DESCRIPTION: Calculates eigenvectors and 
!**            associated eigenvalue parameters of <T3> 
!**         
!**   DEFAULT OUTPUT:  Writes designated parameters to file        
!**   
!**   OPTIONAL OUTPUTS: alpha, beta, phi, beta, gamma (all in degrees 
!**   
!**   ROUTINES CALLED:  none 
!**     
!**   NOTES:  
!**
!**   ASSUMPTIONS:  None 
!**      
!**      
!**   UPDATE LOG:
!** 
!**   Date Changed        Reason Changed        CR # and Version #
!**   ------------       ----------------        -----------------
!**     
!*****************************************************************

   implicit none

   character*500          otstr
   integer                i,j,w,ott,ierr,mm,zrcnt,zrmt(w)
   real*4                 t11a(w),t22a(w),t33a(w),ev1(w),ev2(w),ev3(w)
   real*8                 trt(w)
   complex*8              t12a(w),t13a(w),t23a(w)
   real*4,allocatable::   alp(:),beta(:),gam(:),delta(:),phi(:)
   real*8,allocatable::   nrm(:),p1(:),p2(:),p3(:)
   real*8,allocatable::   alp1(:),alp2(:),alp3(:)
   real*8,allocatable::   beta1(:),beta2(:),beta3(:)
   real*8,allocatable::   del1(:),del2(:),del3(:)
   real*8,allocatable::   phi1(:),phi2(:),phi3(:)
   complex*8,allocatable::    t12c(:),t13c(:),t23c(:)
   complex*16,allocatable::   u11(:),u12(:),u13(:),u21(:),u22(:),u23(:)

   type :: parstruc
      real*8                 ex,ex2,ex2s,sqrt3,log3i,pi,r2d,exp1,dbmul
      complex*16             nmcp,nmcn
   end type

   type(parstruc) parstr

   type :: struc
      character*3         ent,ani,alp,beta,gama,delta,phi,inten
      character*3         ev1,ev2,ev3,shan,sep,sei,trace,detr,inv2
      character*3         prob1,prob2,prob3,paul1,paul2,paul3
      character*3         alp1,alp2,alp3,beta1,beta2,beta3
      integer             geteigvec,getalp,getbeta

   end type

   type(struc) yn


   allocate(t12c(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in eigvec t12c'
   allocate(t13c(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in eigvec t13c'
   allocate(t23c(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in eigvec t23c'

   t12c = conjg(t12a)
   t13c = conjg(t13a)
   t23c = conjg(t23a)

   allocate(u11(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in u11'
   allocate(u21(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in u21'
   allocate(nrm(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in nrm'

   u11 = (ev1-t33a)/t13c + ((ev1-t33a)*t12c+t13c*t23a)*t23c/(((t22a-ev1)*t13c-t12c*t23c)*t13c)
   u21 = ((ev1-t33a)*t12c+t13c*t23a)/((t22a-ev1)*t13c-t12c*t23c)
   nrm = sqrt(u11*conjg(u11) + u21*conjg(u21) + 1.0d0)

   u11 = u11/nrm
   u21 = u21/nrm

   if (zrcnt.ge.1) then
      u11(zrmt(1:zrcnt)) = 0.0d0
      u21(zrmt(1:zrcnt)) = 0.0d0
   endif

   allocate(u12(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in u12'
   allocate(u22(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in u22'

   u12 = (ev2-t33a)/t13c + ((ev2-t33a)*t12c+t13c*t23a)*t23c/(((t22a-ev2)*t13c-t12c*t23c)*t13c)
   u22 = ((ev2-t33a)*t12c+t13c*t23a)/((t22a-ev2)*t13c-t12c*t23c)
   nrm = sqrt(u12*conjg(u12) + u22*conjg(u22) + 1.0d0)

   u12 = u12/nrm
   u22 = u22/nrm
   
   if (zrcnt.ge.1) then
      u12(zrmt(1:zrcnt)) = 0.0d0
      u22(zrmt(1:zrcnt)) = 0.0d0
   endif

   allocate(u13(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in u13'
   allocate(u23(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in u23'

   u13 = (ev3-t33a)/t13c + ((ev3-t33a)*t12c+t13c*t23a)*t23c/(((t22a-ev3)*t13c-t12c*t23c)*t13c)
   u23 = ((ev3-t33a)*t12c+t13c*t23a)/((t22a-ev3)*t13c-t12c*t23c)
   nrm = sqrt(u13*conjg(u13) + u23*conjg(u23) + 1.0d0)

   u13 = u13/nrm
   u23 = u23/nrm

   if (zrcnt.ge.1) then
      u13(zrmt(1:zrcnt)) = 0.0d0
      u23(zrmt(1:zrcnt)) = 0.0d0
   endif

   deallocate(nrm,t12c,t13c,t23c)

   !  Convert eigenvalues to pseudo probabilities    
   allocate(p1(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in eigvec p1'
   allocate(p2(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in eigvec p2'
   allocate(p3(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in eigvec p3'

   p1 = ev1/trt
   p2 = ev2/trt
   p3 = ev3/trt

   if (zrcnt.ge.1) then
      p1(zrmt(1:zrcnt)) = 0.0d0
      p2(zrmt(1:zrcnt)) = 0.0d0
      p3(zrmt(1:zrcnt)) = 0.0d0
   endif


   if (yn%getalp.eq.1.or.yn%getbeta.eq.1) then
      allocate(alp1(w),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in alp1'
      allocate(alp2(w),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in alp2'
      allocate(alp3(w),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in alp3'

      alp1 = acos(abs(u11))*parstr%r2d
      alp2 = acos(abs(u12))*parstr%r2d
      alp3 = acos(abs(u13))*parstr%r2d

      if (yn%alp1.eq.'yes') write(41,rec=i) (sngl(alp1(mm)),mm=1,w)
      if (yn%alp2.eq.'yes') write(42,rec=i) (sngl(alp2(mm)),mm=1,w)
      if (yn%alp3.eq.'yes') write(43,rec=i) (sngl(alp3(mm)),mm=1,w)

      if (yn%alp.eq.'yes') then
         allocate(alp(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in alp'

            alp = alp1*p1 + alp2*p2 + alp3*p3
            write(23,rec=i) (alp(mm),mm=1,w)
            deallocate (alp)
      endif

      if (yn%getbeta.eq.1) then
         allocate(beta1(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in beta1'
         allocate(beta2(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in beta2'
         allocate(beta3(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in beta3'
         allocate(beta(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in beta'

         beta1 = acos(abs(u21)/sin(alp1/parstr%r2d))*parstr%r2d
         beta2 = acos(abs(u22)/sin(alp2/parstr%r2d))*parstr%r2d
         beta3 = acos(abs(u23)/sin(alp3/parstr%r2d))*parstr%r2d

         if (yn%beta1.eq.'yes') write(44,rec=i) (sngl(beta1(mm)),mm=1,w)
         if (yn%beta2.eq.'yes') write(45,rec=i) (sngl(beta2(mm)),mm=1,w)
         if (yn%beta3.eq.'yes') write(46,rec=i) (sngl(beta3(mm)),mm=1,w)

         if (yn%beta.eq.'yes') then
            allocate(beta(w),stat=ierr)
               if(ierr.ne.0)print*,'allocation error in beta'

            beta = p1*beta1 + p2*beta2 + p3*beta3

            write(25,rec=i) (beta(mm),mm=1,w)
            deallocate(beta)
         endif
         deallocate(beta1,beta2,beta3)
      endif
      deallocate(alp1,alp2,alp3)

   endif

   if (yn%delta.eq.'yes'.or.yn%gama.eq.'yes'.or.yn%phi.eq.'yes') then
      allocate(phi1(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in phi1'
      allocate(phi2(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in phi2'
      allocate(phi3(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in phi3'

      phi1 = datan2(aimag(u11),real(u11))*parstr%r2d
      phi2 = datan2(aimag(u12),real(u12))*parstr%r2d 
      phi3 = datan2(aimag(u13),real(u13))*parstr%r2d

      if (yn%phi.eq.'yes') then
         allocate(phi(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in phi'

         phi = p1*phi1 + p2*phi2 + p3*phi3

         write(28,rec=i) (phi(mm),mm=1,w)
         deallocate(phi)
      endif

      if (yn%delta.eq.'yes') then
         allocate(del1(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in delta1'
         allocate(del2(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in delta2'
         allocate(del3(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in delta3'      
         allocate(delta(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in delta'

         del1 = datan2(aimag(u21),real(u21))*parstr%r2d - phi1
         del2 = datan2(aimag(u22),real(u22))*parstr%r2d - phi2
         del3 = datan2(aimag(u23),real(u23))*parstr%r2d - phi3

         delta = p1*del1 + p2*del2 + p3*del3

         write(27,rec=i) (delta(mm),mm=1,w)

         deallocate(del1,del2,del3,delta)
      endif

      if (yn%gama.eq.'yes') then
         allocate(gam(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in gam'

         gam = -(p1*phi1 + p2*phi2 + p3*phi3) 

         write(26,rec=i) (gam(mm),mm=1,w)

         deallocate(gam)
      endif
   
      deallocate(phi1,phi2,phi3)

      endif
         

   deallocate(p1,p2,p3) 
   

   end subroutine




!  ****************************************************************************
!  ****************************************************************************


   subroutine wtdi2(i,w,t11a,t22a,t33a,t12a,t13a,t23a,ott,otstr,yn)

!****************************************************************



!****************************************************************
!**     
!**   Routine NAME:  wtdi2
!**     
!**   DATE WRITTEN: April 2011
!**     
!**   PROGRAMMER: 
!**                        Brent Minchew  
!**               California Institute of Technology 
!**                       bminchew@caltech.edu
!**     
!**   FUNCTIONAL DESCRIPTION: Calculates the trace, determinant, and/or 
!**            the second invariant I2 of <T3> (only called if eigen is not)
!**         
!**   DEFAULT OUTPUT:  Writes designated parameters to file        
!**   
!**   OPTIONAL OUTPUTS: Trace(<T3>), Det(<T3>), I2(<T3>) 
!**   
!**   ROUTINES CALLED:  None 
!**     
!**   NOTES:  
!**
!**   ASSUMPTIONS:  None 
!**      
!**      
!**   UPDATE LOG:
!** 
!**   Date Changed        Reason Changed        CR # and Version #
!**   ------------       ----------------        -----------------
!**     
!*****************************************************************

   implicit none

   character*500          otstr
   integer                i,j,w,ott,ierr,mm
   real*4                 t11a(w),t22a(w),t33a(w)
   complex*8              t12a(w),t13a(w),t23a(w)
   complex*8,allocatable:: t12c(:),t13c(:),t23c(:)
   real*8,allocatable::   trt(:),dett(:),iv2(:)

   type :: struc
      character*3         ent,ani,alp,beta,gama,delta,phi,inten
      character*3         ev1,ev2,ev3,shan,sep,sei,trace,detr,inv2
      character*3         prob1,prob2,prob3,paul1,paul2,paul3
      integer             geteigvec
   end type

   type(struc) yn



   if (yn%detr.eq.'yes'.or.yn%inv2.eq.'yes') then
      allocate(t12c(w),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in t12c'
      allocate(t13c(w),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in t13c'
      allocate(t23c(w),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in t23c'

      t12c = conjg(t12a)
      t13c = conjg(t13a)
      t23c = conjg(t23a)

      if (yn%detr.eq.'yes') then

         allocate(dett(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in dett'
         dett = t11a*t22a*t33a - t33a*t12a*t12c - t22a*t13a*t13c + &
               & t12a*t13c*t23a + t12c*t13a*t23c - t11a*t23a*t23c

         write(62,rec=i) (sngl(dett(mm)),mm=1,w)
         deallocate(dett)
      endif
      if (yn%inv2.eq.'yes') then
         allocate(iv2(w),stat=ierr)
            if(ierr.ne.0)print*,'allocation error in iv2'

         iv2 = t11a*t22a + t11a*t33a + t22a*t33a - t12a*t12c - &
               & t13a*t13c - t23a*t23c

         write(63,rec=i) (sngl(iv2(mm)),mm=1,w)
         deallocate(iv2)
      endif

      deallocate(t12c,t13c,t23c)
   endif

   if(yn%trace.eq.'yes') then

      allocate(trt(w),stat=ierr)
         if(ierr.ne.0)print*,'allocation error in trt'

      trt = t11a + t22a + t33a

      write(61,rec=i) (sngl(trt(mm)),mm=1,w)
      deallocate(trt)
   endif

   end subroutine 





!  ****************************************************************************
!  ****************************************************************************


   subroutine shan(i,w,t11a,t22a,t33a,t12a,t13a,t23a,zrcnt,zrmt,ott,otstr,yn,parstr)

!****************************************************************



!****************************************************************
!**     
!**   Routine NAME:  shan
!**     
!**   DATE WRITTEN: April 2011
!**     
!**   PROGRAMMER: 
!**                        Brent Minchew  
!**               California Institute of Technology 
!**                       bminchew@caltech.edu
!**     
!**   FUNCTIONAL DESCRIPTION: Calculates the Shannon entropy parameters: 
!**            total, polarimetric, and intensity
!**         
!**   DEFAULT OUTPUT:  Writes designated parameters to file        
!**   
!**   OPTIONAL OUTPUTS: Shannon entropy, Shannon polarimetry, Shannon intensity 
!**   
!**   ROUTINES CALLED:  None 
!**     
!**   NOTES:  
!**
!**   ASSUMPTIONS:  None 
!**      
!**      
!**   UPDATE LOG:
!** 
!**   Date Changed        Reason Changed        CR # and Version #
!**   ------------       ----------------        -----------------
!**     
!*****************************************************************

   implicit none

   character*500          otstr
   integer                i,j,w,ott,ierr,mm,zrcnt,zrmt(w)
   real*4                 t11a(w),t22a(w),t33a(w)
   complex*8              t12a(w),t13a(w),t23a(w)
   complex*8,allocatable:: t12c(:),t13c(:),t23c(:)
   real*4,allocatable::   se(:),sei(:),sep(:)
   real*8,allocatable::   trt(:),dett(:)

   type :: parstruc
      real*8                 ex,ex2,ex2s,sqrt3,log3i,pi,r2d,exp1,dbmul
      complex*16             nmcp,nmcn
   end type

   type(parstruc) parstr

   type :: struc
      character*3         ent,ani,alp,beta,gama,delta,phi,inten
      character*3         ev1,ev2,ev3,shan,sep,sei,trace,detr,inv2
      character*3         prob1,prob2,prob3,paul1,paul2,paul3
      integer             geteigvec
   end type

   type(struc) yn


   allocate(t12c(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in t12c in shan'
   allocate(t13c(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in t13c in shan'
   allocate(t23c(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in t23c in shan'

   t12c = conjg(t12a)
   t13c = conjg(t13a)
   t23c = conjg(t23a)

   allocate(dett(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in dett in shan'

   dett = t11a*t22a*t33a - t33a*t12a*t12c - t22a*t13a*t13c + &
         & t12a*t13c*t23a + t12c*t13a*t23c - t11a*t23a*t23c   

   deallocate(t12c,t13c,t23c)

   allocate(trt(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in trt in shan'

   trt = t11a + t22a + t33a

   allocate(se(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in se in shan'
   allocate(sei(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in sei in shan'
   allocate(sep(w),stat=ierr)
      if(ierr.ne.0)print*,'allocation error in sep in shan'

   sei = 3.0d0*log(parstr%pi*parstr%exp1*parstr%ex*trt)
   sep = log(27.0d0*dett/trt**3.0d0)

   if (zrcnt.ge.1) then
      sei(zrmt(1:zrcnt)) = 0.0d0
      sep(zrmt(1:zrcnt)) = 0.0d0
   endif
      
   se  = sei + sep

   if (yn%shan.eq.'yes') write(51,rec=i) (se(mm),mm=1,w)
   if (yn%sep.eq.'yes') write(52,rec=i) (sep(mm),mm=1,w)
   if (yn%sei.eq.'yes') write(53,rec=i) (sei(mm),mm=1,w)

   deallocate(trt,dett,se,sei,sep)

   end subroutine 

!  ==============================================================================
   subroutine print_status(i,j,s)
   implicit none
   integer        i,j,s
   if(mod((100*i)/j,s).eq.0.or.i.eq.j.or.i.eq.1) then
      write(*,FMT="(A1,A,t21,I4,A)",ADVANCE="NO") achar(13), &
         & "  Percent Complete: ", (100*i)/j, " %"
   endif
   end subroutine

