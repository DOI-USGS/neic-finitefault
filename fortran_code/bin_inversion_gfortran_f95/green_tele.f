C	This program is being changed    
C
C
C     SYNHAR.FOR ----- 1996.6   (  DT=0.1  DTS=0.2  )
C     SYNM6.FOR ------- SYNSH.FOR    1994.4
C
C
C
C-----------------------------------------------------------------
c
       implicit none
       integer linst, ldim, ldim2, lgreen, ltde, inpt, inptd, inpth
       integer nll
       PARAMETER ( LINST=2048,LDIM=2048,LDIM2=2*LDIM,LGREEN=64*LDIM2
     *             ,LTDE=320000)
       PARAMETER(INPT=2050,INPTD=4100,INPTH=1025,NLL=45)
       real pi, pi2, pi4
       integer mmsou, max_seg, nnxs, nnys, nnxy, nnpx, nnpy, nnpxy,
     *  nxs_sub, nys_sub, nxys
       PARAMETER(MMSOU=25,Max_seg=11)
       PARAMETER(NNXS=50,NNYS=20,NNXY=NNXS*NNYS,nnpx=25,nnpy=25)
       parameter(nnpxy=nnpx*nnpy)
       real dd, aa, dip_sub, stk_sub, delay_sub
       DIMENSION DD(Max_seg,NNXs,nnYs),AA(Max_seg,NNXs,nnYs)
       dimension dip_sub(Max_seg),stk_sub(Max_seg)
       dimension nxs_sub(Max_seg),nys_sub(Max_seg),delay_sub(Max_seg)
       dimension nxys(max_seg)
       CHARACTER*200 directory,filterfile,risetimefile,pointsourcefile
       CHARACTER*200 velmodelfile,channelsfiles,waveformsfile,responsefile

c     
       real afa, bta, ro, th, qa, qb, vpp, vs, den, thick
       integer jo

       COMMON/MODEL/AFA(NLL),BTA(NLL),RO(NLL),TH(NLL),QA(NLL),QB(NLL)
       COMMON/SMODEL/Jo,vpp(NLL),vs(NLL),den(NLL),THICK(NLL)
       real cr, ci, obse
       DIMENSION CR(LDIM),CI(LDIM),OBSE(6000)
       real rk, alpha, beta, rho, depth, eta, etb
       COMMON/PARA/RK,ALPHA,BETA,RHO,DEPTH,ETA,ETB,LOVE
       CHARACTER*12 FRTDK
       real dt, df, aw
       integer lnpt, nsyn, jf
       COMMON/PARAM/DT,LNPT,NSYN,DF,JF,aw,FRTDK
       COMPLEX CFW1,CFV1,CFW2,CFV2,FINST,SOURCE,CFW0
       complex ws(ldim),ww(ldim),z
       complex wsyn(ldim),z0
       integer msou, nsou
       COMMON/FIGG1/CFW0(INPTH),CFW1(INPTH),FINST(INPTH)
       COMMON/FIGG2/CFW2(INPTH),CFV1(INPTH),CFV2(INPTH)
       COMMON/FIGSS/MSOU,NSOU,SOURCE(inpt,MMSOU)
       real fau_mod, e_del, e_dep, e_dis, zdepth
       real tdel, pmd, pmu, svmu
       integer idts, nos
       COMMON/SOURECTYPE/IDTS
       DIMENSION fau_mod(Max_seg,NNXY,nnpxy,7)
       DIMENSION TDEL(NNXS,NNYS,NNPX,NNPY)
       DIMENSION PMD(3),PMU(3),SVMU(3),NOS(200)

c     fau_mod  1:x, 2:y, 3:dis
      CHARACTER*200 FNAME4,FNAME5
      
C      IDTS=6 FOR MOMENT SOURCE
C      IDTS=3 FOR DISLOCATION SOURCE	

c
c general purpose constants
c
      real aaz, az, azz, cos2aaz, cos2d, cosaaz, cosal, cosd, 
     * costhea, deg_arc, r, oo, d, d2, delt, delta_az, df_final, dp,
     * dto, dxs, dys, dz, e_az, ameve, c_depth, niu,M1,M2,M3,M4,M5,M0,
     * M6, w, omax, sin2aaz, sin2d, sinaaz, sinal, sind,
     * sinthea, all, thea, time, time0, tmax, tmin, v_min, vp, x,
     * xmean, gdp, rpz, hcru, hsou
      real low_freq, high_freq, lat_sta(200), lon_sta(200)
      real fmax, eang, rang, ttvl, lat_p, lon_p, dis, az_s, baz_s
      real ta0, dta, high_freq2
      integer love, mmm, iud, idata, llove
      DIMENSION EANG(200),HCRU(200),MMM(200),HSOU(200)
      DIMENSION IUD(200),IDATA(200)
      character*6 stname(200),sttyp(200)
      dimension rang(200),az(200),TTVL(200),LLOVE(200)
      character*5 EARTH(200)
      character*40 eventname, string
      character*14 fname(200)
      integer i, i_seg, iflag, ll, kxy, kpxy, al, k, nstb,
     *  io_seg, io_v_d, ir, is, iso, ix, ixs, iy, iys, j, no,
     *  nstaon, nnn, nls, jf_final, leng, lnpt_use, mm, n_seg, nb,
     *  nbs, nl, nla, nlay_so, nx_p, ny_p, theai, nxy, nlay_sz,
     *  nsyn_final, nxs0, nys0, nxyp

      include 'table.h'
 21   format(a40)
c 203  format(I2,1X,A5,1X,A4,1X,A6,1X,A13,1X,F5.2,1X,F6.2,1X,F5.2,1X
c     *   ,F3.1,1X,I1,1X,2(F4.1,1X),I1,1X,I2,1X,I1)
C
      write(*,'(/A/)')'PROGRAM TO COMPUTE TELESEISMIC BODY WAVE GF'
      pi=4.0*atan(1.0)
      pi2 = 2.0 * pi
      pi4 = 4.0 * pi
      ALPHA=6.2
      BETA=3.4
      RHO=2.7
      deg_arc=pi/180.0
      Z0=CMPLX(0.0,0.0)

c get the directory files should be in 
      call getarg(1, directory)       
      filterfile=TRIM(directory)//'filter_tele.txt'
      open(1,file=filterfile)
      read(1,*)string, low_freq, high_freq
      read(1,*)string, dt
      close(1)
      
c     
      write(*,*)'Get fault segments data...'
      risetimefile=TRIM(directory)//'fault&rise_time.txt'
      OPEN(22,FILE=risetimefile)
      READ(22,*)NXS0,NYS0,c_depth,lnpt_use
      READ(22,*)n_seg,DXS,DYS,nx_p,ny_p
      READ(22,*)TA0,DTA,MSOU,v_min,io_v_d
      do i_seg=1,n_seg
         read(22,*)io_seg, dip_sub(i_seg),stk_sub(i_seg)
         read(22,*)nxs_sub(i_seg),nys_sub(i_seg),delay_sub(i_seg)
         do iys=1,nys_sub(i_seg)
            do ixs=1,nxs_sub(i_seg)
               read(22,*)dd(i_seg,ixs,iys),aa(i_seg,ixs,iys)
            enddo
         enddo
         nxys(i_seg)=nxy
         nnn=nnn+nxy
      enddo
      close(22)
      pointsourcefile=TRIM(directory)//'point_sources.txt'
      open(22,file=pointsourcefile)
      do i_seg=1,n_seg
         read(22,*)io_seg
         kxy=0
         do IYS=1,NYS_sub(i_seg)
            do IXS=1,NXS_sub(i_seg)
               kxy=kxy+1
               kpxy=0
               do IY=1,ny_p
                  do IX=1,nx_p
                     kpxy=kpxy+1
                     read(22,*)(fau_mod(i_seg,kxy,kpxy,k),k=1,7)
                  enddo
               enddo
            enddo
         enddo
      enddo
      close(22)

      write(*,*)'Get vel_model...'
      velmodelfile=TRIM(directory)//'vel_model.txt'
      open(15,file=velmodelfile,status='old')
      lnpt=lnpt_use!lnpt_use - 1
      aw=0.01
      read(15,*) jo
      do j=1,jo
         read (15,*) vpp(J),vs(J),den(J),thick(J)
      enddo
      nsyn=2**lnpt
      jf=1+nsyn/2
      df=1.0/(nsyn*dt)!df=1./((nsyn-1)*dt)

!      jf_final=nsyn+1
      jf_final = jf
      nsyn_final=2**lnpt_use
      df_final=1.0/(nsyn_final*dt)!df_final=1.0/(nsyn_final-1)/dt

      depth = c_depth

      call cmodel(nl,depth,nbs)
      close(15)
c
c       make the rise time function
c
      do iso=1,msou
         do i=1,nsyn_final
            cr(i)=0.0
            ci(i)=0.0
         enddo
         r=ta0+(iso-1)*dta
         nls=int(r/dt+0.5)+1
         if(io_v_d.eq.0)then
            do i=1,nls
               cr(i)=(1.0-cos(2.*pi*(i-1)*dt/r))/r
            enddo
         else
            do i=1,nls
               cr(i)=2.*pi*sin(2.*pi*(i-1)*dt/r)/r/r
            enddo
         endif
         call fft(cr,ci,lnpt_use,-1.)
         do j=1,jf
            source(j,iso)=cmplx(cr(j),ci(j))
         enddo
      enddo
c
c       End of Rise Time
c
      channelsfiles=TRIM(directory)//'channels_body.txt'
      waveformsfile=TRIM(directory)//'waveforms_body.txt'
      OPEN(9,FILE=channelsfiles,STATUS='OLD')
      OPEN(13,FILE=waveformsfile,STATUS='UNKNOWN')
C
C    1993.6  SYNM6.FOR, SYNHAR.FOR(1996.6), SYNWWW.FOR(1997.1)  BY YAO
C    M1=M33, M2=M13, M3=M23, M4=M11, M5=M22, M6=M12
C
      READ(9,*) DELT,AL,THEAI,M1,M2,M3,M4,M5,M6,M0
      READ(9,*) IDTS
      READ(9,21) EVENTNAME
      READ(9,*) NSTAON
      NSTB=0
      AMEVE=0.0
      TMAX=0.0
      TMIN=1000.0
      IFLAG=0
c	
c	Here, for every station, we define the green functions for every
c	subfault, plus the time delay of the station and subfault.
c	      
      write(*,*)'Begin to compute response for each channel...'
      DO IR=1,NSTAON
         read(9,*)NOS(IR),EARTH(IR),STTYP(IR),STNAME(IR),FNAME(IR),RANG
     * (IR),AZ(IR),lat_sta(ir),lon_sta(ir),EANG(IR),TTVL(IR),MMM(IR),
     * HCRU(IR),HSOU(IR),LLOVE(IR),IUD(IR),IDATA(IR)
C     
         IF(IDATA(IR).gt.0.OR.MMM(IR).EQ.3) cycle !GOTO 200
         NSTB=NSTB+1
         MM=MMM(IR)
         LOVE=LLOVE(IR)
C
C   MM=0 FOR FAR; MM=1 FOR UPPER; MM=3 FOR PNL
C   Because the response of upper mantle and Pnl is not enough clear
C   to be used to study the detail of earthquake, Such data will not use
C
         if (love.eq.0) then
            write(*,*)'Compute P response for station ', STNAME(IR)
            high_freq2 = high_freq
            fname4 = TRIM(directory)//trim(stname(ir))//'.GRE'
            fname5 = TRIM(directory)//trim(stname(ir))//'.TDE'
         else
            write(*,*)'Compute SH response for station ', STNAME(IR)
            high_freq2 = high_freq / 2.0
            fname4 = TRIM(directory)//trim(stname(ir))//'SH.GRE'
            fname5 = TRIM(directory)//trim(stname(ir))//'SH.TDE'
         endif
         OPEN(12,FILE=FNAME4,STATUS='UNKNOWN',ACCESS='DIRECT',
     *     RECL=LGREEN)
         OPEN(32,FILE=FNAME5,STATUS='UNKNOWN',ACCESS='DIRECT',
     *     RECL=LTDE)
      
         OMAX=0.0  
         read(13,*)
         read(13,*)
         read(13,*)string, dto
         read(13,*)string, no
         read(13,*)string, leng
         read(13,*)   
         read(13,*)(OBSE(I),I=1,NO)
         DO I=1,1+LENG
            X=ABS(OBSE(I))
            IF(X.GT.OMAX)OMAX=X
         enddo

         DP=1.0

         IF(MM.EQ.0) THEN
         
!            WRITE(*,313) STNAME(IR),FNAME4
 313        FORMAT(1X,A6,1X,'  FAR ',A14,A14)
            IF(LOVE.GT.0) THEN
               CALL INTERP(DIST,DTDPSH,61,RANG(IR),RK)
               RK=RK/111.1
            ELSE
               if((rang(ir).ge.30).and.(rang(ir).le.90))then
                  CALL INTERP(DIST,DTDP,61,RANG(IR),RK)
                  RK=RK/111.1
               endif
               if(rang(ir).ge.120)then
                  call interp(dist_pkp,dtdp_df,61,rang(ir),rk)
                  rk=rk/111.1
               endif
            END IF
            CALL INTEDP(RANG(IR),DP,MM,c_depth)
            DP=DP*1.0E-5
         END IF


c	 IF(MM.EQ.1) THEN
c
c	   WRITE(*,49) STNAME(IR),EARTH(IR),FNAME4
c49	   FORMAT(1X,A4,1X,A5,1X,A12,A12)
c	   OPEN(8,FILE=FNAME3,STATUS='UNKNOWN')
c	   READ(8,*)NUPP,DTUPP,R,RK
c	   READ(8,*)(UPPER(J),J=1,NUPP)
c           CLOSE(8)
c           IF(ABS(DT-DTUPP).GT.0.1*DT) THEN
c	     WRITE(*,*)' DT NE DTUPP ; STOP ',DT,DTUPP
c             STOP
c	   END IF
c	   IF(NUPP.LT.NSYN) THEN
c	     DO 15 I=NUPP+1,NSYN
c 15          UPPER(I)=UPPER(NUPP)
c           END IF
c           DO 16 I=1,NSYN
c	   CR(I)=UPPER(I)
c 16        CI(I)=0.0
c	   CALL FFT(CR,CI,LNPT_use,-1.)
c
c         END IF

c
c Seismic source?
c
         CALL ABOTM(PMD,PMU,SVMU,RPZ)
         GDP=DP*RPZ
         iflag=0
c	
c	Instrumental response. We give all the green functions the same
c	instrumetnal response.
c	
         CALL INTBHZ(FINST,TTVL(IR),GDP,IR,IFLAG,MM,directory)
c	 IFALG=1
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++         
         DO I=1,NSYN
            WSYN(I)=Z0
         enddo
         LL=0

         IF(MM.EQ.1) THEN
            DO J=1,JF
               FINST(J)=FINST(J)*(CMPLX(CR(J),CI(J))*DT)
            enddo
         endif
C     
         ll=0
         do i_seg=1,n_seg
c	
c	We define the radiation coefficients.
c	
            D=dip_sub(i_seg)*deg_arc
            COSD=COS(D)
            SIND=SIN(D)
            D2=2.0*D
            COS2D=COS(D2)
            SIN2D=SIN(D2)
            AAZ=(AZ(IR)-stk_sub(i_seg))*deg_arc
            COSAAZ=COS(AAZ)
            SINAAZ=SIN(AAZ)
            AAZ=2.0*AAZ
            COS2AAZ=COS(AAZ)
            SIN2AAZ=SIN(AAZ)
            THEA=stk_sub(i_seg)*deg_arc
            SINTHEA=SIN(THEA)
            COSTHEA=COS(THEA)
            
            DO IYS=1,NYS_sub(i_seg)
               do iy=1,ny_p
                  kxy=(iys-1)*nxs_sub(i_seg)+1
                  kpxy=(iy-1)*nx_p+1
                  zdepth=fau_mod(i_seg,kxy,kpxy,3)
                  if(zdepth .le.0.1)then
                     write(*,*)iys,iy, i_seg
                     write(*,*)"ohoh, the subfault is in air"
                     stop
                  endif
                  CALL RTDK(LOVE,RK,ZDEPTH,NBS) 
c		
c		shear modulous
c		
                  niu=bta(nbs)*bta(nbs)*ro(nbs)*(1.e+10)
                  niu=niu/(pi4*ro(nbs))

                  DO I=1,JF_final
c	a trick to remove the dc offset
                     z=cmplx(0.0,-pi2*df_final*(i-1)*20.0)
                     Z=FINST(I)*niu*cexp(z)
c	
c	Green functions in frequency domain, with added instrumental response
c	
                     IF(LOVE.GT.0) THEN
                        WW(I)=(-COS2D*COSAAZ*CFV1(I)-  
     &                       SIN2D*SIN2AAZ*CFV2(I)*0.5)*Z
                        WS(I)=(-COSD*SINAAZ*CFV1(I)+   
     &                       SIND*COS2AAZ*CFV2(I))*Z
                     ELSE
                        WW(I)=(0.5*SIN2D*CFW0(I)-COS2D*SINAAZ*CFW1(I)+ 
     &                       0.5*SIN2D*COS2AAZ*CFW2(I))*Z
                        WS(I)=(COSD*COSAAZ*CFW1(I) 
     &                      +SIND*SIN2AAZ*CFW2(I))*Z
                     END IF
                  enddo
c
                  DO I=1,JF_final
                     CR(I)=REAL(ww(I)) ! CR(I)=REAL(ww(I))
                     CI(I)=AIMAG(ww(I)) ! CI(I)=AIMAG(ww(I))
                  enddo
c	
c	Back to time domain
c
   
                  CALL REALTR(CR,CI,LNPT_use)
                  CALL FFT(CR,CI,LNPT_use,1.)
                  nb=int(20.0/dt)-1
                  xmean=0.0
                  do i=1,nb
                     xmean=xmean+cr(i)
                  enddo
                  xmean=xmean/nb
               
                  do i=1,nsyn_final
                     cr(i)=cr(i)-xmean
                  enddo
                  do i=1,ldim
                     ci(i)=0.0
                     if(i.gt.nsyn)cr(i)=0.0
                  enddo
c	
c	Back to frequency domain.
c
                  call bandpassfilter(cr,ldim,dt,2,2, 
     &            low_freq,10.0)
                  call bandpassfilter(cr,ldim,dt,2,2, 
     &            0.0,high_freq2)
                  CALL FFT(CR,CI,LNPT_use,-1.)
                  do i=1,jf_final
                     w=pi2*df_final*(i-1)*20.0
                     Z=cmplx(0.0,w)
                     ww(i)=cmplx(cr(i),ci(i))*cexp(z)
                  enddo
                  LL=LL+1
                  WRITE(12,REC=LL)i_seg,iys,iy,JF_final,DT,DF_final, 
     &                     (WW(I),I=1,JF_final)
     
                  DO I=1,JF_final
                     CR(I)=REAL(ws(I))  !CR(I)=REAL(ws(I))
                     CI(I)=AIMAG(ws(I))  !CI(I)=AIMAG(ws(I))
                  enddo
                  CALL REALTR(CR,CI,LNPT_use)
                  CALL FFT(CR,CI,LNPT_use,1.)
                  xmean=0.0
                  do i=1,nb
                     xmean=xmean+cr(i)
                  enddo
                  xmean=xmean/nb
                 
                  do i=1,nsyn_final
                     cr(i)=cr(i)-xmean
                  enddo
                  do i=1,ldim
                     ci(i)=0.0
                     if(i.gt.nsyn)cr(i)=0.0
                  enddo
                  call bandpassfilter(cr,ldim,dt,2,2,  
     &            low_freq,10.0)
                  call bandpassfilter(cr,ldim,dt,2,2, 
     &            0.0,high_freq2)
                  CALL FFT(CR,CI,LNPT_use,-1.)
                  do i=1,jf_final
                     w=pi2*df_final*(i-1)*20.0
                     Z=cmplx(0.0,w)
                     ws(i)=cmplx(cr(i),ci(i))*cexp(z)
                  enddo
                  LL=LL+1
!                  write(*,*)'parameters of GF', jf_final,lnpt_use, 
!     &            nsyn_final
!                  write(*,*)'nonzero values of GF', jf,lnpt,nsyn
                  WRITE(12,REC=LL)i_seg,iys,iy,JF_final,LNPT_use, 
     &                        NSYN_final,(WS(I),I=1,JF_final)
c
c              write(*,*)'NSYN=',nsyn,ll,iys,iy,zdepth
c 	
                  IF(LOVE.EQ.0)VP=AFA(NBS)
                  IF(LOVE.GT.0)VP=BTA(NBS)
                  ETA=SQRT(ABS(RK*RK-1./VP/VP))

c         the time delay caused by shift in vertical position
c
c         TIME=(DEPTH-zdepth)*ETA
c	  Modify this by the layer model
c

                  do nla=1,jo+1
                     if((th(nla).le.depth) 
     &               .and.(th(nla+1).gt.depth)) nlay_so=nla
                     if((th(nla).le.zdepth) 
     &               .and.(th(nla+1).gt.zdepth)) then
                        nlay_sz=nla
                     endif
                  enddo
c
                  if(nlay_sz.eq.nlay_so)then
                     if(llove(ir).ne.2)then
                        eta=sqrt(1/afa(nlay_sz)/afa(nlay_sz)-rk*rk)
                        TIME=-(ZDEPTH-DEPTH)*ETA
                     else
                        eta=sqrt(1/bta(nlay_sz)/bta(nlay_sz)-rk*rk)
                        TIME=-(ZDEPTH-DEPTH)*ETA
                     endif
                  endif
c	
                  if(nlay_sz.lt.nlay_so)then
c
c     delta_t > 0
c
                     time=0.0
                     dz=depth-th(nlay_so)
                     if(llove(ir).ne.2)then
                        eta=sqrt(1./afa(nlay_so)/afa(nlay_so)-rk*rk)
                     else
                        eta=sqrt(1./bta(nlay_so)/bta(nlay_so)-rk*rk)
                     endif
                     time=time+dz*eta

                     do nla=nlay_so,nlay_sz+2,-1
                        dz=th(nla)-th(nla-1)
                        if(llove(ir).ne.2)then
                           eta=sqrt(1./afa(nla-1)/afa(nla-1)-rk*rk)
                        else
                           eta=sqrt(1./bta(nla-1)/bta(nla-1)-rk*rk)
                        endif
                        time=time+dz*eta
                     enddo
                     dz=th(nlay_sz+1)-zdepth
                     if(llove(ir).ne.2)then
                        eta=sqrt(1./afa(nlay_sz)/afa(nlay_sz)-rk*rk)
                     else
                        eta=sqrt(1./bta(nlay_sz)/bta(nlay_sz)-rk*rk)
                     endif
                     time=time+dz*eta
                  endif

                  if(nlay_sz.gt.nlay_so)then
c     delta_t < 0
                     time=0.0
                     dz=depth-th(nlay_so+1)
            
                     if(llove(ir).ne.2)then
                        eta=sqrt(1./afa(nlay_so)/afa(nlay_so)-rk*rk)
                     else
                        eta=sqrt(1./bta(nlay_so)/bta(nlay_so)-rk*rk)
                     endif
                     time=time+dz*eta

                     do nla=nlay_so+1,nlay_sz-1
           
                        dz=-th(nla+1)+th(nla)
                        if(llove(ir).ne.2)then
                           eta=sqrt(1./afa(nla)/afa(nla)-rk*rk)
                        else
                           eta=sqrt(1./bta(nla)/bta(nla)-rk*rk)
                        endif
                        time=time+dz*eta
                     enddo
                     dz=th(nlay_sz)-zdepth
                     if(llove(ir).ne.2)then
                        eta=sqrt(1./afa(nlay_sz)/afa(nlay_sz)-rk*rk)
                     else
                        eta=sqrt(1./bta(nlay_sz)/bta(nlay_sz)-rk*rk)
                     endif
                     time=time+dz*eta            
                  endif
                  time0=time
                  DO IXS=1,NXS_sub(i_seg)
c                     ALL=AA(i_seg,IXS,IYS)*0.017453
                     all = aa(i_seg,ixs,iys)*pi/180.0
                     COSAL=COS(ALL) 
                     SINAL=SIN(ALL)
                     do ix=1,nx_p
                        kxy=(iys-1)*nxs_sub(i_seg)+ixs
                        kpxy=(iy-1)*nx_p+ix
                        e_dep=fau_mod(i_seg,kxy,kpxy,3)
                        e_del=delay_sub(i_seg) 
     &                  +fau_mod(i_seg,kxy,kpxy,4)/v_min
                        e_dis=fau_mod(i_seg,kxy,kpxy,6)
                        e_az =fau_mod(i_seg,kxy,kpxy,7)
c                        delta_az=(e_az-az(ir)*0.017453
                        delta_az=(e_az-az(ir))*pi/180.0
                        time =time0-rk*e_dis*cos(delta_az)  

                        TDEL(IXS,IYS,IX,IY)=TIME     
                        time=time+e_del
                        if(Time.lt.-1.e-2)then
                           write(*,*)'ohoh, later event comes first'
                           write(*,*)ixs,iys,ix,iy,time
                        endif
                        IS=1
                        IF(IS.LT.1)IS=1
                        IF(IS.GT.MSOU)IS=MSOU
                        DO I=1,JF_final
                            W=-PI2*(I-1)*DF_final*TIME
                            z=cmplx(0.0,w)
                            z=cexp(z)*source(i,is)*dt
                            WSYN(I)=WSYN(I) 
     &                +(SINAL*WW(I)+COSAL*WS(I))*DD(i_seg,IXS,IYS)*Z
                        enddo
                     enddo
                  enddo
               enddo
            enddo
c	
c	Writing the time delay for the station and for every point source
c	in a binary file.
c	
            NXY=NXS_sub(I_seg)*NYS_sub(i_seg)
            NXYP=NXY*NX_P*NY_P
!            write(*,*)nx_p,ny_p,nxs_sub(i_seg),nys_sub(i_seg) 
            WRITE(32,REC=i_seg)nxs_sub(i_seg),nys_sub(i_seg), 
     &          nx_p,ny_p, ((((TDEL(IXS,IYS,IX,IY),Ix=1,Nx_P), 
     &         Iy=1,Ny_P),IxS=1,NxS_sub(i_seg)),IyS=1,NyS_sub(i_seg))
         enddo
         CLOSE(12)
         CLOSE(32)
         DO I=1,JF_final
            CR(I)=REAL(WSYN(I))
            CI(I)=AIMAG(WSYN(I))
         enddo
         CALL REALTR(CR,CI,LNPT_use)
         CALL FFT(CR,CI,LNPT_use,1.)
         NL=NSYN_final
         FMAX=0.0
         DO I=1,1+LENG
            X=ABS(CR(I))
            IF(X.GT.FMAX)FMAX=X
         enddo
         OO=OMAX/FMAX
         AMEVE=AMEVE+OO
!         write(*,*)'Useful for debugging: '
!         write(*,*)'station', ir
!         write(*,*)'maximum value of observed trace', omax
!         write(*,*)'maximum value of preliminary synthetic trace', fmax
!         write(*,*)'ratio max_observed/max_synthetic',oo 
C
c200   CONTINUE
      enddo
      AMEVE=AMEVE/NSTAON
      AMEVE=AMEVE*3.39*1.0E16
!      WRITE(*,*)' M0= ',AMEVE
      CLOSE(9)
      CLOSE(14)
      CLOSE(10)
c      write(*,*)NEW_LINE('A')
      write(*,'(/A/)')'END PROGRAM TO COMPUTE TELESEISMIC BODY WAVE GF'
c      write(*,*)NEW_LINE('A')
      STOP
      END



        SUBROUTINE ABOTM(P,P2,SV,RPZ)
        implicit none
        real p, p2, sv, rpz
        real rk, alpha, beta, rho, depth, a, b
        integer love, idts
        COMMON/PARA/RK,ALPHA,BETA,RHO,DEPTH,A,B,LOVE
        COMMON/SOURECTYPE/IDTS
C      IDTS=6 FOR MOMENT SOURCE
C      IDTS=3 FOR DISLOCATION SOURCE	
        DIMENSION P(3),P2(3),SV(3)
c
c general purpose variables
c
        integer i
        real ab, omega, omega2, rk2, rka, rka2, rkb, rkb2, rpp, rsp
        real y, y1
C  FOR SH WAVES, GOTO 30
        RK2 = RK*RK
        RKA=1./ALPHA
        RKB=1./BETA
        RKA2 = RKA*RKA
        RKB2 = RKB*RKB
C         WRITE(*,*)'RK,RKA,RKB= ',RK,RKA,RKB
        B=SQRT((RKB-RK)*(RKB+RK))
        IF(LOVE.GT.0) GOTO 30
        A=SQRT((RKA-RK)*(RKA+RK))
        AB=A/B
        OMEGA=-RK2+0.5*RKB2
        OMEGA2=OMEGA*OMEGA
        Y1=A*B*RK2
        Y = Y1+OMEGA2
        RPP=(Y1-OMEGA2)/Y
        RSP=-2.*RK*B*OMEGA/Y
        RPZ=-OMEGA*RKB2/Y*A
        P(2)=2.*RK*A
        P(3)=-RK2
        IF(IDTS.EQ.6) THEN
           P(1)=-2.*A*A
           SV(1)=-2.*RK*B
        ELSE
           P(1)=RK2-2.*A*A
           SV(1)=-3.*RK*B
        END IF
        SV(2)=RKB2-2.*RK2
        SV(3)=RK*B
        P2(1)=P(1)
        P2(2)=-P(2)
        P2(3)=P(3)
        DO I=1,3
           P2(I)=P2(I)*RPP
           SV(I)=SV(I)*RSP*AB
        enddo
        RETURN
C NEXT FOR SH WAVES
 30     RPZ=2.0*RK
        P(1)=0.
        P(2)=-RKB2*B/RK
        P(3)=RKB2
        P2(1)=0.
        P2(2)=-P(1)
        P2(3)=P(3)
        RETURN
        END


c	SUBROUTINE FAULTT(THEAT,M1,M2,M3,M4,M5,M6,M00,AM)
c	DIMENSION AM(6)
c	 REAL M1,M2,M3,M4,M5,M6,M00
c	X=0.0174533
c	T=THEAT*X
c	T2=2.0*T
c	 COST=COS(T)
c	 SINT=SIN(T)
c         SINT2=SIN(T2)
c	 COST2=COS(T2)
c	 AM(1)=M1/2.0
c	 AM(2)=-M2*COST-M3*SINT
c	 AM(3)=M4*COST*COST+M5*SINT*SINT+M6*SIN(T2)
c         AM(4)=M2*SINT-M3*COST
c	 AM(5)=-(M4-M5)*SINT2/2.0+M6*COST2
c	 DO 10 I=1,6
c	 AM(I)=AM(I)*M00
c 10      CONTINUE	 
c	RETURN
c	END


c	SUBROUTINE STIME2(TA0,DTA)
c        PARAMETER (MMSOU=50)
c	 PARAMETER ( LINST=2800,LGREEN=3660,LDIM=2048,inpt=2050)
c	COMMON/CO1/ DR(LDIM),DI(LDIM)
c        CHARACTER*12 FRTDK
c        COMMON/PARAM/DT,LNPT,NPT,DF,JF,AW,FRTDK
c	COMPLEX SOURCE
c	COMMON/FIGSS/MSOU,NSOU,SOURCE(inpt,MMSOU)
c        TIME=10.8
cC        TIME=4.4
c        K=INT(TIME/DT+0.1)
cC	WRITE(*,*)' K= ',K
c        IF(NPT.GT.1024) THEN
c	WRITE(*,*)' IN STIME2 NPT>1024, STOP!'
c	END IF
c        DT2=0.0
c        DO 10 IS=1,MSOU
c        DO 1 I=1,2050
c        DI(I)=0.0
c 1      DR(I)=0.0	      	
c	TA=TA0+(IS-1)*DTA
c	DT1=TA/2.0
c	DT3=DT1
c	N1=INT(TA/DT/2.+0.1)
c	NN=2*(N1+1)
c	NSOU=2*N1+1
c	DO 11 I=1,N1
c	DR(K+I)=(I-1.0)/N1
c11	DR(K+NN-I)=(I-1.0)/N1
c	DR(K+N1+1)=1.0
c        NSOU=NSOU+K
cC	WRITE(*,*)'NSOU= ',NSOU
c	SUM=0.0
c	DO 21 I=1,NSOU
c21	SUM=SUM+DR(I)
c        SUM=SUM*DT
c	DO 22 I=1,NSOU
c22	DR(I)=DR(I)/SUM
c	DO 23 I=NSOU+1,NPT
c23	DR(I)=DR(NSOU)
c	CALL FFT(DR,DI,LNPT,-1.)
c        DO 24 J=1,JF
c	SOURCE(J,IS)=CMPLX(DR(J),DI(J))
c 24     CONTINUE	
c10      CONTINUE
c	RETURN
c	END


c	SUBROUTINE LOG2FD(NP,N,L2N)
c	N1=0
c	N=1
c	L2N=0
c 1	CONTINUE
c	IF(NP.GT.N1.AND.NP.LE.N) GOTO 2
c	N1=N1*2
c	N=N*2
c	L2N=L2N+1
c	GOTO 1
c 2	RETURN
c	END


c
c Related to Fast Fourier transform
c
        SUBROUTINE REALTR(XR,XI,N)
        implicit none
        real xr, xi
        integer n
        DIMENSION XR(*),XI(*)
        integer lx, lh, lb, i, i1, i2
        LX=2**N
        LH=LX/2
        LB=LH-1
        LH=LH+1
        DO I=1,LB
           I1=LH+I
           I2=LH-I
           XR(I1)=XR(I2)
           XI(I1)=-XI(I2)
        enddo
        XI(LH)=0.0
        RETURN
        END


        SUBROUTINE FFT(XR,XI,N,SN)
        implicit none
        real xr, xi, sn
        integer n, m
        DIMENSION XR(*),XI(*),M(25)
        integer lx, i, l, nb, lb, lbh, k, ib, fk, flx, ist, j, jh, ii
        integer j1
        real wkr, wki, qr, qi, holdr, holdi
        real pi, v
        pi = 4.0*atan(1.0)
        LX = 2**N
        DO I=1,N
           M(I)=2**(N-I)
        enddo
        DO L=1,N
           NB=2**(L-1)
           LB=LX/NB
           LBH=LB/2
           K=0
           DO IB=1,NB
              FK=K
              FLX=LX
              V=SN*2.0*pi*FK/FLX
              WKR=COS(V)
              WKI=SIN(V)
              IST=LB*(IB-1)
              DO I=1,LBH
                 J=IST+I
                 JH=J+LBH
                 QR=XR(JH)*WKR-XI(JH)*WKI
                 QI=XR(JH)*WKI+XI(JH)*WKR
                 XR(JH)=XR(J)-QR
                 XI(JH)=XI(J)-QI
                 XR(J)=XR(J)+QR
                 XI(J)=XI(J)+QI
              enddo
              DO I=2,N
                 II=I
                 IF(K.LT.M(I)) GOTO 4
                 K=K-M(I)
              enddo
 4            K=K+M(II)
           enddo
        enddo
        K=0
        DO J=1,LX
           IF(K.Lt.J) GOTO 7
           HOLDR=XR(J)
           HOLDI=XI(J)
           J1=K+1
           XR(J)=XR(J1)
           XI(J)=XI(J1)
           XR(J1)=HOLDR
           XI(J1)=HOLDI
 7         DO I=1,N
              II=I
              IF(K.Lt.M(I)) GOTO 6
              K=K-M(I)
           enddo
 6         K=K+M(II)
        enddo
        IF(SN.LT.0.) GOTO 9 !IF(SN.LT.0.) GOTO 9
        DO I=1,LX
           XR(I)=XR(I)/FLX
           XI(I)=XI(I)/FLX
        enddo
 9      CONTINUE
        RETURN
        END


c      SUBROUTINE INTERW (U,NU,DTU,F,NF,DTF)
cC--- DTU=0.1
cC--- DTF=0.2
c      DIMENSION U(*),F(*)
c      K=0
c      DO 2 I=1,NU-5,2
c	K=K+1
c      J1=I
c      J2=I+2
c      J3=I+3
c      J4=I+4
c      J5=I+5
c       F(K)=0.50*U(J3)+(U(J2)+U(4))*0.25
c    2 CONTINUE
c	 NF=K
c      RETURN
c      END
c      SUBROUTINE CORLAT(U,Y,DT,AC,NPS)
c      DIMENSION U(*),Y(*),S(1000)
c      NP=NPS-1
c      DO 1 I=1,NP
c      S(I)=.5*DT*(U(I)*Y(I)+U(I+1)*Y(I+1))
c1     CONTINUE
c      AC=0.
c      DO 2 I=1,NP
c      AC=AC+S(I)
c2     CONTINUE
c      RETURN
c      END


c	SUBROUTINE FAULDM(THEAT,DELT,AL,M1,M2,M3,M4,M5,M6,MM)
c	 REAL M1,M2,M3,M4,M5,MM,M6
c	X=0.0174533
c	T=THEAT*X
c	B=AL*X
c	D=DELT*X
c	T2=2.0*T
c	D2=2.0*D
c	SINB=SIN(B)
c	COSB=COS(B)
c	SINT=SIN(T)
c	COST=COS(T)
c	SIND=SIN(D)
c	COSD=COS(D)
c	SIND2=SIN(D2)
c	COSD2=COS(D2)
c	SINT2=SIN(T2)
c	COST2= COS(T2)
c	 M1=MM*SINB*SIND2
c	 M2=-MM*(COSB*COSD*COST+SINB*COSD2*SINT)
c	 M3=-MM*(COSB*COSD*SINT-SINB*COSD2*COST)
c	 M4=-MM*(SINB*SIND2*SINT*SINT-COSB*SIND*SINT2)
c	 M5=-MM*(SINB*SIND2*COST*COST-COSB*SIND*SINT2)
c	 M6= MM*(SINB*SIND2*SINT2+COSB*SIND*COST2)
c	RETURN
c	END


c	SUBROUTINE FAULD(AZ,THEAT,DELT,AL,MO,AM)
c	DIMENSION AM(6)
c	REAL MO
c	X=0.0174533
c	T=(AZ-THEAT)*X
c	B=AL*X
c	D=DELT*X
c	T2=2.0*T
c	D2=2.0*D
c	SINT2=SIN(T2)
c	COST2= COS(T2)
c	SINB=SIN(B)
c	COSB=COS(B)
c	SIND=SIN(D)
c	COSD=COS(D)
c	SIND2=SIN(D2)
c	COSD2=COS(D2)
c	COST=COS(T)
c	SINT=SIN(T)
c        AM(1)=0.5*MO*SINB*SIND2
c	AM(2)=MO*(COSB*COSD*COST-SINB*COSD2*SINT)
c	AM(3)=MO*(0.5*SINB*SIND2*COST2+COSB*SIND*SINT2)
c	AM(4)=-MO*(SINB*COSD2*COST+COSB*COSD*SINT)
c	AM(5)=MO*(COSB*SIND*COST2-0.5*SINB*SIND2*SINT2)
c	AM(6)=0.
c	RETURN
c	END


      SUBROUTINE INTEDP(X,Y,M,depth)
      implicit none
      real x, y, depth
      integer m
      real XP(16),YP(16),YPP100(16),ypp0(16),ypp200(16)
      real YPP300(16),ypp400(16),ypp500(16),ypp600(16),ypp700(16)
      REAL DIF1,DIF2,DIFY,DR, ratio
      integer k, n, i
      DATA YPP0/8.4,7.95,7.65,7.4,7.15,6.9,6.7,6.35,5.95,5.6,5.3,
     *5.0,4.7,4.4,4.1,3.85/
      DATA YPP100/9.32,8.38,7.43,7.38,7.32,7.24,7.16,6.96,6.39,5.93,
     *5.47,5.02,5.07,5.20,4.71,3.56/
      DATA YPP200/9.21,8.20,8.15,8.17,8.08,7.97,7.59,6.78,6.56,6.34,
     *6.06,5.56,5.40,4.89,4.35,3.25/
      DATA YPP300/10.18,9.49,9.11,8.94,8.80,8.75,8.40,8.05,7.53,6.89,
     *6.47,6.38,6.10,5.65,5.18,3.99/
      DATA YPP400/11.04,10.53,10.51,10.64,10.26,9.64,9.37,9.09,8.51,
     *7.44,7.24,6.98,6.77,6.37,5.55,4.09/
      DATA YPP500/12.65,12.02,11.49,11.11,10.69,10.64,9.99,9.13,8.51,
     *8.02,7.32,7.02,6.83,6.64,5.61,4.13/
      DATA YPP600/13.75,13.50,13.29,12.93,12.50,11.86,11.46,10.44,9.56,
     *9.12,9.18,8.80,8.41,8.22,7.14,5.12/
      DATA YPP700/17.54,16.66,16.08,15.84,15.25,14.56,12.21,6.40,6.10,
     *10.94,10.31,9.89,9.31,8.37,8.11,5.89/
      DATA XP/30.,32.5,35.,37.5,40.,42.5,45.,50.,55.,60.,65.,70.
     * ,75.,80.,85.,90./ 
      IF(M.EQ.1)THEN
         Y=1.0
         RETURN
      endif
      if(x.gt.90)then
         y=2.25
         return
      END IF
      if(depth.lt.33)then
         do k=1,16
            yp(k)=ypp0(k)
         enddo
      endif 

      if(depth.ge.33. and. depth.lt.100)then
         do k=1,16
            yp(k)=ypp0(k)+(ypp100(k)-ypp0(k))*depth/100
         enddo
      elseif(depth.ge.100.and.depth.lt.200)then
         do k=1,16
            yp(k)=ypp100(k)+(ypp200(k)-ypp100(k))*(depth-100.)/100.0
         enddo
      elseif(depth.ge.200.and.depth.lt.300)then 
         do k=1,16
            yp(k)=ypp200(k)+(ypp300(k)-ypp200(k))*(depth-200.)/100.0
         enddo
      elseif(depth.ge.300.and.depth.lt.400)then 
         do k=1,16
            yp(k)=ypp300(k)+(ypp400(k)-ypp300(k))*(depth-300.)/100.0
         enddo
      elseif(depth.ge.400.and.depth.lt.500)then 
         do k=1,16
            yp(k)=ypp400(k)+(ypp500(k)-ypp400(k))*(depth-400.)/100.0
         enddo
      elseif(depth.ge.500.and.depth.lt.600)then 
         do k=1,16
            yp(k)=ypp500(k)+(ypp600(k)-ypp500(k))*(depth-500.)/100.0
         enddo
      elseif(depth.ge.600.and.depth.lt.700)then 
         do k=1,16
            yp(k)=ypp600(k)+(ypp700(k)-ypp600(k))*(depth-600.)/100.0
         enddo
      elseif(depth.ge.700)then 
         do k=1,16
            yp(k)=ypp700(k)
         enddo
      endif

      IF (X.GT.XP(16)) then
         y = yp(16)
         return
      endif
      IF (X .LT. XP(1)) then
         y = yp(1)
         return
      endif
      N=16
      DO I=1,N
         if ((xp(i) - x) .le. -1.e-6) cycle
         if (abs(xp(i) - x) .lt. 1.e-6) then
            y = yp(i)
         else
            K=I-1
            DIF1=XP(I)-XP(K)
            DIF2=XP(I)-X
            RATIO = DIF2/DIF1
            DIFY=ABS(YP(I) - YP(K))
            DR = DIFY*RATIO
            IF (YP(I) .GT. YP(K)) then
               y = yp(i)-dr
            else
               Y= YP(I)+DR
            endif
         endif
         return
      enddo
c      IF (X.GT.XP(16))  GO TO 66
c      IF (X .LT. XP(1))  GO TO 6
c      N=16
c      DO 10  I=1,N
c         IF (XP(I)-X) 10,102,3
c   10       CONTINUE
c    3       K=I-1
c            DIF1=XP(I)-XP(K)
c            DIF2=XP(I)-X
c            RATIO = DIF2/DIF1
c            DIFY=ABS(YP(I) - YP(K))
c            DR = DIFY*RATIO
c            IF (YP(I) .GT. YP(K))  GO TO 4
c            Y=YP(I) + DR
c            GOTO 70
c    4       Y=YP(I)-DR
c            GOTO 70
c  102       Y=YP(I)
c            GOTO 70
c    6 Y=YP(1)
c      GOTO 70
c  66  Y=YP(16)
c  70  CONTINUE
      RETURN     
      END


      SUBROUTINE RTDK(LOVE,PO,DEPTH,NBS)
      implicit none
      real po, depth
      integer love, nbs
      integer inpt, inptd, inpth, nll
      PARAMETER(INPT=2050,INPTD=4100,INPTH=1025,NLL=45)
      integer idts
      COMMON/SOURECTYPE/IDTS
      real dt, df, aw
      integer lnpt, npt, jf, m
      CHARACTER*12 FNAME1
C      IDTS=6 FOR MOMENT SOURCE
C      IDTS=3 FOR DISLOCATION SOURCE	
      COMMON/PARAM/DT,LNPT,NPT,DF,JF,AW,FNAME1
      COMPLEX CFW1,CFV1,CFW2,CFV2,FINST,CFW0,CTIMED(INPTH)
      COMMON/FIGG1/CFW0(INPTH),CFW1(INPTH),FINST(INPTH)
      COMMON/FIGG2/CFW2(INPTH),CFV1(INPTH),CFV2(INPTH)
c	COMMON/FIGSS/MSOU,NSOU,SOURCE(1024,MMSOU)
      real afa, bta, ro, thk, qa, qb
      COMMON/MODEL/AFA(NLL),BTA(NLL),RO(NLL),THK(NLL),QA(NLL),QB(NLL)
      COMPLEX RESPN,RDSL,RESPSH,RDSLSH,PD,PU,SVD,SVU,SHD,SHU
      COMMON/RESSAV/RESPN(2,2),RDSL(2,2),RESPSH,RDSLSH
     *,PU(4),PD(4),SVD(4),SVU(4),SHD(3),SHU(3)

      COMPLEX RDPS,TDPS,RUPS,TUPS,A,B,TERM,RSUR
      COMMON/SAVE/RDPS(2,2,NLL),TDPS(2,2,NLL),RUPS(2,2,NLL),
     * TUPS(2,2,NLL),A(NLL),B(NLL),TERM(2,2),RSUR(2,2)
c
c extra variables
c
      COMPLEX CZERO,CI,UWV(4,4), cw(inpt), z, w, w4, w2, CVK1,CVK2,
     * CWK0,CWK1,CWK2 
      real ws(inpt), tc(inpt), twopi, pi, pid2, fk, fk2, tstart, wi,
     * tpsr
      integer i, iflag, ifr, jf2, nl, ib

c	write(*,*)'inside RTDK'     
      pi = 4.0*atan(1.0) 
      TWOPI=2.0*pi
      PID2=.5*PI
      CZERO=CMPLX(0.,0.)
      CI=CMPLX(0.,1.)
      M=1
      CALL CMODEL(NL,DEPTH,NBS)
      CALL FIND2(NBS,NL,PO,TPSR,LOVE)
C
C     IDTS=3 DISLOCATION SOURCE ; IDTS=6 FOR Mij
C 
C        READ(15,*)DT,TSTART
      TSTART=0.0
C
C      NFR=NPT/2
C      TW=NPT*DT
C      DF=1./TW
      WI=-PI*DF*AW
C      JF=NFR+1
      JF2=2*JF
      DO I=1,JF
         WS(I)=TWOPI*(I-1)*DF
         WS(1)=0.000001
         CW(I)=CMPLX(WS(I),WI)
         Z=CI*CW(I)*TPSR
         Z=CEXP(Z)
         CTIMED(I)=Z
      enddo
C
      DO I=1,JF
c      IF(LOVE.GT.0) THEN
         CFV1(I)=CZERO
         CFV2(I)=CZERO
c      ELSE
         CFW0(I)=CZERO
         CFW1(I)=CZERO
         CFW2(I)=CZERO
c      END IF
      enddo
      DO I=1,NPT
         TC(I)=(I-1)*DT
      enddo

      IB=2
      IFLAG=0
      DO IFR=IB,JF
         I=IFR
         FK=WS(I)*PO
         FK2=FK*FK
         W=CW(I)
         W2=W*W
         W4=W2*W2
         CALL RESPON(FK,W,NL,NBS,LOVE,PO,IFLAG)
c         CALL SOUR(FK,FK2,NBS) 
         CALL RESUWV(UWV)
c	if (ifr.eq.3) then
c	write(*,*)'common/ressav'
c	write(*,*)((RESPN(j,k),j=1,2),k=1,2)
c	write(*,*)((RDSL(j,k),j=1,2),k=1,2)
c	write(*,*)RESPSH,RDSLSH
c	write(*,*)(PU(j),j=1,4)
c	write(*,*)(PD(j),j=1,4)
c	write(*,*)(SVD(j),j=1,4)
c	write(*,*)(SVU(j),j=1,4)
c	write(*,*)(SHD(j),j=1,3)
c	write(*,*)(SHU(j),j=1,3)
c	write(*,*)'uwv'
c	write(*,*)((uwv(j,k),j=1,4),k=1,4)
c	endif
         IFLAG=IFLAG+1
C
c      IF(LOVE.GT.0) THEN
c      CVK1=UWV(3,2)*(0.0,1.0)
c      CVK2=-UWV(3,3)
         CVK1=-UWV(3,2)*(0.0,1.0)
         CVK2=UWV(3,3)
         CFV1(I)=CVK1*CTIMED(I)
         CFV2(I)=CVK2*CTIMED(I)
c      ELSE
         CWK0=-UWV(1,1)
         CWK1=-UWV(1,2)*(0.0,1.0)
c      CWK2=UWV(1,3)
         CWK2=UWV(1,3)*sqrt(2.)/2. !*0.707
         CFW0(I)=CWK0*CTIMED(I)
         CFW1(I)=CWK1*CTIMED(I)
         CFW2(I)=CWK2*CTIMED(I)
c      END IF
      enddo
      CLOSE(15)
      RETURN
      END


      SUBROUTINE RETRSH(XK,XK2,N)
      implicit none
      real xk, xk2
      integer n
      integer inpt, inptd, inpth, nll
      PARAMETER(INPT=2050,INPTD=4100,INPTH=1025,NLL=45)
      real afa, bta, ro, thk, qa, qb
      COMMON/MODEL/AFA(NLL),BTA(NLL),RO(NLL),THK(NLL),QA(NLL),QB(NLL)
      complex cka2, ckb2
      real cmu, roh
      COMMON/MOD/CKA2(NLL),CKB2(NLL),CMU(NLL),ROH(NLL)
      complex rd, td, ru, tu, bs
      COMMON/SAVESH/RD(NLL),TD(NLL),RU(NLL),TU(NLL),BS(NLL)
      complex rdps, tdps, rups, tups, a, b, term, rsur
      COMMON/SAVE/RDPS(2,2,NLL),TDPS(2,2,NLL),RUPS(2,2,NLL),
     * TUPS(2,2,NLL),A(NLL),B(NLL),TERM(2,2),RSUR(2,2)
c
c otras variables
c
      complex z11, z0, b1, b2, z, z1, z2, cr2
      integer i1, i2, j
      Z11=CMPLX(1.,0.)
      Z0=CMPLX(0.,0.)
      DO J=1,N
         I1=J-1
         I2=J
         Z=CKB2(J)
         B(J)=CR2(XK2,Z)
         IF(J.EQ.1) cycle!GOTO 100
         IF(ABS(BTA(I2)-BTA(I1)).gt.0.001) then
            B1=B(I1)
            B2=B(I2)
            Z1=CMU(I1)*B1
            Z2=CMU(I2)*B2
            Z=Z1+Z2
            Z1=Z1/Z
            Z2=Z2/Z
            RD(I2)=Z1-Z2
            RU(I2)=-RD(I2)
            TD(I2)=2.*Z1
            TU(I2)=2.*Z2
         else
            RD(I2)=Z0
            RU(I2)=Z0
            TD(I2)=Z11
            TU(I2)=Z11
         endif
      enddo
      RETURN
      END


      SUBROUTINE LAYRSH(XK,RUFS,RESPN,N,NBS)
      implicit none
      real xk
      complex rufs, respn
      integer n, nbs
      integer inpt, inptd, inpth, nll
      PARAMETER(INPT=2050,INPTD=4100,INPTH=1025,NLL=45)
      complex rd, td, ru, tu, bs
      COMMON/SAVESH/RD(NLL),TD(NLL),RU(NLL),TU(NLL),BS(NLL)
      complex rdps, tdps, rups, tups, a, b, term, rsur
      COMMON/SAVE/RDPS(2,2,NLL),TDPS(2,2,NLL),RUPS(2,2,NLL),
     * TUPS(2,2,NLL),A(NLL),B(NLL),TERM(2,2),RSUR(2,2)
      real afa, bta, ro, thk, qa, qb
      COMMON/MODEL/AFA(NLL),BTA(NLL),RO(NLL),THK(NLL),QA(NLL),QB(NLL)
      complex RD1,TD1,RU1,TU1,RD2,TD2,RU2,TU2
      COMMON/ATSH/RD1,TD1,RU1,TU1,RD2,TD2,RU2,TU2
c
c otras variables
c
      complex TDRS,RURS, RDSL,TDSL,RDRS,TURS, RDPH(NLL),TDPH(NLL),
     * RUPH(NLL),TUPH(NLL),Z1,Z0, RD3,TD3,TU3,RU3,AA, B1,EXPB,EXP2B,Z
      integer j, i1, i2, iz, k2, nbr, nn, nz
      real theck

      Z1=CMPLX(1.0,0.)
      Z0=CMPLX(0.0,0.)
      DO J=2,N
         I1=J-1
         I2=J
         B1=B(I1)
         THECK=THK(I2)-THK(I1)
         Z=-THECK*B1
         EXPB=CEXP(Z)
         EXP2B=EXPB*EXPB
         RDPH(J)=RD(J)*EXP2B
         TDPH(J)=TD(J)*EXPB
         TUPH(J)=TU(J)*EXPB
         RUPH(J)=RU(J)
C     IF(IOLYSH.GT.0)WRITE(*,666)EXPB,EXP2B
c666   FORMAT(1X,'EXPB,EXP2B=',4E13.4)
C     IF(IOLYSH.GT.0)WRITE(*,661)RDPH(J),RD(J),TDPH(J),TD(J)
c661   FORMAT(1X,'RDPH(J),RD(J),TDPH(J),TD(J)=',/8E13.4)
      enddo
C
C     FOR RDSL
C
      IF(NBS.EQ.N) THEN
         tdsl = z0
         RDSL=Z0
         GOTO 401
      END IF
      NN=N
      K2=N-2
      RD3=RDPH(NN)
      TD3=TDPH(NN)
C      IF((N-NBS).EQ.1)GOTO 401
      DO IZ=NBS,K2
         NZ=K2+NBS+1-IZ
         RD1=RDPH(NZ)
         TD1=TDPH(NZ)
         RU1=RUPH(NZ)
         TU1=TUPH(NZ)
         RD2=RD3
         TD2=TD3
         AA=Z1-RU1*RD2
         RD3=RD1+TU1*RD2*TD1/AA
         TD3=TD2*TD1/AA
      enddo
      RDSL=RD3
      TDSL=TD3
401   CONTINUE      
C
C     FOR RDRS,TDRS,RURS,TURS
C
      NN=NBS
      K2=NBS-2
      RD3=RDPH(NN)
      TD3=TDPH(NN)
      RU3=RUPH(NN)
      TU3=TUPH(NN)
C
      NBR=1
      DO IZ=NBR,K2
         NZ=K2+NBR+1-IZ
         RD1=RDPH(NZ)
         TD1=TDPH(NZ)
         RU1=RUPH(NZ)
         TU1=TUPH(NZ)
         CALL MATSH(RD3,TD3,RU3,TU3)
      enddo
C
      RDRS=RD3
      TURS=TU3
      RURS=RU3
      TDRS=TD3
C
C     RUFS
C
      AA=Z1-RDRS
      RUFS=RURS+TDRS*TURS/AA
C      
C     RESPN 
C 
      AA=Z1-RUFS*RDSL
      RESPN=TDSL/AA

C     IF(IOLYSH.LT.1) GOTO 300
C     WRITE(*,610)RDSL
C     WRITE(*,630)RDRS
C     WRITE(*,640)TURS
C     WRITE(*,650)RUFS
C     IF(IOLYSH.LT.1) GOTO 950
c660   FORMAT(1X,'RESPONS=',/4E13.4)
C      WRITE(*,660)RESPN
c950   CONTINUE
      RETURN
      END


      SUBROUTINE MATSH(RD3,TD3,RU3,TU3)
      implicit none
      complex rd3, td3, ru3, tu3
      COMPLEX RD1,TD1,RU1,TU1,RD2,TD2,RU2,TU2
      COMMON/ATSH/RD1,TD1,RU1,TU1,RD2,TD2,RU2,TU2
c
c variables extra
c
      COMPLEX AA,BB,Z1
      Z1=CMPLX(1.0,0.)
      RD2=RD3
      TD2=TD3
      RU2=RU3
      TU2=TU3
      AA=Z1-RU1*RD2
      BB=TD1/AA
      TD3=TD2*BB
      AA=TU1*RD2
      AA=AA*BB
      RD3=RD1+AA
      AA=Z1-RD2*RU1
      BB=TU2/AA
      TU3=TU1*BB
      AA=TD2*RU1
      AA=BB*AA
      RU3=RU3+AA
      RETURN
      END


      SUBROUTINE SOUR(XK,XK2,NBS)
      implicit none
      real xk, xk2
      integer nbs
      integer inpt, inptd, inpth, nll
      PARAMETER(INPT=2050,INPTD=4100,INPTH=1025,NLL=45)
      COMPLEX RESPN,RDSL,RESPSH,RDSLSH,PD,PU,SVD,SVU,SHD,SHU
      COMMON/RESSAV/RESPN(2,2),RDSL(2,2),RESPSH,RDSLSH
     $,PD(4),PU(4),SVD(4),SVU(4),SHD(3),SHU(3)
      real afa, bta, ro, thk, qa, qb
      COMMON/MODEL/AFA(NLL),BTA(NLL),RO(NLL),THK(NLL),QA(NLL),QB(NLL)
      complex ckasq, ckbsq
      real cmu,roh
      COMMON/MOD/CKASQ(NLL),CKBSQ(NLL),CMU(NLL),ROH(NLL)
      COMPLEX RD,TD,RU,TU,ACOM,BCOM,TERM,RSUR
      COMMON/SAVE/RD(2,2,NLL),TD(2,2,NLL),RU(2,2,NLL),TU(2,2,NLL),ACOM(
     * NLL),BCOM(NLL),TERM(2,2),RSUR(2,2)
      integer idts
        COMMON/SOURECTYPE/IDTS
c
c variables extras
c
      COMPLEX CKA2,CKB2,A,B,Z,Z1,z0
C      IDTS=6 FOR MOMENT SOURCE
C      IDTS=3 FOR DISLOCATION SOURCE	
      Z0=CMPLX(0.,0.)
      Z1=CMPLX(1.0,0.)
      CKA2=CKASQ(NBS)
      CKB2=CKBSQ(NBS)
      A=ACOM(NBS)
      B=BCOM(NBS)
C      IF(M.EQ.0) GOTO 100
C      Z=-XK2/A
      Z=-XK2
      PD(3)=Z
      PU(3)=Z
C       Z=Z1*2.*XK
      Z=Z1*2.*XK*A
      PD(2)=Z
      PU(2)=-Z
C
C     FOR DISLOCATION SOURCE
C
      IF(IDTS.EQ.3) THEN
C      Z=(2.*CKA2-3.*XK2)/A
         Z=(2.*CKA2-3.*XK2)
         PD(1)=Z
         PU(1)=Z
C      Z=-3.*XK*Z1
         Z=-3.*XK*Z1*A
         SVD(1)=Z
         SVU(1)=-Z
      END IF
C
C     FOR Mij SOURCE
C
      IF(IDTS.EQ.6) THEN
C      Z=-2.0*A
         Z=-2.0*A*A
         PD(1)=Z
         PU(1)=Z
C      Z=-2.*XK*Z1
         Z=-2.*XK*Z1*A
         SVD(1)=Z
         SVU(1)=-Z
      END IF
C
C      Z=(2.*XK2-CKB2)/B
      Z=(2.*XK2-CKB2)*A/B
      SVD(2)=Z
      SVU(2)=Z
C      Z=-Z1*XK
      Z=-Z1*XK*A
      SVD(3)=Z
      SVU(3)=-Z
C      Z=CKB2/B
      Z=CKB2
      SHD(3)=Z
      SHU(3)=Z
C       Z=-CKB2/XK
      Z=-B*CKB2/XK
      SHD(2)=Z
      SHU(2)=-Z
C     IF(IOSOUR.EQ.0) GOTO 300
C     DO 400 I=1,3
C400   WRITE(*,410)PD(I),PU(I),SVD(I),SVU(I),SHD(I),
C    * SHU(I)
      RETURN
      END


      SUBROUTINE RESUWV(UWV)
      implicit none
      complex uwv(4,4)
      integer inpt, inptd, inpth, nll
      PARAMETER(INPT=2050,INPTD=4100,INPTH=1025,NLL=45)
      COMPLEX RESPN,RUFS,RESPSH,RUFSSH,PD,PU,SVD,SVU,SHD,SHU
      COMMON/RESSAV/RESPN(2,2),RUFS(2,2),RESPSH,RUFSSH
     $,PD(4),PU(4),SVD(4),SVU(4),SHD(3),SHU(3)
c
c variables extras
c
      COMPLEX U(2),D(2),B(2),Z,Z1,Z0
      integer i, j, k, m2, m3

      Z0=CMPLX(0.,0.)
      Z1=CMPLX(1.,0.)
      M3=3
      M2=2
C      IF(M.EQ.0)M3=1
C      IF(M.EQ.4)M3=4
      DO I=1,M3
         D(1)=PD(I)
         D(2)=SVD(I)
         U(1)=PU(I)
         U(2)=SVU(I)
         DO J=1,M2
            Z=Z0
            DO K=1,M2
               Z=Z+RUFS(J,K)*U(K)
            enddo
            B(J)=Z+D(J)
         enddo
         DO J=1,M2
            Z=Z0
            DO K=1,M2
               Z=Z+RESPN(J,K)*B(K)
            enddo
            UWV(J,I)=Z
         enddo
         IF(I.EQ.1) cycle!GOTO 100
C********************************
         Z=RUFSSH*SHU(I)+SHD(I)
         UWV(3,I)=RESPSH*Z
      enddo
C     IF(IOTUWV.EQ.0) GOTO 200
C     WRITE(*,4)
C     WRITE(*,1)(UWV(1,J),J=1,3)
C     WRITE(*,2)(UWV(2,J),J=1,3)
C     WRITE(*,3)(UWV(3,J),J=1,3)
c200   CONTINUE
      RETURN
      END


      SUBROUTINE REFTRA(XK,XK2,N)
      implicit none
      real xk, xk2
      integer n
      integer inpt, inptd, inpth, nll
      PARAMETER(INPT=2050,INPTD=4100,INPTH=1025,NLL=45)
      real afa, bta, ro, thk, qa, qb
      COMMON/MODEL/AFA(NLL),BTA(NLL),RO(NLL),THK(NLL),QA(NLL),QB(NLL)
      complex cka2, ckb2
      real cmu, roh
      COMMON/MOD/CKA2(NLL),CKB2(NLL),CMU(NLL),ROH(NLL)
      complex rd, td, ru, tu, a, b, term, rsur
      COMMON/SAVE/RD(2,2,NLL),TD(2,2,NLL),RU(2,2,NLL),TU(2,2,NLL)
     $,A(NLL),B(NLL),TERM(2,2),RSUR(2,2)
c
c variables extras
c
      complex CR2,A1,A2,B1,B2,AB1,AB2,AB12,Z1,ROH1,ROH2,CD,A1B2,B1A2,
     * Z,CQ1,CQ2,CQ3,CQ4,Q1,Q2,Q3,Q4,Q5,Q6,D,OMEGA,ROH12
      real x
      integer i1, i2, j, isafe, k1, k2
C
      Z=CKA2(1)
      A(1)=CR2(XK2,Z)
      Z=CKB2(1)
      B(1)=CR2(XK2,Z)
      DO 100 J=1,N
         IF(J.EQ.1) GOTO 199
         Z=CKA2(J)
         A(J)=CR2(XK2,Z)
         Z=CKB2(J)
         B(J)=CR2(XK2,Z)
         I1=J-1
         I2=J
         A1=A(I1)
         B1=B(I1)
         A2=A(I2)
         B2=B(I2)
         IF(ABS(AFA(I2)-AFA(I1)).LT.0.001) GOTO 200
         AB1=A1*B1
         roh2=ckb2(i2)*cmu(i2)
         AB2=A2*B2
         AB12=AB1*AB2
         CQ1=CMU(I2)-CMU(I1)
         ROH1=CKB2(I1)*CMU(I1)
         ROH12=0.5*ROH1
         CQ4=-CQ1*XK2-ROH12
         Z=0.5*ROH2
         CQ2=CQ4+Z
         CQ3=CQ1*XK2-Z
         Q1=XK2*AB12*CQ1*CQ1
         Q2=XK2*CQ2*CQ2
         Q3=AB1*CQ3*CQ3
         Q4=AB2*CQ4*CQ4
         ROH12=ROH12*Z
         A1B2=A1*B2
         Q5=ROH12*A1B2
         B1A2=B1*A2
         Q6=ROH12*B1A2
         D=Q1+Q2-Q3-Q4-Q5-Q6
C
         X=CABS(D)
         IF(X.lt.1.0E-15) then
            ISAFE=1
            WRITE(*,605)
605         FORMAT(1X,'D=0.0! STOP IN REFTRA')
            STOP
         endif
c         IF(X.GT.1.0E-15) GOTO 606
c         ISAFE=1
c         WRITE(*,605)
c605      FORMAT(1X,'D=0.0! STOP IN REFTRA')
c         STOP
c606      CONTINUE
C
         Z=CQ3*B1+CQ4*B2
         Z=CD(Z,D)
         TD(1,1,I2)=ROH1*A1*Z
         TU(1,1,I2)=ROH2*A2*Z
         Z=XK*(CQ1*A1B2+CQ2)
         Z=CD(Z,D)
         TD(1,2,I2)=ROH1*B1*Z
         TU(2,1,I2)=ROH2*A2*Z
         Z=XK*(CQ1*B1A2+CQ2)
         Z=CD(Z,D)
         TD(2,1,I2)=ROH1*A1*Z
         TU(1,2,I2)=ROH2*B2*Z
         Z=A1*CQ3+CQ4*A2
         Z=CD(Z,D)
         TD(2,2,I2)=ROH1*B1*Z
         TU(2,2,I2)=ROH2*B2*Z
         Z=Q1-Q2-Q3+Q4
         Z1=Z-Q5+Q6
         RD(1,1,I2)=CD(Z1,D)
         Z1=Z+Q5-Q6
         RD(2,2,I2)=CD(Z1,D)
         Z=Q1-Q2+Q3-Q4
         Z1=Z+Q5-Q6
         RU(1,1,I2)=CD(Z1,D)
         Z1=Z-Q5+Q6
         RU(2,2,I2)=CD(Z1,D)
         Z=-2.*XK*(CQ3*CQ2-AB2*CQ1*CQ4)
         Z=CD(Z,D)
         RD(1,2,I2)=B1*Z
         RD(2,1,I2)=A1*Z
         Z=-2.*XK*(CQ4*CQ2-AB1*CQ1*CQ3)
         Z=CD(Z,D)
         RU(1,2,I2)=B2*Z
         RU(2,1,I2)=A2*Z
         GOTO 100
C
199      CONTINUE
         A1=A(1)
         B1=B(1)
         AB1=A1*B1
         OMEGA=XK2-0.5*CKB2(1)
         Z=XK2*AB1-OMEGA*OMEGA
C
         X=CABS(Z)
         IF(X.lt.1.0E-20) then
            ISAFE=1
            WRITE(*,907)
907         FORMAT(1X,'Z=0.0! STOP IN REFTRA FOR FIRST LAYER')
            STOP
         endif
c         IF(X.GT.1.0E-20) GOTO 906
c         ISAFE=1
c         WRITE(*,907)
c907      FORMAT(1X,'Z=0.0! STOP IN REFTRA FOR FIRST LAYER')
c         STOP
c906      CONTINUE
         Z1=CKB2(1)*XK*AB1
         Z1=CD(Z1,Z)
         TERM(1,1)=Z1
         TERM(2,2)=Z1
         Z1=CKB2(1)*OMEGA
         Z1=CD(Z1,Z)
         TERM(1,2)=B1*Z1
         TERM(2,1)=A1*Z1
         Z1=XK2*AB1+OMEGA*OMEGA
         Z1=CD(Z1,Z)
         RSUR(1,1)=Z1
         RSUR(2,2)=Z1
         Z1=2.*XK*OMEGA
         Z1=CD(Z1,Z)
         RSUR(1,2)=B1*Z1
         RSUR(2,1)=A1*Z1
         GOTO 100
200      Z=CMPLX(0.,0.)
         Z1=CMPLX(1.,0.)
         DO K1=1,2
            DO K2=1,2
               RD(K1,K2,I2)=Z
               RU(K1,K2,I2)=Z
               TD(K1,K2,I2)=Z
               TU(K1,K2,I2)=Z
            enddo
         enddo
         TD(1,1,I2)=Z1
         TD(2,2,I2)=Z1
         TU(1,1,I2)=Z1
         TU(2,2,I2)=Z1
100   CONTINUE
      RETURN
      END


      SUBROUTINE LAYER(RUFS,RESPN,N,NBS)
      implicit none
      complex rufs, respn
      integer n, nbs
      integer inpt, inptd, inpth, nll
      PARAMETER(INPT=2050,INPTD=4100,INPTH=1025,NLL=45)
      complex RD1,TD1,RU1,TU1,RD2,TD2,RU2,TU2
      COMMON/ATPSV/RD1(2,2),TD1(2,2),RU1(2,2),TU1(2,2),
     *  RD2(2,2),TD2(2,2),RU2(2,2),TU2(2,2)
      real afa, bta, ro, thk, qa, qb
      COMMON/MODEL/AFA(NLL),BTA(NLL),RO(NLL),THK(NLL),QA(NLL),QB(NLL)
      complex rdph, tdph, ruph, tuph
      COMMON/TIMERT/RDPH(2,2,NLL),TDPH(2,2,NLL),RUPH(2,2,NLL)
     *,TUPH(2,2,NLL)
      integer iflag
      COMMON/TESTIN/IFLAG
      complex rd, td, ru, tu, a, b, term, rsur
      COMMON/SAVE/RD(2,2,NLL),TD(2,2,NLL),RU(2,2,NLL),TU(2,2,NLL)
     $,A(NLL),B(NLL),TERM(2,2),RSUR(2,2)
      complex cgmau, cgmad
      COMMON/SOURCE/CGMAU(2),CGMAD(2)

      COMPLEX RD3,TD3,TU3,RU3,AA,BB,RDSL
      DIMENSION RESPN(2,2)
      DIMENSION RDSL(2,2),RUFS(2,2)
      DIMENSION AA(2,2),BB(2,2)
      DIMENSION RD3(2,2),TD3(2,2),RU3(2,2),TU3(2,2)
c
c variables extras
c
      complex CC(2,2),Z0,Z1,DD(2,2),TDSL(2,2),RDRS(2,2),TDRS(2,2),
     * RURS(2,2),TURS(2,2),A1,B1,EXPA,EXPB,EXP2A,EXP2B,EXPAB,Z
      integer i, i1, i2, j, iz, k1, k2, nbr, nsave, nn, nz
      real theck

c	write(*,*)'----'
c	write(*,*)N, nbs
      NSAVE=N
      Z0=CMPLX(0.,0.)
      Z1=CMPLX(1.,0.)
      DO J=2,N
         I1=J-1
         I2=J
         A1=A(I1)
         B1=B(I1)
         THECK=THK(I2)-THK(I1)
         Z=-THECK*A1
         EXPA=CEXP(Z)
         Z=-THECK*B1
         EXPB=CEXP(Z)
         EXP2A=EXPA**2
         EXP2B=EXPB**2
         EXPAB=EXPA*EXPB
         TDPH(1,1,J)=TD(1,1,J)*EXPA
         TDPH(1,2,J)=TD(1,2,J)*EXPB
         TDPH(2,1,J)=TD(2,1,J)*EXPA
         TDPH(2,2,J)=TD(2,2,J)*EXPB
C      WRITE(*,662)TDPH(1,1,J),TD(1,1,J)

C     WRITE(*,661)TDPH(1,2,J),TD(1,2,J)
C     WRITE(*,663)EXPA,EXPB
c4000    CONTINUE
         RDPH(1,1,J)=RD(1,1,J)*EXP2A
         RDPH(1,2,J)=RD(1,2,J)*EXPAB
         RDPH(2,1,J)=RD(2,1,J)*EXPAB
         RDPH(2,2,J)=RD(2,2,J)*EXP2B
         DO K1=1,2
            DO K2=1,2
               RUPH(K1,K2,J)=RU(K1,K2,J)
            end do
         end do
         TUPH(1,1,J)=TU(1,1,J)*EXPA
         TUPH(1,2,J)=TU(1,2,J)*EXPA
         TUPH(2,1,J)=TU(2,1,J)*EXPB
         TUPH(2,2,J)=TU(2,2,J)*EXPB
C     IF(CABS(EXPB).GT.1.0E-10)NSAVE=NSAVE+1
      enddo
C
C     FOR RDSL
C      
      NN=N
      K2=N-2
      IF(NBS.EQ.N) THEN
         DO i=1,2
            DO j=1,2
               tdsl(i,j)=z0 ! correcting a bug previously found!!!!!!!!!!!!!!!
c
c here we correct a bug. previously, we didnt define this thing when nbs==n, thus it could take
c any value. thus our green functions had a random component
c
               rdsl(i,j)=z0
            enddo
         enddo
         GOTO 401
      END IF
      DO I=1,2
         DO J=1,2
            RD3(I,J)=RDPH(I,J,NN)
            TD3(I,J)=TDPH(I,J,NN)
         enddo
      enddo
      DO IZ=NBS,K2
         NZ=K2+NBS+1-IZ
         DO I=1,2
            DO J=1,2
               RD1(I,J)=RDPH(I,J,NZ)
               TD1(I,J)=TDPH(I,J,NZ)
               RU1(I,J)=RUPH(I,J,NZ)
               TU1(I,J)=TUPH(I,J,NZ)
               RD2(I,J)=RD3(I,J)
               TD2(I,J)=TD3(I,J)
            enddo
         enddo
         CALL INVMAT(RU1,RD2,AA)
         CALL MATPRO(AA,TD1,BB)
         CALL MATPRO(TU1,RD2,CC)
         CALL MATPRO(CC,BB,AA)
         DO I=1,2
            DO J=1,2
               RD3(I,J)=RD1(I,J)+AA(I,J)
            enddo
         enddo
         CALL MATPRO(TD2,BB,TD3)
      enddo
      DO I=1,2
         DO J=1,2
            TDSL(I,J)=TD3(I,J)
            RDSL(I,J)=RD3(I,J)
         enddo
      enddo
401   CONTINUE
C
C     FOR RDRS,TURS,RURS,TURS
C      
      NN=NBS
      K2=NBS-2
      DO I=1,2
         DO J=1,2
            RD3(I,J)=RDPH(I,J,NN)
            TD3(I,J)=TDPH(I,J,NN)
            RU3(I,J)=RUPH(I,J,NN)
            TU3(I,J)=TUPH(I,J,NN)
         enddo
      enddo
      NBR=1
      DO IZ=NBR,K2
         NZ=K2+NBR+1-IZ
         DO I=1,2
            DO J=1,2
               RD1(I,J)=RDPH(I,J,NZ)
               TD1(I,J)=TDPH(I,J,NZ)
               RU1(I,J)=RUPH(I,J,NZ)
               TU1(I,J)=TUPH(I,J,NZ)
            enddo
         enddo
         CALL MATPSV(RD3,TD3,RU3,TU3)
      enddo
      DO I=1,2
         DO J=1,2
            RDRS(I,J)=RD3(I,J)
            TDRS(I,J)=TD3(I,J)
            RURS(I,J)=RU3(I,J)
            TURS(I,J)=TU3(I,J)
         enddo
      enddo
C      WRITE(*,*)'RDRS= ',RD3

C
C     FOR RUFS
C      
C
      CALL INVMAT(RDRS,RSUR,AA)
      CALL MATPRO(AA,TURS,BB)
      CALL MATPRO(TDRS,RSUR,CC)
      CALL MATPRO(CC,BB,DD)
      DO I=1,2
         DO J=1,2
            RUFS(I,J)=RURS(I,J)+DD(I,J)
         enddo
      enddo

C      WRITE(*,*)'RUFS= ',RUFS
C
C     FOR RESPN
C
      CALL INVMAT(RUFS,RDSL,AA)
      CALL MATPRO(TDSL,AA,RESPN)

C      WRITE(*,*)'RESPN= ',RESPN

C      IF(IOLYPS.LT.1) GOTO 950
C      WRITE(*,660)RESPN(1,1),RESPN(1,2)
C      WRITE(*,620)RESPN(2,1),RESPN(2,2)
c950   CONTINUE
      RETURN
      END


      SUBROUTINE MATPSV(RD3,TD3,RU3,TU3)
      implicit none
      complex rd3(2,2), td3(2,2), ru3(2,2), tu3(2,2)
      COMPLEX RD1,TD1,RU1,TU1,RD2,TD2,RU2,TU2
      COMMON/ATPSV/RD1(2,2),TD1(2,2),RU1(2,2),TU1(2,2),
     *  RD2(2,2),TD2(2,2),RU2(2,2),TU2(2,2)
c
c variables extra
c
      COMPLEX CC(2,2),AA(2,2),BB(2,2)
      integer i, j
      DO I=1,2
         DO J=1,2
            RD2(I,J)=RD3(I,J)
            TD2(I,J)=TD3(I,J)
            RU2(I,J)=RU3(I,J)
            TU2(I,J)=TU3(I,J)
         enddo
      enddo
      CALL INVMAT(RU1,RD2,AA)
      CALL MATPRO(AA,TD1,BB)
      CALL MATPRO(TD2,BB,TD3)
      CALL MATPRO(TU1,RD2,CC)
      CALL MATPRO(CC,BB,AA)
      DO I=1,2
         DO J=1,2
            RD3(I,J)=RD1(I,J)+AA(I,J)
         enddo
      enddo
      CALL INVMAT(RD2,RU1,AA)
      CALL MATPRO(AA,TU2,BB)
      CALL MATPRO(TU1,BB,TU3)
      CALL MATPRO(TD2,RU1,CC)
      CALL MATPRO(CC,BB,AA)
      DO I=1,2
         DO J=1,2
            RU3(I,J)=RU2(I,J)+AA(I,J)
         enddo
      enddo
      RETURN
      END


      SUBROUTINE RESPON(XK,W,N,NBS,LOVE,PO,IFLAG)
      implicit none
      real xk, po
      complex w
      integer n, nbs, love, iflag
      integer inpt, inptd, inpth, nll
      PARAMETER(INPT=2050,INPTD=4100,INPTH=1025,NLL=45)
      complex shrd, shtd, shru, shtu, bsh
      COMMON/SAVESH/SHRD(NLL),SHTD(NLL),SHRU(NLL),SHTU(NLL),BSH(NLL)
      complex cgmau, cgmad
      COMMON/SOURCE/CGMAU(2),CGMAD(2)
      complex rdph, tdph, ruph, tuph
      COMMON/TIMERT/RDPH(2,2,NLL),TDPH(2,2,NLL),RUPH(2,2,NLL)
     *,TUPH(2,2,NLL)
      complex cka2, ckb2
      real cmu, roh
      COMMON/MOD/CKA2(NLL),CKB2(NLL),CMU(NLL),ROH(NLL)
      complex rd, td, ru, tu, a, b, term, rsur
      COMMON/SAVE/RD(2,2,NLL),TD(2,2,NLL),RU(2,2,NLL),TU(2,2,NLL),
     $ A(NLL),B(NLL),TERM(2,2),RSUR(2,2)
      real afa, bta, ro, thk, qa, qb
      COMMON/MODEL/AFA(NLL),BTA(NLL),RO(NLL),THK(NLL),QA(NLL),QB(NLL)
      complex respn, rufs, respsh, shrufs, pd, pu, svd, svu, shd, shu
      COMMON/RESSAV/RESPN(2,2),RUFS(2,2),RESPSH,SHRUFS
     *,PD(4), PU(4),SVD(4),SVU(4),SHD(3),SHU(3)
c
c variables extras
c
      complex CW2,CR2,z
      integer i
      real pa, po2, xk2
      CW2=W*W
      XK2=XK*XK
      PO2=PO*PO
      z=cmplx(0.0,1.0)
      PA=4.0*atan(1.0)!3.1415926
      IF(IFLAG.EQ.0) THEN
         DO I=1,N
            CKA2(I)=1.0/AFA(I)/AFA(I)
            CKB2(I)=1.0/BTA(I)/BTA(I)
            CMU(I)=BTA(I)*BTA(I)*RO(I)
         enddo
         IF(LOVE.EQ.0) THEN
            CALL REFTRA(PO,PO2,N)
         ELSE
            CALL RETRSH(PO,PO2,N)
         END IF
         CALL SOUR(PO,PO2,NBS)
      END IF
      DO I=1,N
         CKA2(I)=CW2/AFA(I)/AFA(I)
         CKB2(I)=CW2/BTA(I)/BTA(I)
         CMU(I)=BTA(I)*BTA(I)*RO(I)
         Z=CKA2(I)
         A(I)=CR2(XK2,Z)
         Z=CKB2(I)
         B(I)=CR2(XK2,Z)
      enddo

      IF(LOVE.EQ.0) THEN
C      CALL REFTRA(XK,XK2,N)
         CALL LAYER(RUFS,RESPN,N,NBS)
      ELSE
C      CALL RETRSH(XK,XK2,N)
         CALL LAYRSH(XK,SHRUFS,RESPSH,N,NBS)
      END IF

      RETURN
      END


      SUBROUTINE INVMAT(A,B,C)
      implicit none
      complex a, b, c
      DIMENSION A(2,2),B(2,2),C(2,2)
c
c extra variables
c
      COMPLEX CMAD,ZZ1,ZZ2,Z1,Z0,Z,CMM
C
C     C(I,J)=(I-A*B)**(-1)
C
      Z0=CMPLX(0.,0.)
      Z1=CMPLX(1.,0.)
      C(1,1)=CMAD(A(1,1),B(1,1),A(1,2),B(2,1))
      C(1,2)=CMAD(A(1,1),B(1,2),A(1,2),B(2,2))
      C(2,1)=CMAD(A(2,1),B(1,1),A(2,2),B(2,1))
      C(2,2)=CMAD(A(2,1),B(1,2),A(2,2),B(2,2))
      Z=CMM(C(1,1),C(2,2),C(1,2),C(2,1))
      Z=Z-C(1,1)-C(2,2)+Z1

C      WRITE(*,50)((A(I,J),J=1,2),I=1,2)
C     WRITE(*,50)((B(I,J),J=1,2),I=1,2)
C     WRITE(*,50)((C(I,J),J=1,2),I=1,2)
C     WRITE(*,50)Z

      C(1,2)=C(1,2)/Z
      C(2,1)=C(2,1)/Z
      ZZ1=(Z1-C(2,2))/Z
      ZZ2=(Z1-C(1,1))/Z
      C(1,1)=ZZ1
      C(2,2)=ZZ2

C     WRITE(*,50)((C(I,J),J=1,2),I=1,2)
      RETURN
      END


       SUBROUTINE MATPRO(A,B,C)
       implicit none
       complex a,b,c
       DIMENSION A(2,2),B(2,2),C(2,2)
c
c variables extra
c
       COMPLEX CM,Z,Z1
       integer i, j, k
       Do I=1,2
          DO J=1,2
             Z=CMPLX(0.,0.)
             DO K=1,2
                Z1=CM(A(I,K),B(K,J))
                Z=Z+Z1
             enddo
             C(I,J)=Z
          enddo
       enddo
       RETURN
       END


       COMPLEX FUNCTION CR2(P2,C2)
       implicit none
       real p2
       integer iflag1
       COMMON/TESTIN/IFLAG1
       DOUBLE PRECISION U,X,R,W1,W2,R1,R2
       COMPLEX CZ,C2
       IFLAG1=2
       CZ=C2-P2
C     IF(IFLAG1.EQ.1) WRITE(*,100) C,P,CZ
c 100  FORMAT (' C= ',E16.8,'P= ',2E16.8,'CZ= ',2E16.8)
       U=REAL(CZ)
       X=AIMAG(CZ)
       R=DSQRT(X*X+U*U)
C     IF (IFLAG1.EQ.1) WRITE(*,101) U,X,R
c 101  FORMAT (' U= ',E16.8,'  X= ',E16.8,'  R= ',E16.8)
       W1=DABS(R+U)/2.
       W2=DABS(R-U)/2.
       R1=DSQRT(W1)
       R2=DSQRT(W2)
C     IF (IFLAG1.EQ.1) WRITE(*,102) W1,W2,R1,R2
c 102  FORMAT(1X,' W1=',E16.8,'  W2= ',E16.8,/' R1=',E16.8,'  R2=',E16.8)
C     CR1=R11-R22*(0.,1.)
       CR2=R2+R1*CMPLX(0.,1.)
C     IF (IFLAG1.EQ.1) WRITE(*,103) CR1
       RETURN
       END


       COMPLEX FUNCTION CMM(Z1,Z2,Z3,Z4)
       implicit none
       COMPLEX Z1,Z3,Z4,Z2
       DOUBLE PRECISION X1,X2,Y1,Y2,W1,W2,X3,X4,Y3,Y4,W3,W4
C
C     CMM=Z1*Z2-Z3*Z4
C
       X1=REAL(Z1)
       Y1=AIMAG(Z1)
       X2=REAL(Z2)
       Y2=AIMAG(Z2)
       W1=X1*X2-Y1*Y2
       W2=X1*Y2+X2*Y1
       X3=REAL(Z3)
       Y3=AIMAG(Z3)
       X4=REAL(Z4)
       Y4=AIMAG(Z4)
       W3=X3*X4-Y3*Y4
       W4=X3*Y4+X4*Y3
       W1=W1-W3
       W2=W2-W4
       CMM=W1+W2*(0.,1.)
       RETURN
       END


c       COMPLEX FUNCTION CMINUS(Z1,Z2)
c       COMPLEX Z1,Z2
c       DOUBLE PRECISION X1,X2,Y1,Y2,W1,W2
c       X1=REAL(Z1)
c       Y1=AIMAG(Z1)
c       X2=REAL(Z2)
c       Y2=AIMAG(Z2)
c      W1=X1-X2
c      W2=Y1-Y2
c       W11=W1
c       W22=W2
c       CMINUS=W11+W22*(0.,1.)
c       RETURN
c       END


       COMPLEX FUNCTION CMAD(Z1,Z2,Z3,Z4)
       implicit none
       COMPLEX Z1,Z3,Z4,Z2
       DOUBLE PRECISION X1,X2,Y1,Y2,W1,W2,X3,X4,Y3,Y4,W3,W4
C
C     CMM=Z1*Z2+Z3*Z4
C
       X1=REAL(Z1)
       Y1=AIMAG(Z1)
       X2=REAL(Z2)
       Y2=AIMAG(Z2)
       W1=X1*X2-Y1*Y2
       W2=X1*Y2+X2*Y1
       X3=REAL(Z3)
       Y3=AIMAG(Z3)
       X4=REAL(Z4)
       Y4=AIMAG(Z4)
       W3=X3*X4-Y3*Y4
       W4=X3*Y4+X4*Y3
       W1=W1+W3
       W2=W2+W4
       CMAD=W1+W2*(0.,1.)
       RETURN
       END


c       COMPLEX FUNCTION CADD(Z1,Z2)
c       COMPLEX Z1,Z2
c       DOUBLE PRECISION X1,X2,Y1,Y2,W1,W2
c       X1=REAL(Z1)
c       Y1=AIMAG(Z1)
c       X2=REAL(Z2)
c       Y2=AIMAG(Z2)
c      W1=X1+X2
c      W2=Y1+Y2
c       W11=W1
c       W22=W2
c       CADD=W11+W22*(0.,1.)
c       RETURN
c       END


       COMPLEX FUNCTION CM(Z1,Z2)
       implicit none
       COMPLEX Z1,Z2
       DOUBLE PRECISION X1,X2,Y1,Y2,W1,W2
       X1=REAL(Z1)
       Y1=AIMAG(Z1)
       X2=REAL(Z2)
       Y2=AIMAG(Z2)
       W1=X1*X2-Y1*Y2
       W2=X1*Y2+X2*Y1
       CM=W1+W2*(0.,1.)
       RETURN
       END


       COMPLEX FUNCTION CD(Z1,Z2)
       implicit none
       COMPLEX Z1,Z2
       DOUBLE PRECISION X1,X2,Y1,Y2,W1,W2,W3
       X1=REAL(Z1)
       Y1=AIMAG(Z1)
       X2=REAL(Z2)
       Y2=AIMAG(Z2)
       W3=X2*X2+Y2*Y2
       W1=(X1*X2+Y1*Y2)/W3
       W2=(Y1*X2-Y2*X1)/W3
       CD=W1+W2*(0.,1.)
       RETURN
       END


c      SUBROUTINE CONJ(C,NCENT)
c      DIMENSION C(*)
c      NHM1=NCENT-2
c      NCENT2=NCENT*2
c      DO 1 I=1,NHM1
c      L2=2*I
c      L1=L2-1
c      L3=L2+1
c      K1=NCENT2+L1
c      K2=NCENT2+L2
c      J1=NCENT2-L3
c      J2=NCENT2-L2
c      C(K1)=C(J1)
c1     C(K2)=-C(J2)
c      RETURN
c      END
C

       SUBROUTINE CMODEL(JO,DEPTH,NB)
       implicit none
       integer jo, nb
       real depth
c
c relacionada con el modelo de velocidad
c
       integer inpt, inptd, inpth, nll
       PARAMETER(INPT=2100,INPTD=4200,INPTH=1100,NLL=45)
       integer joo
       real cc, ss, dd, tth, c, s, d, th, qa, qb
       common/smodel/JOO,CC(NLL),SS(NLL),DD(NLL),TTH(NLL)
       COMMON/MODEL/C(NLL),S(NLL),D(NLL),TH(NLL),QA(NLL),QB(NLL)
c
c variables extras
c
       integer j, jj, j1
       real x
c       OPEN(39,FILE=FNAME)
c       READ(39,*) JO
c       READ (39,*) (CC(J),SS(J),DD(J),TTH(J),J=1,JO)
c       CLOSE(39)
       JO=JOO
       X=0.0
       TH(1)=0.0
       DO JJ=1,JO
          X=X+TTH(JJ)
          TH(JJ+1)=X
          C(JJ)=CC(JJ)
          S(JJ)=SS(JJ)
          D(JJ)=DD(JJ)
       enddo
      
       DO J=1,JO
          IF(DEPTH.GT.TH(J).AND.DEPTH.LE.TH(J+1)) THEN
             NB=J+1
             exit
          END IF
       enddo
       DO J=NB,JO
          JJ=JO+NB-J
          J1=JJ+1
          TH(J1)=TH(JJ)
          C(J1)=C(JJ)
          S(J1)=S(JJ)
          D(J1)=D(JJ)
       enddo
       TH(NB)=DEPTH
       if(abs(th(nb) - th(nb+1)) .lt. 1.e-6)th(nb)=th(nb)-0.01
       C(NB)=C(NB-1)+0.0001
       S(NB)=S(NB-1)+0.0001
       D(NB)=D(NB-1)+0.0001
       JO=JO+1
       RETURN
       END


      SUBROUTINE FIND2 (NB,NL,PO,TO,LOVE)
      implicit none
      integer nb, nl, love
      real po, to
      integer inpt, inptd, inpth, nll, nhh
      PARAMETER(INPT=2050,INPTD=4100,INPTH=1025,NLL=45,NHH=2)
      real c, bta, ro, thk, qa, qb
      COMMON/MODEL/C(NLL),BTA(NLL),RO(NLL),THK(NLL),QA(NLL),QB(NLL)
c
c variables extras
c
      integer j
      real yy, xx, e(30)
      TO = 0.0
      DO J = NB,NL-1
         IF(LOVE.EQ.0) THEN
            YY=ABS(1./C(J)+PO)
            XX=ABS(1./C(J)-PO)
         else
            YY=ABS(1./BTA(J)+PO)
            XX=ABS(1./BTA(J)-PO)
         endif
         E(J) =SQRT(XX*YY)
         TO=TO+(THK(J+1)-THK(J))*E(J)
      enddo
C      WRITE (*,17) PO, TO
      RETURN
      END


c      SUBROUTINE SMOOTH(U,NX,NY)
c      DIMENSION U(*)
c      DO 10 IY=1,NY
c      K1=3*(IY-1)*NX+1
c      Y1=U(K1)
c      DO 11 IX=2,NX-1
c      K=K1+3*(IX-1)
c      Y2=U(K)
c      K5=K+3
c      U(K)=0.5*Y2+0.25*(Y1+U(K5))
c11    Y1=Y2
c      K2=K1+3
c      U(K1)=U(K2)*0.4
c      U(K5)=U(K )*0.4
c10    CONTINUE            
c      DO 20 IX=1,NX
c      K1=3*(IX-1)+1
c      Y1=U(K1)
c      DO 21 IY=2,NY-1
c      K=K1+3*(IY-1)*NX
c      Y2=U(K)
c      K5=K1+3*IY*NX
c      U(K)=0.5*Y2+0.25*(Y1+U(K5))
c21    Y1=Y2
c      U(K5)=U(K)*0.1
c20    CONTINUE            
c      RETURN
c      END


      SUBROUTINE INTERP(XP,YP,N,X,Y)
      implicit none
      real xp, yp, x, y
      integer n
      DIMENSION XP(N),YP(N)
c
c variables extras
c
      REAL DIF1,DIF2,DIFY,DR,ratio
      integer i,k
      IF ((X .GT. XP(N)) .or. (X .LT. XP(1))) then
         y = 0.0
         return
      endif
      DO I=1,N
         if ((xp(i) - x) .le. -1.e-6) cycle
         if (abs(xp(i) - x) .lt. 1.e-6) then
            y = yp(i)
         endif
         if ((xp(i) - x) .ge. 1.e-6) then
            k = i-1
            dif1 = xp(i)-xp(k)
            dif2 = xp(i)-x
            ratio = dif2/dif1
            dify = abs(yp(i) - yp(k))
            dr = dify*ratio
            if (yp(i) .gt. yp(k)) then
               y = yp(i)-dr
            else
               y = yp(i)+dr
            endif
         endif
         return
      enddo
c      IF (X.GT.XP(N))  GO TO 6
c      IF (X .LT. XP(1))  GO TO 6
c      DO 10  I=1,N
c         IF (XP(I)-X) 10,102,3
c   10     CONTINUE
c    3      K=I-1
c           DIF1=XP(I)-XP(K)
c           DIF2=XP(I)-X
c           RATIO = DIF2/DIF1
c           DIFY=ABS(YP(I) - YP(K))
c           DR = DIFY*RATIO
c           IF (YP(I) .GT. YP(K))  GO TO 4
c           Y=YP(I) + DR
c           RETURN
c    4      Y=YP(I)-DR
c           RETURN
c  102      Y=YP(I)
c           RETURN
c    6 Y=0.0
c      RETURN
      END


c      SUBROUTINE INTERWSH (U,NU,DTU,F,NF,DTF)
cC--- DTU=0.2
cC--- DTF=0.5
c      DIMENSION U(*),F(*),W(4096)
c      MM = 2*NU
c      DO 1 I=1,NU-1
c	 I1=(I-1)*2+1
c	 I2=2*I
c	 W(I1)=U(I)
c	 W(I2)=(U(I)+U(I+1))*0.5
c    1 CONTINUE
c	 NF=INT(NU*DTU/DTF)-1
c      DO 2 I=1,NF-1
c      J=5*I
c      J10=J-1
c      J11=J+1
c      J21=J+2
c      J20=J-2
c      J30=J-3
c      J31=J+3
c      J41=J+4
c      J40=J-4
c       F(I)=0.20*W(J)+(W(J10)+W(J11))*0.16+(W(J20)+W(J21))*0.12
c       F(I)=F(I)+(W(J30)+W(J31))*0.08+(W(J40)+W(J41))*0.04
c    2 CONTINUE
c	 NF=NF-1
c      RETURN
c      END


        SUBROUTINE INTBHZ(FINST,TVL,GDP,IR,IFLAG,MM,DIRECTORY)
        implicit none
        complex finst(1024)
        real tvl, gdp
        integer iflag, mm, ir
        CHARACTER*200 DIRECTORY
        CHARACTER*12 FRTDK
        integer lnpt, jf, nsyn
        real dt, df, aw
        COMMON/PARAM/DT,LNPT,NSYN,DF,JF,AW,FRTDK
        complex zeroes, poles
        real ff
        integer npoles, nzeroes
        COMMON/PARINT/ff(500),NPOLES(500),NZEROES(500),POLES(100,500)
     *               ,ZEROES(100,500)
c
c variables extras
c
        character*80 pole_inf
        character*200 responsefile
        real rr, ri, pi, w
        COMPLEX S,FLT(1024)
        integer i,ii,irr,jff,jj,nsta 
        IF(IFLAG.EQ.0) THEN
           responsefile=TRIM(DIRECTORY)//'instrumental_response.txt'
           open(14,file=responsefile,status='unknown')
           READ(14,*)NSTA
           DO IRR=1,nsta
c	 READ(14,12)A0(IRR),DS(IRR),DTT,NPOLES(IRR),NZEROES(IRR)
c12	 FORMAT(2(1X,E12.6),F4.1,1X,2(1X,I2))
              read(14,13)pole_inf
!              write(*,13)pole_inf
 13           format(a)
              read(14,*)nzeroes(irr)
              if(nzeroes(irr).gt.0)then
                 do i=1,nzeroes(irr)
                    read(14,*)rr,ri
                    zeroes(i,irr)=cmplx(rr,ri)
                 enddo
              endif
              read(14,*)npoles(irr)
              if(npoles(irr).gt.0)then
                 do i=1,npoles(irr)
                    read(14,*)rr,ri
                    poles(i,irr)=cmplx(rr,ri)
                 enddo
                 read(14,*)ff(irr)            
                 i=nzeroes(irr)
                 if(cabs(zeroes(i,irr)).le.1e-8)then
                    nzeroes(irr)=nzeroes(irr)-1
                 else
                    npoles(irr)=npoles(irr)+1
                    i=npoles(irr)
                    poles(i,irr)=cmplx(1.0e-8,0.0)
                 endif 
              endif
           enddo
           CLOSE(14)
        END IF         
        PI = 4.*ATAN(1.)
c        FF=A0(IR)*DS(IR)

        if(iflag.eq.3)then
           do i=1,jf
              FLT(I)=CMPLX(1.0,0.0)
              FINST(I)=CMPLX(1.,0.)
           enddo
           FINST(1)=CMPLX(0.,0.)
           JFF=JF-1
           IF(TVL.GT.0.01)CALL FTTQ(JFF,DF,TVL,FLT)
           do i=1,jf
              FINST(I)=FINST(I)*FLT(I)*GDP/DT
           enddo
           return
        endif

        DO I=1,JF
           FLT(I)=CMPLX(1.0,0.0)
           FINST(I)=CMPLX(1.,0.)
        enddo
        FINST(1)=CMPLX(0.,0.)
        JFF=JF-1
        IF(TVL.GT.0.01)CALL FTTQ(JFF,DF,TVL,FLT)
        if(nzeroes(ir).le.0 .and. npoles(ir).le.0)then
           do i=1,jf
              finst(i)=finst(i)*flt(i)*GDP/dt
           enddo
        else
           DO I=1,JF
              W=2*PI*(I-1)*DF
              S=CMPLX(0.0,W)
              DO JJ=1,NZEROES(IR)
                 FINST(I)=FINST(I)*(S-ZEROES(JJ,IR))
              enddo
              DO II=1,NPOLES(IR)
                 FINST(I)=FINST(I)/(S-POLES(II,IR))
              enddo
              FINST(I)=FINST(I)*FLT(I)*FF(ir)*GDP/DT
           enddo
        endif
        RETURN
        END


        SUBROUTINE FTTQ(NF,FDF,RQ,FLT)
        implicit none
        integer nf
        real fdf, rq
        complex flt(1024)
C     FUTTERMAN'S Q OPERATOR MODIFIED.  CONTAINS NONE OF THE PHASE
C   SHIFTS AND MADE TO USE THE RATIO T/Q.
C     ND= NUMBER OF POINTS IN FREQUENCY DOMAIN
C     FDF=FREQENCY SEPARATION IN THE FREQENCY DOMAIN
C     RQ=RATIO T/Q
        real fqw
        DIMENSION FQW(2048)
        real a, b, c, d, ck, cr, ep, si, vlg, w0, pi
        integer i, j, nd, ndmi2

C     NUMBER OF POINTS
        ND = 2*NF
C     FUNDEMENTAL FREQUENCY
C     FDF=1.0/((ND-1)*DT)
        A = .25
        C = -3./55.
        D = .991
        pi = 4.0*atan(1.0)
        VLG=A+D*RQ+C*RQ*RQ
C     SAMPLING FREQUENCY
        W0 = 2.0*pi*FDF
        FQW(1)=1.0
        FQW(2)=0.0
        FLT(1)=CMPLX(1.0,0.0)
C
C     THE REAL PART OF FQW IS CONTAINED IN THE ODD NUMBERED LOCATIONS
C     AND THE IMAGINARY PART IS IN THE EVEN LOCATIONS.
C
C     FOLDING FREQUENCY POINT
C     NF=ND/2
        DO I=1,NF
           CK=W0*I*0.5*RQ
           IF (CK.GE.60.) then
              ep = 0.
           else
              EP=EXP(-W0*I*0.5*RQ)
           endif
c           IF (CK.GE.60.) GO TO 30
c           EP=EXP(-W0*I*0.5*RQ)
c           GO TO 31
c 30        EP=0.
c 31        CONTINUE
           A=W0*I*RQ*0.31830989*ALOG(W0*I)
C     B CONTAINS AN ARBITRARY 2 SECOND PHASE SHIFT.
C     B=W0*I*2.0
           B=W0*I*VLG
           CR=COS(A)*COS(B)+SIN(A)*SIN(B)
           SI=SIN(B)*COS(A)-SIN(A)*COS(B)
           FQW(2*I+1)=EP*CR
           FQW(2*I+2)=EP*SI
           FLT(I+1)=CMPLX(FQW(2*I+1),FQW(2*I-1))
           NDMI2=(ND-I)*2+1
           FQW(NDMI2)=FQW(2*I+1)
           NDMI2=NDMI2+1
           FQW(NDMI2)=-FQW(2*I+2)
        enddo
        FQW(ND+2)=0.0
C     MAKING CORRECTION FOR CHANGE IN DEFN'S FOR FOURIER TRANSFORMS.
        DO I=1,ND
           FQW(2*I)=-FQW(2*I)
        enddo
        J=0
        DO I=1,ND,2
           J=J+1
           FLT(J)=CMPLX(FQW(I),FQW(I+1))
        enddo      
        RETURN
        END
