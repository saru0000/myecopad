    program TECO_MCMC
    implicit none 
!#
!#
!#Short names:
!#ORC - ORiginal Code
!#NC - New Code
!#
!
! USE IFPORT !! Was: USE DFPORT
!   for parameter file
    real lat,longi,wsmax,wsmin                 
    real LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax
    real SapR,SapS,SLA,GLmax,GRmax,Gsmax,stom_n
    real a1,Ds0,Vcmax0,extkU,xfang,alpha
    real Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C
    real shcap_s,condu_s,shcap_snow,condu_snow,albedo_snow,resht,thd_snow_depth,b_bound
    real infilt_rate
    real Tau_Micro,Tau_slowSOM,Tau_Passive
    real gddonset,Q10,Rl0,Rs0,Rr0
    real r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi       !parameters added for MEMCMC
    real obs_CH4_MEMCMC(2,1000),std_CH4(2,1000)
    
    real fa,fsub,rho_snow,decay_m
    character(len=50) parafile
    
!   for climate file
    integer, parameter :: ilines=90000
    integer, parameter :: iiterms=7
    integer,dimension(ilines):: year_seq,doy_seq,hour_seq
    real forcing_data(iiterms,ilines)
    character(len=50) climatefile
    integer yr_length,lines
    integer,dimension(ilines):: year_seq2,doy_seq2,hour_seq2
    integer yr_length2,lines2
    
    ! water table 
    real water_table(ilines),snow_in(ilines)
    character(len=50) watertablefile,snowdepthfile
    
!   for observation files    
    real obs_spruce(12,1000),std(12,1000),obs_soilwater(5,40000),obs_soilt(11,40000)
    real obs_watertable_d(2,4000),obs_soilt_d(11,4000)
    real treatment
    character(len=80) obsfile1,obsfile2,covfile,obsfile3,obsfile4,obsfile5,obsfile6
    integer len1,len2,len3,len6
    
!   for MCMC
    integer,parameter:: MCMC=0 ! 0:run model simulation; 1:data assimilation 
    integer IDUM,upgraded,isimu,tmp_up
    integer, parameter :: npara=5       ! Number of parameters to be estimated
    real search_length
    real J_last
    real Simu_dailyflux(12,80000),Simu_soilwater(10,40000),Simu_soiltemp(11,40000)
    real Simu_dailysoilt(11,80000),Simu_dailywatertable(1,80000),Simu_dailyice(10,80000)
    real Simu_watertable(1,40000),Simu_snowdepth(1,80000)
    real Simu_dailyCH4(6,80000)
    real coef(npara),coefac(npara),coefnorm(npara)
    real coefmax(npara),coefmin(npara)
    real r,fact_rejet
    real gamma(npara,npara),gamnew(npara,npara)     ! covariance matrix
    integer k1,k2,rejet,paraflag,k3,i,reject
    integer, parameter :: nc=100
    integer, parameter :: ncov=500
!! nc: the multiplicative constant will be adjusted !!
!!        every nc iterations, to preserve an adequate!!
!!        acceptation rate   
!! ncov: the covariance matrix gamma will be updated  !!
!!        every ncov iterations		
    real coefhistory(ncov,npara)
    character(len=80) outfile
    integer new,covexist,num,j
    
!   for consts parameteres
    real,dimension(3):: tauL,rhoL,rhoS
    real pi,emleaf,emsoil,Rconst,sigma,cpair,Patm,Trefk
    real H2OLv0,airMa,H2OMw,chi,Dheat,wleaf,gsw0,eJmx0,theta
    real conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,Edjm
    real Entrpy,gam0,gam1,gam2
    
!   for initialize
    real fwsoil,topfws,omega,Storage,nsc
    real fwsoil_initial,topfws_initial,omega_initial
    real Storage_initial,nsc_initial
    real wcl(10),QC(8)
    real wcl_initial(10),QC_initial(8)
    integer yrs_eq
    real diff_snow,diff_s,condu_b
    real depth_ex
    
!!!!!!!!!!!! user supplied configuration !!!!!!!!!!!!!!!!!!!!    
    logical, parameter :: do_soilphy     = .True.  ! if true, model simulate layered soil temperature 
!    logical, parameter :: do_watertable   = .True. ! if true, model simulate watertable; else read snow from observation files
    logical, parameter :: do_snow   = .TRUE. ! if true, model simulate snow; else read snow from observation files
    
    logical, parameter :: do_soilt_da     = .False.
    logical, parameter :: do_snow_da     = .False.
    logical, parameter :: do_watertable_da     = .False.
    logical, parameter :: do_methane_da   = .True.
    logical, parameter :: do_co2_da   = .False.
    
!   Read parameters from file
!ORC
 !   parafile='./input/f0_SPRUCE_parameterfile.txt'
!ORC
!NC
call getarg(1,parafile)
!NC
    call Getparameters(lat,longi,wsmax,wsmin,                   &              
    &   LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax,                   &
    &   SapR,SapS,SLA,GLmax,GRmax,Gsmax,stom_n,                 &
    &   a1,Ds0,Vcmax0,extkU,xfang,alpha,                        &
    &   Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C,                 &
    &   Tau_Micro,Tau_slowSOM,Tau_Passive,                      &
    &   gddonset,Q10,RL0,Rs0,Rr0,parafile,                      &
    &   r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi)    !added for MEMCMC
    
! Parameters for soil physical part Yuanyuan 
    shcap_snow=690000.  ! refer to csm 4.0 physical constant ! 2090 J/kg/°K with density of 330 kg/m3
    condu_snow=0.1
    condu_b = 0.2
    depth_ex=0.05
    
!    shcap_snow=700246.3125  ! refer to csm 4.0 physical constant ! 2090 J/kg/°K with density of 330 kg/m3
!    condu_snow=0.0202
!    condu_b = 0.0797
!    depth_ex=0.0472
    
    diff_s=1.
    diff_snow =1.8
    albedo_snow=0.7
    resht=40.
    thd_snow_depth=4.0
    b_bound=100.
    
    infilt_rate= 0.0017
    fa = 1
    fsub=0.1
    rho_snow=100.
    decay_m=2.2192
    
!   Read climatic forcing
call getarg(2,climatefile)
!    climatefile='./input/f1_EM1forcing2011-2015.txt'
!call getarg(2,climatefile)
!COMMENT Recive ERROR from getarg(3,w..)  At line 4371 of file /source/frc_01.f90 (unit = 111, file = 'fort.111')
call getarg(3,watertablefile)
!    watertablefile='./input/f2_SPRUCE_Water_Table_Level_2011-2014.txt'
call getarg(4,snowdepthfile)
!    snowdepthfile='./input/f3_SPRUCE_Snow_Depth_2011-2014.txt'
       
    call Getclimate(year_seq,doy_seq,hour_seq,          &
    &   forcing_data,climatefile,lines,yr_length)
    call Getwatertable(year_seq,doy_seq,hour_seq,          &
    &   water_table,watertablefile,lines,yr_length)
    call Getsnowdepth(year_seq,doy_seq,hour_seq,          &
    &   snow_in,snowdepthfile,lines,yr_length)
    
!   Read observation data
    obsfile1='./input/f4_obs_SPRUCE.txt'
    obsfile2='./input/f5_SPRUCE_hummock_toplayer.txt'
    obsfile3='./input/f6_SPRUCE_soilt.txt'
    obsfile4='./input/f7_daily_soilt_2011-2014.txt'
    obsfile5='./input/f8_daily_watertable_2011-2014.txt'
    obsfile6='./input/f9_obs_CH4_for_MEMCMC.txt'
    treatment=0.    ! Ambient temperature
    call GetObsData(obs_spruce,std,obs_soilwater,     &
            &   len1,len2,obsfile1,obsfile2,obs_soilt,obsfile3,obs_watertable_d,obs_soilt_d,&
            &   obsfile4,obsfile5,len3)  
            !***************************************************************
!   this region is modified for MEMCMC      Shuang 07/18/2016
    call GetObsCH4Data(obs_CH4_MEMCMC,std_CH4,len6,obsfile6)	!all are outputs
!***************************************************************    
        
!   initiations for canopy model, including canopy traits variation in a year
    call consts(pi,tauL,rhoL,rhoS,emleaf,emsoil,&
     &    Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat,&
     &    wleaf,gsw0,Vcmax0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci,&
     &    Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2)
    fwsoil=1.0
    topfws=1.0
    omega=1.0
    do i=1,10
        wcl(i)=wsmax/100.
    enddo 
    Storage=32.09           !g C/m2
    nsc=85.35
!         put the values into a matrix
!   QC=(/300.,6300.,300.,119.,300.,322.,8834.,312./)
!   QC=(/100.,800.,100.,39.,100.,122.,834.,12./)
    QC=(/440.,700.,300.,119.,300.,42.,38340.,23120./)
    !QC=(/500.,700.,300.,209.,300.,82.,38340.,23120./)  !changed for a better Rh initial value Jun.30th Shuang
       
!   Start main loop
    open(61,file='./output/SPRUCE_yearly.txt')
    !open(62,file='Simu_dailyflux.txt')
    !open(63,file='Simu_soilwater.txt')
    !open(633,file='Simu_soiltemp.txt')
    !open(81,file='TECO_ME_water_table0.95.csv')
    !open(82,file='TECO_ME_Methane.csv')
    !open(83,file='TECO_output.csv')
    !open(664,file='Simu_dailyCH4.txt')
    open(88,file='./output/snow.csv')
    write(88,*) "melt,snow_dsim,snow_in,ta"
    
    !write(81,*)'vtot,zwt,wsc(1),wsc(2),wsc(3),wsc(4),wsc(5),wsc(6),wsc(7),wsc(8),wsc(9),wsc(10)'  
    !write(82,*)'zwt,Ebu_sum_sat,Ebu_sum_unsat,Pla_sum,simuCH4,  &
    !      & Rh(1),Rh(2),Rh(3),Rh(4),Rh(5),Rh(6),Rh(7),Rh(8),Rh(9),Rh(10),   &
    !      & consum,Pro_sum, &
    !      & ProCH4(1),ProCH4(2),ProCH4(3),ProCH4(4),ProCH4(5),ProCH4(6),ProCH4(7),ProCH4(8),ProCH4(9),ProCH4(10),   &
    !      & CH4(1),CH4(2),CH4(3),CH4(4),CH4(5),CH4(6),CH4(7),CH4(8),CH4(9),CH4(10), &
    !      & Fdifu(1),Fdifu(2),Fdifu(3),Fdifu(4),Fdifu(5),Fdifu(6),Fdifu(7),Fdifu(8),Fdifu(9),Fdifu(10), &
    !      & OxiCH4(1),OxiCH4(2),OxiCH4(3),OxiCH4(4),OxiCH4(5),OxiCH4(6),OxiCH4(7),OxiCH4(8),OxiCH4(9),OxiCH4(10),   &
    !      & wsc(1),wsc(2),wsc(3),wsc(4),wsc(5),wsc(6),wsc(7),wsc(8),wsc(9),wsc(10), &
    !      & Rh_pools(1),Rh_pools(2),Rh_pools(3),Rh_pools(4),Rh_pools(5)'
    !write(83,*)'zwt,Tsoil,Rh_pools(1),Rh_pools(2),Rh_pools(3),Rh_pools(4),Rh_pools(5),&
    !      & wsc(1),wsc(2),wsc(3),wsc(4),wsc(5),wsc(6),wsc(7),wsc(8),wsc(9),wsc(10)'
           
    open(71,file='./output/Paraest.txt')
    
    if(MCMC.eq.1) GOTO 100
    yrs_eq=yr_length*32  ! spin up length 
    num=0
    call TECO_simu(MCMC,Simu_dailyflux,Simu_soilwater,      &
     &        obs_spruce,obs_soilwater,              &
     &        forcing_data,yr_length,year_seq,doy_seq,hour_seq,lines,   &
     &        fwsoil,topfws,omega,wcl,Storage,nsc,yrs_eq,QC,    &
     &        lat,longi,wsmax,wsmin,LAIMAX,LAIMIN,rdepth,     &
     &        Rootmax,Stemmax,SapR,SapS,SLA,GLmax,GRmax,Gsmax,      &
     &        stom_n,a1,Ds0,Vcmax0,extkU,xfang,alpha,               &
     &        tau_Leaf,tau_Wood,tau_Root,tau_F,tau_C,tau_Micro,     &   ! the unit is year
     &        tau_SlowSOM,tau_Passive,gddonset,                     &
     &        Q10,Rl0,Rs0,Rr0,pi,tauL,rhoL,rhoS,emleaf,emsoil,&
     &    Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat,&
     &    wleaf,gsw0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci,&
     &    Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,Simu_soiltemp,water_table,snow_in,&
     &    obs_soilt,diff_s,diff_snow,albedo_snow,resht,thd_snow_depth,b_bound,&
     &    Simu_watertable,infilt_rate,Simu_dailysoilt,Simu_dailywatertable,&
     &    Simu_snowdepth,fa,fsub,rho_snow,decay_m,Simu_dailyice,shcap_snow,condu_snow,condu_b,&
     &    depth_ex,r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi,Simu_dailyCH4,do_snow,do_soilphy)    
    write(*,*)'run simulation'
    return
100 continue

!   spin up
!    yrs_eq=yr_length*50  ! spin up length 
!    call TECO_simu(MCMC,Simu_dailyflux,Simu_soilwater,                &
!     &        obs_spruce,obs_soilwater,len1,len2,              &
!     &        forcing_data,yr_length,year_seq,doy_seq,hour_seq,lines,   &
!     &        fwsoil,topfws,omega,wcl,Storage,nsc,yrs_eq,QC,    &
!     &        lat,longi,wsmax,wsmin,LAIMAX,LAIMIN,rdepth,     &
!     &        Rootmax,Stemmax,SapR,SapS,SLA,GLmax,GRmax,Gsmax,      &
!     &        stom_n,a1,Ds0,Vcmax0,extkU,xfang,alpha,               &
!     &        tau_Leaf,tau_Wood,tau_Root,tau_F,tau_C,tau_Micro,     &   ! the unit is year
!     &        tau_SlowSOM,tau_Passive,gddonset,                     &
!     &        Q10,Rl0,Rs0,Rr0,pi,tauL,rhoL,rhoS,emleaf,emsoil,&
!     &    Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat,&
!     &    wleaf,gsw0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci,&
!     &    Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2)    
    fwsoil_initial = fwsoil
    topfws_initial = topfws
    omega_initial = omega
    wcl_initial = wcl
    Storage_initial = Storage
    nsc_initial = nsc
    QC_initial = QC
!   end of spin up
       
!	simulation begins here
    J_last=9000000.0
    IDUM = 542
    upgraded=0
    new=0
    k3=0
    
if (do_co2_da) then     
!    coef(1)=Tau_Leaf
!    coef(2)=Tau_Wood
!    coef(3)=Tau_Root
!    coef(4)=Tau_F
!    coef(5)=Tau_C
!    coef(6)=Tau_Micro
!    coef(7)=Tau_slowSOM
!    coef(8)=Tau_Passive
!    coef(9)=GLmax
!    coef(10)=GRmax
!    coef(11)=Gsmax
!    coef(12)=Vcmax0
!    coef(13)=RL0
!    coef(14)=Rs0
!    coef(15)=Rr0
!    coef(16)=Q10
!    coef(17)=SLA
!    coef(18)=gddonset
!    coef(19)=diff_s
!    coef(20)=diff_snow
!    coef(21)=albedo_snow
!    coef(22)=resht
!    coef(23)=thd_snow_depth
!    coef(24)=b_bound
!    coef(25)=infilt_rate
!!
!    coefmax=(/3.0, 800.0, 2.0, 0.5, 20.0, 0.5, 1000.0, 4000.0,   &
!        &   50.0, 30.0, 30.0, 180.0, 45.0, 10.5, 45.0, 4.0, 200., 160., &
!        &   587.,3600.,1.0,400.,25.,1000.,0.02/)
!    coefmin=(/0.5, 5.0, 0.3, 0.1, 1.0, 0.05, 5.0, 500.0,    &
!        &   10.0, 10.0, 10.0, 14.0, 10.0, 4.5,  10.0, 1.0, 10., 100., &
!        &   0.01,0.01,0.5,40.,0.,20.,0.0001/)
        
endif        
   

if (do_watertable_da) then     
!    coef(1)=infilt_rate
!    coefmax=(/0.02/)
!    coefmin=(/0.0001/)
endif     

!!!!!!!!!!!!! soil physical DA !!!!!!!!!!!!!!!!!    
if (do_soilt_da) then 
!    coef(1)=shcap_snow
!    coef(2)=condu_snow
!    coef(3)=condu_b
!    coef(4)=depth_ex
!    coefmax=(/1150000.,0.56,0.23,0.1/)
!    coefmin=(/20909.,0.01,0.023,0./)  
    
!    coef(1)=shcap_s
!    coef(2)=condu_s
!    coef(3)=shcap_snow
!    coef(4)=condu_snow
!    coef(5)=albedo_snow
!    coef(6)=resht
!    coef(7)=thd_snow_depth
!    coef(8)=b_bound
!       
!    coefmax=(/4200000.0,2.29,1150000.0,0.56,1.0,400.,25.,1000./)
!    coefmin=(/1225.,0.02,20909.,0.02,0.5,40.,0.,20./)    
!     
    ! shcap_snow with density 10 kg/m3 550 kg/m3
!    coef(1)=diff_s
!    coef(2)=diff_snow
!    coef(3)=infilt_rate
!    coefmax=(/5.,5.,0.03/)
!    coefmin=(/0.1,0.01,0.001/)
    
endif
    
!!!!!!!! methane DA !!!!!!!!!!!!!!!!!!!!!!!!!! 
if (do_methane_da) then   
    coef(1)=r_me                                                !para in methane production
    coef(2)=Q10pro    
    coef(3)=Omax                                            	!para in methane oxidation
    coef(4)=Tveg  
    coef(5)=CH4_thre
    coefmax=(/0.7, 16.0, 40.0,2.0,500.0/)             !5 paras          
    coefmin=(/0.01, 0.1, 10.0,0.001,350.0/)  
endif

if (do_snow_da) then
!    coef(1)=fa
!    coef(2)=fsub
!    coef(3)=rho_snow
!    coef(4)=decay_m
!    coefmax=(/1.,1.,500.,10./)
!    coefmin=(/0.1,0.1,100.,0.1/)
endif
        
    ! initialize covariance matrix    
    covexist=0
    if(covexist.eq.1)then      ! If prior covariance exists, read from file
        covfile='covariance.txt'
        call getCov(gamma,covfile,npara)

        call racine_mat(gamma,gamnew,npara)      ! square root of covariance matrix
        gamma=gamnew
        do k1=1,npara
!           coefnorm(k1)=(coef(k1)-coefmin(k1))/(coefmax(k1)-coefmin(k1))
            coefnorm(k1)=0.5
            coefac(k1)=coefnorm(k1)
        enddo
    else
        coefac=coef
    endif
 
    fact_rejet=2.4/sqrt(real(npara))
    search_length=0.05
    rejet = 0
    
    do isimu=1,10000
!       generate parameters
        if(covexist.eq.1)then 
            paraflag=1
            do while(paraflag.gt.0)
                call gengaussvect(fact_rejet*gamma,coefac,coefnorm,npara)
                paraflag=0
                do k1=1,npara
                    if(coefnorm(k1).lt.0. .or. coefnorm(k1).gt.1.)then
                    paraflag=paraflag+1
                    write(*,*)'out of range',paraflag
                    endif
                enddo
            enddo
            do k1=1,npara
                coef(k1)=coefmin(k1)+coefnorm(k1)*(coefmax(k1)-coefmin(k1))
            enddo
        else
            call coefgenerate(coefac,coefmax,coefmin,coef,search_length,npara)
        endif

!         update parameters
      if (do_co2_da) then        
!        Tau_Leaf=coef(1)
!        Tau_Wood=coef(2)
!        Tau_Root=coef(3)
!        Tau_F=coef(4)
!        Tau_C=coef(5)
!        Tau_Micro=coef(6)
!        Tau_slowSOM=coef(7)
!        Tau_Passive=coef(8)
!        GLmax=coef(9)
!        GRmax=coef(10)
!        Gsmax=coef(11)
!        Vcmax0=coef(12)
!        RL0=coef(13)
!        Rs0=coef(14)
!        Rr0=coef(15)
!        Q10=coef(16)
!        SLA=coef(17)
!        gddonset=coef(18)
!        
!        diff_s = coef(19)
!        diff_snow = coef(20)
!        albedo_snow=coef(21)
!        resht=coef(22)
!        thd_snow_depth=coef(23)
!        b_bound=coef(24)
!        infilt_rate=coef(25)
      endif
     
      if (do_soilt_da) then 
!        shcap_snow=coef(1)
!        condu_snow=coef(2)
!        condu_b=coef(3)
!        depth_ex=coef(4)
      endif

!!!!!!!! CH4 DA !!!!!!!!!!!!!!!!!!!!!!!!!!
     if (do_methane_da) then     
        r_me=coef(1)                                                !para in methane production
        Q10pro=coef(2)    
        Omax=coef(3)                                            	!para in methane oxidation
        Tveg=coef(4)  
        CH4_thre=coef(5)
     endif

    if (do_watertable_da) then     
        !infilt_rate=coef(1)
    endif  

    if (do_snow_da) then
!        fa=coef(1)
!        fsub=coef(2)
!        rho_snow=coef(3)
!        decay_m=coef(4)
    endif
              
        yrs_eq = 0
        call TECO_simu(MCMC,Simu_dailyflux,Simu_soilwater,            &
        &        obs_spruce,obs_soilwater,              &
        &        forcing_data,yr_length,year_seq,doy_seq,hour_seq,lines,   &
        &        fwsoil,topfws,omega,wcl,Storage,nsc,yrs_eq,QC,         &
        &        lat,longi,wsmax,wsmin,LAIMAX,LAIMIN,rdepth,     &
        &        Rootmax,Stemmax,SapR,SapS,SLA,GLmax,GRmax,Gsmax,      &
        &        stom_n,a1,Ds0,Vcmax0,extkU,xfang,alpha,               &
        &        tau_Leaf,tau_Wood,tau_Root,tau_F,tau_C,tau_Micro,     &   ! the unit is year
        &        tau_SlowSOM,tau_Passive,gddonset,                     &
        &        Q10,Rl0,Rs0,Rr0,pi,tauL,rhoL,rhoS,emleaf,emsoil,&
     &    Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat,&
     &    wleaf,gsw0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci,&
     &    Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,Simu_soiltemp,water_table,snow_in,&
     &    obs_soilt,diff_s,diff_snow,albedo_snow,resht,thd_snow_depth,b_bound,&
     &    Simu_watertable,infilt_rate,Simu_dailysoilt,Simu_dailywatertable,&
     &    Simu_snowdepth,fa,fsub,rho_snow,decay_m,Simu_dailyice,shcap_snow,condu_snow,condu_b,&
     &    depth_ex,r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi,Simu_dailyCH4,do_snow,do_soilphy)      

!         M-H algorithm
        tmp_up=upgraded
        
        call costFObsNee(Simu_dailyflux,Simu_soilwater,   &
            &   obs_spruce,std,obs_soilwater,len1,len2,   &
            &   J_last,upgraded,Simu_soiltemp,obs_soilt,Simu_watertable,water_table,&
            &   Simu_dailysoilt,Simu_dailywatertable,obs_soilt_d,obs_watertable_d,len3,Simu_dailyice, &
            &   Simu_dailyCH4,obs_CH4_MEMCMC,std_CH4,len6,do_soilt_da,do_snow_da,do_watertable_da, &
            &   do_methane_da,do_co2_da)

        if(upgraded.gt.tmp_up)then 
            new=new+1
            if(covexist.eq.1)then
                coefac=coefnorm
                coefhistory(new,:)=coefnorm
            else
                coefac=coef
                do k1=1,npara
                    coefnorm(k1)=(coef(k1)-coefmin(k1))/(coefmax(k1)-coefmin(k1))
                enddo
            endif
            coefhistory(new,:)=coefnorm
            if(new.ge.ncov)new=0
            write(71,701)upgraded,(coef(i),i=1,npara)
            if(upgraded.gt.1000 .and. k3.lt.500)then
                CALL random_number(r)
                if(r.gt.0.95)then
                    k3=k3+1
                    write(outfile,"(A21,I3.3,A4)") "output/Simu_dailyflux",k3,".txt"
                    outfile=trim(outfile)
                    open(62,file=outfile)
                    do i=1,1461
                        write(62,602)i,(Simu_dailyflux(j,i),j=1,12)
                    enddo
602                 format((i7),",",12(f15.4,","))
                    close(62)
                    
                    write(outfile,"(A21,I3.3,A4)") "output/Simu_dailyCH4",k3,".txt"
                    outfile=trim(outfile)
                    open(664,file=outfile)

                    do i=1,1461
                        write(664,604)i,(Simu_dailyCH4(j,i),j=1,6)
                    enddo
604                 format((i7),",",6(f15.4,","))
                    
                    write(outfile,"(A21,I3.3,A4)") "output/Simu_soilwater",k3,".txt"
                    outfile=trim(outfile)
                    open(63,file=outfile)
                    do i=1,35040
                        write(63,603)i,(Simu_soilwater(j,i),j=1,10)
                    enddo
603                 format((i7),",",10(f15.4,","))
                    close(63)
                    
                    write(outfile,"(A26,I3.3,A4)") "output/Simu_soiltemp_daily",k3,".txt"
                    outfile=trim(outfile)
                    write(*,*)'outfile',outfile
                    open(633,file=outfile)
                    do i=1,1460
                        write(633,613)i,(Simu_dailysoilt(j,i),j=1,11)
                    enddo
                    close(633)
613                 format((i7),",",11(f15.4,","))
                    
                    write(outfile,"(A21,I3.3,A4)") "output/Simu_ice_daily",k3,".txt"
                    outfile=trim(outfile)
                    write(*,*)'outfile',outfile
                    open(6633,file=outfile)
                    do i=1,1460
                        write(6633,603)i,(Simu_dailyice(j,i),j=1,10)
                    enddo
                    close(6633)
                    					
		   write(outfile,"(A27,I3.3,A4)") "output/Simu_dailywatertable",k3,".txt"
                    outfile=trim(outfile)
                    write(*,*)'outfile',outfile
                    open(6333,file=outfile)
                    do i=1,1460
                        write(6333,603)i,(Simu_dailywatertable(j,i),j=1,1)
                    enddo
                    close(6333)
                    
                   write(outfile,"(A21,I3.3,A4)") "output/Simu_snowdepth",k3,".txt"
                    outfile=trim(outfile)
                    write(*,*)'outfile',outfile
                    open(63333,file=outfile)
                    do i=1,1461
                        write(63333,603)i,(Simu_snowdepth(j,i),j=1,1)
                    enddo
                    close(63333)
                                        
                endif
            endif
        else
            reject=reject+1
        endif

        fwsoil = fwsoil_initial
        topfws = topfws_initial
        omega = omega_initial
        wcl = wcl_initial
        Storage=Storage_initial
        nsc=nsc_initial
        QC = QC_initial
701     format(I12,",",25(F15.4,","))           
        write(*,*)isimu,upgraded
    
    	! updates of the multiplicative constant
        if(covexist.eq.1)then
            if(mod(isimu,nc).eq.0)then
                if ((1. - real(rejet)/real(nc)) < 0.23) then
                !    fact_rejet = fact_rejet*0.9
                else
                    if ((1. - real(rejet)/real(nc)) > 0.44) then
                !    fact_rejet = fact_rejet * 1.1
                    endif
                endif
            rejet=0
            write(*,*)'search length is', search_length
            endif
        else
            if(mod(isimu,nc).eq.0)then
                if(real(upgraded)/real(isimu) .lt. 0.23)then
                        search_length=search_length*0.9
                else
                    if(real(upgraded)/real(isimu) .gt. 0.44)then
                        search_length=search_length*1.1
                    endif
                endif
                reject=0
                write(*,*)'search length is', search_length
            endif
        endif

	! updates of the covariance matrix
        if(covexist.eq.0 .and. mod(upgraded,ncov).eq.0)then
            covexist=1       !  Yuanyuan
            coefac=coefnorm
            call varcov(coefhistory,gamnew,npara,ncov)
            if (.not.(all(gamnew==0.))) then
                gamma=gamnew
                call racine_mat(gamma,gamnew,npara)
                gamma=gamnew
            endif
        endif
	if (mod(upgraded,ncov).eq.0 .and. covexist.eq.1) then
            call varcov(coefhistory,gamnew,npara,ncov)
            if (.not.(all(gamnew==0.))) then
                gamma=gamnew
                call racine_mat(gamma,gamnew,npara)
                gamma=gamnew
            endif
	endif

    enddo !isimu

    open(72,file='covvariance_temp.txt')
    do i=1,npara
        write(72,*) (gamma(j,i),j=1,npara)
    enddo
    close(72)    
    close(61)

    close(71)
    write(*,*)'run MCMC'
    end

! ====================================================================
    subroutine TECO_simu(MCMC,Simu_dailyflux,Simu_soilwater,      &
     &        obs_spruce,obs_soilwater,              &
     &        forcing_data,yr_length,year_seq,doy_seq,hour_seq,lines,   &
     &        fwsoil,topfws,omega,wcl,Storage,nsc,yrs_eq,QC,            &
     &        lat,longi,wsmax,wsmin,LAIMAX,LAIMIN,rdepth,     &
     &        Rootmax,Stemmax,SapR,SapS,SLAx,GLmx,GRmx,Gsmx,      &
     &        stom_n,a1,Ds0,Vcmax0,extkU,xfang,alpha,               &
     &        tau_L,tau_W,tau_R,tau_F,tau_C,tau_Micr,     &   ! the unit is year
     &        tau_Slow,tau_Pass,gddonset,                     &
     &        Q10,Rl0,Rs0,Rr0,pi,tauL,rhoL,rhoS,emleaf,emsoil,&
     &    Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat,&
     &    wleaf,gsw0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci,&
     &    Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,Simu_soiltemp,water_table,snow_in,&
     &    obs_soilt,diff_s,diff_snow,albedo_snow,resht,thd_snow_depth,b_bound,&
     &    Simu_watertable,infilt_rate,Simu_dailysoilt,Simu_dailywatertable,&
     &    Simu_snowdepth,fa,fsub,rho_snow,decay_m,Simu_dailyice,shcap_snow,condu_snow,condu_b,&
     &    depth_ex,r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi,Simu_dailyCH4,do_snow,do_soilphy)
      
      implicit none    
   !==================== test variables
    real fwsoil_yr,omega_yr,topfws_yr, difference,diff_yr,diff_d
            
      integer, parameter :: iiterms=7            ! 9 for Duke forest FACE
      integer, parameter :: ilines=90000         ! the maxmum records of Duke Face, 1998~2007
      real, parameter:: times_storage_use=720.   ! 720 hours, 30 days
      integer  lines,idays,MCMC
      integer,dimension(ilines):: year_seq,doy_seq,hour_seq
      real forcing_data(iiterms,ilines),input_data(iiterms,ilines)
      real water_table(ilines),snow_in(ilines)
      real Simu_dailyflux(12,80000),Simu_soilwater(10,40000),Simu_soiltemp(11,40000)
      real Simu_watertable(1,40000),Simu_dailysoilt(11,80000),Simu_dailywatertable(1,80000)
      real Simu_dailyice(10,80000)
      real Simu_snowdepth(1,80000)
      real obs_spruce(12,1000),obs_soilwater(5,40000),obs_soilt(11,40000)
      integer pheno,phenoset,day_mod,num
!      site specific parameters
      real lat,longi,rdepth,LAIMAX,LAIMIN
      real wsmax,wsmin,co2ca
      real tau_L,tau_W,tau_R
      real tau_F,tau_C,tau_Micr,tau_Slow,tau_Pass
      real TauC(8)
!      the variables that should be initialized in the begining
      real Q_soil
      real QC(8) !  leaf,wood,root,fine lit.,coarse lit.,Micr,Slow,Pass
      real Pool1,Pool2,Pool3,Pool4,Pool5,Pool6,Pool7,Pool8
      real out1_yr,out2_yr,out3_yr,out4_yr,out5_yr,out6_yr,out7_yr,out8_yr
      real OutC(8)
      real Rh_pools(5)
!      for soil conditions
      real WILTPT,FILDCP,infilt
      real Rsoilabs
      real fwsoil,topfws,omega
!      for plant growth and allocation
      real NSC,NSCmin,NSCmax,add               ! none structural carbon pool
      real Growth,Groot,Gshoot,GRmax           ! growth rate of plant,root,shoot,and max of root
      real St,Sw,Ss,Sn,Srs,Sps,fnsc,Weight     ! scaling factors for growth
!      variables for canopy model
      real evap,transp,ET
      real wind,eairp,esat,rnet
      real Pa_air 
      real gpp,gpp_ra,NPP,NEE,NEP,gpp_d,NPP_d
      real evap_d,transp_d
      real,dimension(3):: tauL,rhoL,rhoS,reffbm,reffdf,extkbm,extkdm
      real,dimension(2):: Radabv
      real Qcan(3,2)
!      parameters for photosynthesis model
      real stom_n,a1,Ds0,Vcmx0,Vcmax0,extkU,xfang,alpha
      real pi,emleaf,emsoil
      real Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat
      real wleaf,gsw0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci
      real Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2
!     for nitrogen sub-model
      real CNmin,CNmax,NSNmax,NSNmin
      real NSN
!      QNleaf,QNwood,QNroot,QNfine,QNcoarse,QNmicr,QNslow,QNpass
!      CN_leaf,CN_wood,CN_root,CN_fine,CN_coarse,CN_micr,CN_slowC,CN_pass
      real QN(8),CN0(8),CN(8),OutN(8),QNplant,QNminer
      real N_leaf,N_wood,N_root,N_deficit
      real N_LF,N_WF,N_RF
      real N_uptake,N_leach,N_vol,N_fixation,N_deposit,N_fert
      real N_up_d,N_fix_d,N_dep_d,N_leach_d,N_vol_d
      real N_up_yr,N_fix_yr,N_dep_yr,N_leach_yr,N_vol_yr
      real N_miner,alphaN
      real SNvcmax,SNgrowth,SNRauto,SNrs
      real diff_s,diff_snow,albedo_snow,resht,thd_snow_depth,shcap_snow,condu_snow,depth_ex
      real infilt_rate
      real b_bound,fa,fsub,rho_snow,decay_m,condu_b

!      additional arrays to allow output of info for each layer
      real,dimension(5):: RnStL,QcanL,RcanL,AcanL,EcanL,HcanL
      real,dimension(5):: GbwcL,GswcL,hG,hIL
      real,dimension(5):: Gaussx,Gaussw,Gaussw_cum 
!      for phenology
      real LAI,bmroot,bmstem,bmleaf,bmplant,totlivbiom,ht
      real SLA,SLAx,L_fall,L_add,litter,seeds
      real GDDonset,GDD5,accumulation,storage,stor_use,store
      real RaL,RaS,RaR  !allocation to respiration
      real alpha_L,alpha_W,alpha_R ! allocation ratio to Leaf, stem, and Root
      real Q10,Rl0,Rs0,Rr0         ! parameters for auto respiration
      real Rgrowth,Rnitrogen,Rmain,Rauto !respirations
      real RmLeaf,RmStem,RmRoot          ! maintanence respiration
      real RgLeaf,RgStem,RgRoot          ! growth respiration
      real RaLeaf,RaStem,RaRoot
      real Rsoil,Rhetero,Rtotal
      real Ra_Nfix,Rh_Nfix
      real gpp_yr,NPP_yr,NEE_yr,RaL_yr,RaR_yr,RaS_yr,Rh_yr
      real Rh4_yr,Rh5_yr,Rh6_yr,Rh7_yr,Rh8_yr,Ra_yr
      real R_Ntr_yr
      real NPPL_yr,NPPR_yr,NPPS_yr,NPP_L,NPP_R,NPP_W
      real Rootmax,Stemmax,SapS,SapR,StemSap,RootSap
      REAL ws,wdepth
!      climate variables for every day
      real Ta,Tair,Ts,Tsoil,water_table_depth,snow_depth
      real doy,hour,Dair,Rh,radsol
      real PAR
!      output daily means of driving variables
      real CO2air_d_avg,SWdown_d_avg,Psurf_d_avg
      real Rain_d_avg,Tair_d_avg,Wind_d_avg
      real Simu_dailyCH4(6,80000)

!      output from canopy model
      real evap_yr,transp_yr
      real,dimension(10):: thksl,wupl,evapl,wcl,FRLEN
      real runoff,runoff_d,runoff_yr,rain,rain_d,rain_yr
      real wsc(10),ws1,ws2,dws,net_dws
      real Esoil,Hcrop,ecstot,Anet,DEPH2O,Acanop,G
      real Hcanop,Hcanop_d
      real Raplant,Glmax,Gsmax,Rh_d
      real GLmx,Gsmx,GRmx
!      output for ORNL model comparison
      real CO2h,PARh,ATh,STh,VPDh,SWh
      real RECOh
      real ETh,Th,Eh,INTh,ROh,DRAINh,LEh,SHh
      real LWH,Rgrowth_d,abvLitter,blLitter
!     daily output
      real PAR_d,AT_d,ST_d,VPD_d
      real SW_d,NEP_d,NEE_d,RECO_d
      real Ra_d,RLEAV_d,RWOOD_d,RROOT_d,RHET_d,RSOIL_d,ET_d,T_d
      real E_d,INT_d,RO_d,DRAIN_d,LE_d,SH_d,CL_d,CW_d,CFR_d,TNC_d
      real CSOIL_d,GL_d,GW_d,GR_d,LFALL_d,LMA_d,NCAN_d,NWOOD_d
      real GL_yr,GR_yr,GW_yr
      real NFR_d,NSOIL_d,NUP_d,NMIN_d,NVOL_d,NLEACH_d
      real N_LG_d,N_WG_d,N_RG_d
      real N_LF_d,N_WF_d,N_RF_d
      real WFALL_D,RFALL_D
      real Simu_lit,zwt,phi
      ! for soil temp
      real sftmp,Tsnow,Twater,Tice,ice_tw,water_tw 
      real,dimension(10):: Tsoill,ice,liq_water
      real,dimension(11):: testout
      real soilt_d_simu(11),soilt_d_obs(7),watertable_d_obs,ice_d_simu(10)
      integer obs_counter(7) 
      real zwt_d,snow_depth_e,snow_dsim,melt,dcount,dcount_soil
      character(len=80) outfile

!      NEE observation
      real NEE_annual,Cumol2gram
      real NEE_annual_array(30)
      integer year_array(30),year_obs
!     for loops
      integer jrain,W_flag(7)
      integer onset !flag of phenological stage
      integer year,yr,days,i,j,m,n,yrs_eq,hoy,iyr,daily
      integer k1
      integer lines_NEE,yr_NEE
      integer istat1,istat2,istat3,istat4
      integer dtimes,yr_length
      integer num_scen,isite
      integer idoy,ihour,ileaf,first_year,dlayer
      
      !*****for methane subroutine      MS
      integer, parameter :: nlayers=10
      real CH4(nlayers)                                                      !MS
      real ProCH4(nlayers),Pro_sum,Pro_sum_d,Pro_sum_yr
      real OxiCH4(nlayers),Oxi_sum,Oxi_sum_d,Oxi_sum_yr
      real simuCH4,simuCH4_d,simuCH4_yr
      real Fdifu(nlayers),Fdifu1_d,Fdifu1_yr
      real Ebu_sum,Ebu_sum_d,Ebu_sum_yr
      real Pla_sum,Pla_sum_d,Pla_sum_yr
      real S_omega
      real r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi
      logical do_snow,do_soilphy
!*******************************************

!     Default C/N ratios of Oak Ridge FACE
!      data CN0 /45.,350.,60.,40.,300.,10.,15.,8./
!     Default C/N ratios of Duke FACE
      data CN0 /50.,350.,60.,40.,300.,10.,20.,12./
!     thickness of every soil layer
      !thksl=(/5.,15.,20.,40.,100.,20.,b_bound,40.,40.,40./)
      data thksl /20.,20.,20.,20.,20.,20.,20.,20.,20.,20./
      
!     ratio of roots in every layer, Oak Ridge FACE
      data FRLEN /0.1,0.25,0.25,0.2,0.1,0.05,0.025,0.015,0.005,0.005/
     
      sftmp =-0.
      Tsnow = -20.
      Twater=0.0
      Tice =0.0
      
      G=20.5
      Esoil=0.5*G
      snow_dsim =0.575
      dcount=50.
      dcount_soil=50.
    
      ice_tw =0.0
      Tsoill=(/ -0.09, 0.73, 1.95, 3.89, 5.88, 6.27, 3.89,    &
        &   5.99, 6.27,7./)
        !-0.24
      ice=(/0., 0., 0., 0., 0.0, 0.0, 0.0, 0.0,    &
        &   0.0, 0.0/) 
        
      liq_water=(/0., 0., 0., 0., 0., 0., 0.,0.,0.,0./) 
      zwt=0.0
      water_tw=zwt*0.001
          
!     Nitrogen input
!      N_deposit=0.000144634702 !(gN/h/m2, 1.2+0.067 gN/yr/m2,Oak ridge)
!     0.7 gN/yr/m2, 13.4 kg N ha-1 yr-1, 2000, Dentener et al. 2006, GBC, Duke FACE
      N_deposit=2.34/8760. !(gN/h/m2, )

!      N_fert=0. ! (20.0 gN m-2 yr-1, in spring, from 2004, Oak Ridge)
      N_fert=0. !5.6 ! (11.2 gN m-2 yr-1, in spring, Duke Forest FACE)

!         the unit of residence time is transformed from yearly to hourly
          tauC=(/tau_L,tau_W,tau_R,tau_F,tau_C,&
     &           tau_Micr,tau_Slow,tau_Pass/)*8760.
     
          SLA=SLAx/10000.         ! Convert unit from cm2/g to m2/g
!         growth rates of plant
          GLmax=GLmx/8760.
          GRmax=GRmx/8760.
          Gsmax=GSmx/8760.
!         end of setting parameters
 
          input_data=forcing_data
!         end of reading forcing data
!     ===============================================================
!         cycle  ! skip the following blocks, read input data only.
!         ===================================================
!         Initialize parameters and initial state:

          WILTPT=wsmin/100.0
          FILDCP=wsmax/100.0
!         define soil for export variables for satisfying usage of canopy submodel first time
!          wscontent=WILTPT
          infilt=0.
!         gddonset=320.0
          phi=wsmax/100.
          
          stor_use=Storage/times_storage_use
          accumulation=0.0
          SNvcmax=1.0

          LAI=LAIMIN
          !LAI=4.
          bmleaf=QC(1)/0.48
          bmstem=QC(2)/0.48
          bmroot=QC(3)/0.48
          bmplant=bmstem+bmroot+bmleaf

!         initial values of Nitrogen pools and C/N ratio
          alphaN=0.0    ! the transfer of N before littering

          NSN=6.0
          QNminer= 1.2
          N_deficit=0
          CN=CN0
          QN=QC/CN0
          QNplant  =QN(1) + QN(2) + QN(3)
!=============================================================
          m=1
          n=1
          k1=1
          iyr=0
          idays=365
          daily=0
          first_year=2011
          do i=1,nlayers                                                          !MS
              CH4(i)=0
          enddo
          
          !yr_length=1
          do yr=1,yrs_eq+yr_length  ! how many years
           !do yr=1,1    ! Yuan test
              if(yr.gt.3)then
!                  write(*,*)'One year done'
              endif
!            using ambient data to run equilibiurm, elevated only for the last cycle
             iyr=iyr+1
             if(iyr>yr_length)iyr=1

!!          leap year
           if(MOD(first_year+iyr-1,4).eq.0)then
                 idays=366
           else
                 idays=365
           endif
                     
             GDD5=0.0
             onset=0
             phenoset=0
             diff_yr=0.0
             gpp_yr=0.0
             R_Ntr_yr=0.
             NPP_yr=0.0
             Rh_yr =0.0
             Rh4_yr=0.0
             Rh5_yr=0.0
             Rh6_yr=0.0
             Rh7_yr=0.0
             Rh8_yr=0.0
             Ra_yr =0.0
             GL_yr=0.0
             GW_yr=0.0
             GR_yr=0.0
             Pool1=0.0
             Pool2=0.0
             Pool3=0.0
             Pool4=0.0
             Pool5=0.0
             Pool6=0.0
             Pool7=0.0
             Pool8=0.0
             out1_yr=0.0
             out2_yr=0.0
             out3_yr=0.0
             out4_yr=0.0
             out5_yr=0.0
             out6_yr=0.0
             out7_yr=0.0
             out8_yr=0.0             
             NEE_yr=0.0
!            water fluxes
             rain_yr=0.0
             transp_yr=0.0
             evap_yr=0.0
             runoff_yr=0.0
             Simu_lit=0.
             
!            Nitrogen fluxes
             N_up_yr=0
             N_fix_yr=0.
             N_dep_yr=0.
             N_leach_yr=0.
             N_vol_yr=0.
             
        !============================== test variable
             fwsoil_yr=0.
             omega_yr=0.
             topfws_yr=0.
             
             hoy=0
!!     end of leap year
             do days=1,idays !the days of a year
             ! do days=1,1   ! Yuan test
                 if(days.gt.72)then
!                     write(*,*)'pause'
                 endif
!             Nitrogen fertilization since 2004 in Oak Ridge
!              if(yr>yrs_eq+5.and.days==135)then
!                  QNminer=QNminer+N_fert     !(20 gN/yr/m2,N fertiliztion in Spring)
!              endif

!             Nitrogen fertilization since 1999 in Duke
              if(yr>yrs_eq+1.and.(days==75.OR.days==105))then
                  QNminer=QNminer+N_fert     !(5.6 gN/yr/m2,N fertiliztion in March and Apr)
              endif

              StemSap=AMIN1(Stemmax,SapS*bmStem)   ! Stemmax and SapS were input from parameter file, what are they? Unit? Maximum stem biomass? -JJJJJJJJJJJJJJJJJJJJJJ 
              RootSap=AMIN1(Rootmax,SapR*bmRoot)
              NSCmin=5. 
              NSCmax=0.05*(StemSap+RootSap+QC(1))
              
              soilt_d_simu=(/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./) 
              ice_d_simu=(/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./) 
              
              soilt_d_obs=(/0.,0.,0.,0.,0.,0.,0./) 
              zwt_d=0.0
              obs_counter = (/0,0,0,0,0,0,0/) 
              
              if(Ta.gt.5.0)GDD5=GDD5+Ta
              
              if (do_snow) then 
                 call snow_d(rain_d,lat,days,ta,snow_dsim,fa,fsub,rho_snow,melt,dcount,decay_m)                            
                 snow_depth_e=snow_dsim
              endif 
              
!             THE FIRST PART:  coupled canopy and soil model
              diff_d = 0.0
              gpp_d   =0.0   ! daily
              gpp_ra  =0.0   ! daily
              NPP_d   =0.0   ! daily
              NEP_d=0.0
              NEE_d = 0.0
!             rain_d,transp_d,evap_d
              transp_d=0.0   ! daily
              Hcanop_d=0.0   ! daily
              evap_d  =0.0   ! daily
              ta=0.0         ! daily 
              Ts=0.0         ! daily
              rain_d=0.0     ! daily
              runoff_d=0.0    ! daily
              LE_d=0.0
              RaL=0.0
              RaS=0.0
              RaR=0.0
              Rauto=0.0
              Rh_d=0.0
              N_up_d=0.
              N_fix_d=0.
              N_dep_d=0.
              N_leach_d=0.
              N_vol_d=0.
              PAR_d=0.
              VPD_d=0.0
              RECO_d=0.0
              RLEAV_d=0.0 
              RWOOD_d=0.0
              RROOT_d=0.0
              GL_d   =0.0
              GW_d   =0.0
              GR_d   =0.0
              LFALL_d=0.0
              NUP_d=0.0
              NVOL_d=0.
              NLEACH_d=0.0
              NMIN_d=0.0
              N_LG_d=0.0
              N_WG_d=0.0
              N_RG_d=0.0
              N_LF_d=0.0
              N_WF_d=0.0
              N_RF_d=0.0
              WFALL_d=0.0
              RFALL_d=0.0
              
              !*********** This session added for initial value of methane daily variables
              simuCH4_d=0.0
              Pro_sum_d=0.0
              Oxi_sum_d=0.0
              Fdifu1_d=0.0
              Ebu_sum_d=0.0
              Pla_sum_d=0.0
!*****************************             
              dtimes=24 !how many times a day,24 means every hour
              do i=1,dtimes
                  if(i.gt.11)then
!                      write(*,*)'pause'
                  endif
                    
!                 input data
                  if(m > lines)then    ! Repeat forcing data for the whole time period
                !  if(yr.ge.1 .and. mod(days,idays).eq.0     &
                !  &     .and. mod(i,dtimes).eq.0)then   ! Repeat forcing data for specific time period
                     m=1		! m is row sequence in the file of input_data
                     n=1
                     hoy=0
                  endif
                  
                  year =year_seq(m)
                  doy  =doy_seq(m)
                  hour =hour_seq(m)+1
                  
                  !!       for Duke Forest
                  Tair=input_data(1,m)+0.0   ! Tair
                  Tsoil=input_data(2,m)    ! SLT
                  RH=input_data(3,m)
                  Dair=input_data(4,m)       !air water vapour defficit? Unit Pa
                  co2ca=380.0*1.0E-6 ! CO2 concentration,ppm-->1.0E-6
!                  if(isite==2.and.yr>yrs_eq)then
!                      co2ca=(input_data(5,m)+200.)*1.0E-6
!                  endif
                  day_mod=mod(m,24)
                  rain=input_data(5,m)    ! rain fal per hour
                  wind=ABS(input_data(6,m))     ! wind speed m s-1
                  PAR=input_data(7,m)             ! Unit ? umol/s/m-2
                  radsol=input_data(7,m)        ! unit ? PAR actually
                  
                 ! water_table_depth=water_table(m)
                 ! if (water_table_depth .gt. -9999) water_table_depth = water_table_depth!*0.1   !! change mm to cm
                  !water_tw=water_table_depth*0.01
!!       endof Duke Forest
                  if (do_snow) then
                      snow_depth=snow_depth_e
                  else
                      snow_depth=snow_in(m)
                  endif

                  if (snow_depth .lt. 0.0) snow_depth = 0.0   
                  snow_depth = snow_depth*100.   ! change from m to cm
                  
!                  Rnet=input_data(9,m)
!                 Ajust some unreasonable values
                  RH=AMAX1(0.01,AMIN1(99.99,RH))
                  eairP = esat(Tair)*RH/100.             ! Added for SPRUCE, due to lack of VPD data
                  Dair=esat(Tair)-eairP
                  radsol=AMAX1(radsol,0.01)
                  hoy=hoy+1
                m=m+1											

!                  if(radsol.gt.10.0) then
!                      G=-25.0
!                  else
!                      G=20.5											! What's G? -JJJJJJJJJJJJJJJJJJJJJJ
!                  endif
!                  Esoil=0.05*radsol
!                  if(radsol.LE.10.0) Esoil=0.5*G						! If solar radiation less than 10, Esoil = 10.25 -JJJJJJJJJJJJJJJJJJJJJ
                  
                  Hcrop=0.1  ! never used in routine
                  Ecstot=0.1 ! never used in routine
                  Anet=0.1 ! never used in routine
                  DepH2O=0.2
!                 for daily mean conditions
                  VPD_d=VPD_d+Dair/24./1000.                ! Water vapour pressure or depecit ?? -JJJJJJJJJJJJJJJJJ, Change unit from Pa to KPa
                  PAR_D=PAR_D+radsol/dtimes                 ! umol photons m-2 s-1   
                  ta= ta + tair/24.0             ! sum of a day, for calculating daily mean temperature
                  Ts=Ts+Tsoil/24.0
                  rain_d=rain_d+rain
!                 calculating scaling factor of NSC
                  if(NSC.le.NSCmin)fnsc=0.0
                  if(NSC.ge.NSCmax)fnsc=1.0
                  if((NSC.lt.NSCmax).and.(NSC.gt.NSCmin))then 
                     fnsc=(NSC-NSCmin)/(NSCmax-NSCmin)
                  endif
!                 update vcmx0 and eJmx0 according to C/N of leaves
                  Vcmx0 = Vcmax0*SNvcmax*1.0e-6
!                  eJmx0 = 2.7*Vcmx0  ! original
                  eJmx0 = 1.67*Vcmx0 ! Weng 02/21/2011 Medlyn et al. 2002
                  if(yr.eq.2 .and. days.eq.124 .and. i.eq.1)then
!                      write(*,*)'pause'
                  endif
                  
                  call canopy(gpp,evap,transp,Acanop,Hcanop,Rsoilabs,  & ! outputs
     &              fwsoil,topfws, &                   ! from soil model
     &              LAI,Sps,&
     &              doy,hour,radsol,tair,dair,eairP,    &        ! from climate data file,including 
     &              wind,rain,&
     &              Rnet,G,Esoil,Hcrop,Ecstot,Anet,&
     &              Tsoil,DepH2O,&
     &              wsmax,wsmin,  &                              !constants specific to soil and plant
     &              lat,co2ca,a1,Ds0,Vcmx0,extkU,xfang,alpha,&
     &              stom_n,pi,tauL,rhoL,rhoS,emleaf,emsoil,&
     &              Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,&
     &              H2OMw,chi,Dheat,wleaf,gsw0,eJmx0,theta,&
     &              conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,&
     &              Edjm,Entrpy,gam0,gam1,gam2,wcl,gddonset,sftmp,Tsoill,testout,ice,&
     &              water_table_depth,snow_depth,Tsnow,Twater,Tice,water_tw,ice_tw,diff_s,&
     &              diff_snow,albedo_snow,resht,thd_snow_depth,THKSL,zwt,liq_water,shcap_snow,&
     &              condu_snow,condu_b,depth_ex,dcount_soil,do_soilphy)

                 call soilwater(wsmax,wsmin,rdepth,FRLEN,THKSL,    &   !constants specific to soil/plant
     &                rain,transp,evap,wcl,runoff,infilt,   &   !inputs
     &                fwsoil,topfws,omega,wsc,zwt,phi,liq_water,infilt_rate,melt,ta,day_mod, &
     &                do_soilphy,snow_depth,ice)               !outputs

                  ET=evap+transp
                  rain_yr=rain_yr+rain
                  transp_yr=transp_yr+transp
                  evap_yr=evap_yr+evap
                  runoff_yr=runoff_yr+runoff

                  call respiration(LAIMIN,GPP,Tair,Tsoil,DepH2O,&
     &                       Q10,Rl0,Rs0,Rr0,SNRauto,&
     &                       LAI,SLA,bmstem,bmroot,bmleaf,&
     &                       StemSap,RootSap,NSC,fnsc,&
     &                       RmLeaf,RmStem,RmRoot,Rmain)

!                 THE Third Part: update LAI
                  call plantgrowth(Tair,omega,GLmax,GRmax,GSmax,&
     &                    LAI,LAIMAX,LAIMIN,SLA,TauC(1),         &    !Tau_L,
     &                    bmleaf,bmroot,bmstem,bmplant,&
     &                    Rootmax,Stemmax,SapS,SapR,&
     &                    StemSap,RootSap,Storage,GDD5,&
     &                    stor_use,onset,accumulation,gddonset,&
     &                    Sps,NSC,fnsc,NSCmin,NSCmax,&
     &                    NSN,CN,CN0,SNgrowth,N_deficit,&
     &                    store,add,L_fall,ht,&
     &                    NPP,alpha_L,alpha_W,alpha_R,&
     &                    RgLeaf,RgStem,RgRoot,Rgrowth)

!                 THE Fourth PART: simulating C influx allocation in pools
                  call TCS_CN(Tair,Tsoil,omega,runoff,&
     &               NPP,alpha_L,alpha_W,alpha_R,L_fall,&
     &               tauC,QC,OutC,Rh_pools,Rnitrogen,NSC,&
     &               CNmin,CNmax,NSNmax,NSNmin,alphaN,   &         ! nitrogen
     &               NSN,N_uptake,N_miner,QN,QNminer,&
     &               CN,CN0,fnsc,rdepth,N_deficit,&
     &               N_leaf,N_wood,N_root,N_LF,N_WF,N_RF,&
     &               N_deposit,N_fixation,N_leach,N_vol,&
     &               SNvcmax,SNgrowth,SNRauto,SNrs,Q10,tsoill,testout,do_soilphy)
     
                call methane(Rh_pools,Tsoil,zwt,wsc,     &      !update single value in a hourly loop when MEMCMC=0
                &           phi,LAIMIN,LAIMAX,           &
                &           ProCH4,Pro_sum,OxiCH4,Oxi_sum,Fdifu,Ebu_sum,Pla_sum,simuCH4,CH4,   &
                &           r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi,  &
                &           testout,do_soilphy)       !update single value of Rh_pools,Tsoil,zwt,wsc 
        !          

!                 update NSC
                  Rauto  =Rmain+Rgrowth+Rnitrogen
                  NSC    =NSC+GPP-Rauto-(NPP-add)-store
                  Difference = GPP-Rauto-NPP
                  if(NSC<0)then
                      bmstem=bmstem+NSC/0.48
                      NPP=NPP+NSC
                      NSN=NSN-NSC/CN(2)
                      NSC=0.
                  endif
                  GL_d   =GL_d+NPP*alpha_L
                  GW_d   =GW_d+NPP*alpha_W
                  GR_d   =GR_d+NPP*alpha_R
                  LFALL_d=LFALL_d+L_fall
!                 update
                  RaLeaf = RgLeaf + RmLeaf
		  RaStem = RgStem + RmStem
		  RaRoot = RgRoot + RmRoot + Rnitrogen
                  WFALL_d=WFALL_d+OutC(2) !_wood
                  RFALL_d=RFALL_d+OutC(3) !_root
                  N_LG_d=N_LG_d+N_leaf
                  N_WG_d=N_WG_d+N_wood
                  N_RG_d=N_RG_d+N_root
                  N_LF_d=N_LF_d+N_LF
                  N_WF_d=N_WF_d+N_WF
                  N_RF_d=N_RF_d+N_RF

                  N_up_d=N_up_d+N_uptake
                  N_fix_d=N_fix_d+N_fixation
                  N_dep_d=N_dep_d+N_deposit
                  N_leach_d=N_leach_d+N_leach
                  N_vol_d=N_vol_d+N_vol

                  N_up_yr=N_up_yr+N_uptake
                  N_fix_yr=N_fix_yr+N_fixation
                  N_dep_yr=N_dep_yr+N_deposit
                  N_leach_yr=N_leach_yr+N_leach
                  N_vol_yr=N_vol_yr+N_vol

                  R_Ntr_yr=R_Ntr_yr + Rnitrogen
                  do dlayer=1,10
                      ice_d_simu(dlayer)=ice_d_simu(dlayer)+ice(dlayer) 
                  enddo       
                   do dlayer=1,11
                      soilt_d_simu(dlayer)=soilt_d_simu(dlayer)+testout(dlayer)  
                      ! first = surface soil temperature 2:11=1:10 layer soil temperatures 
                  enddo                    
                  
                  zwt_d=zwt_d+zwt                  
           ! ==================== test variables
                  topfws_yr = topfws_yr+topfws/8760.
                  omega_yr=omega_yr+omega/8760.
                  fwsoil_yr=fwsoil_yr+fwsoil/8760.

!                 Rhetero=Rh_f + Rh_c + Rh_Micr + Rh_Slow + Rh_Pass
                  Rhetero= Rh_pools(1)+Rh_pools(2)+Rh_pools(3) &
     &                    +Rh_pools(4)+Rh_pools(5)
                  Rsoil  =Rhetero+RmRoot+RgRoot+Rnitrogen
                  NEE=Rauto+Rhetero - GPP
                  Q_soil=QC(6) + QC(7) + QC(8)

                  bmleaf=QC(1)/0.48
                  bmstem=QC(2)/0.48
                  bmroot=QC(3)/0.48
                  bmplant=bmleaf+bmroot+bmstem
                  LAI=bmleaf*SLA
                  NMIN_d = NMIN_d+N_miner
!                 output hourly
                  Recoh=Rhetero+Rauto
                  ETh =ET !*1000.
                  Th  =transp !*1000.
                  Eh  =evap !*1000.
                  INTh=-9999
                  VPDh=Dair/1000.
                  ROh =runoff !*1000.
                  DRAINh=-9999
                  LEh =ETh*((2.501-0.00236*Tair)*1000.0)/3600.
                  SHh =-9999
                  LWh =-9999
                  NEP=-NEE
!                 sums of a day
                  diff_d=diff_d+difference
                  gpp_d=gpp_d + GPP
                  gpp_ra=gpp_ra+Rauto
                  NPP_d   =NPP_d+NPP
                  NEP_d=NEP_d+NEP
                  NEE_d=NEE_d+NEE
                  RECO_d=RECO_d+Recoh
                  Rh_d=  Rh_d + Rhetero
                  Ra_d=Reco_d-Rh_d
                  RLEAV_d=RLEAV_d+RmLeaf+RgLeaf
                  RWOOD_d=RWOOD_d+RmStem+RgStem
                  RROOT_d=RROOT_d+RmRoot+RgRoot+Rnitrogen
                  Rsoil_d=Rh_d+RROOT_d
                  NUP_d=NUP_d+N_uptake
                  NVOL_d=NVOL_d+N_vol
                  NLEACH_d=NLEACH_d+N_leach
                  transp_d=transp_d + transp*(24./dtimes)
                  evap_d=evap_d + evap*(24./dtimes)
                  ET_d=transp_d + evap_d
                  LE_d=LE_d+LEh/24.
                  Hcanop_d=Hcanop_d+Hcanop/(24./dtimes)
                  runoff_d=runoff_d+runoff                 
                  ! added for MEMCMC also for generation of daily methane emission                  
                  simuCH4_d=simuCH4_d+simuCH4
                  Pro_sum_d=Pro_sum_d+Pro_sum
                  Oxi_sum_d=Oxi_sum_d+Oxi_sum
                  Fdifu1_d=Fdifu1_d+Fdifu(1)
                  Ebu_sum_d=Ebu_sum_d+Ebu_sum
                  Pla_sum_d=Pla_sum_d+Pla_sum                                              
                  
!                 sum of the whole year
                  diff_yr = diff_yr+difference
                  gpp_yr=gpp_yr+gpp
                  NPP_yr=NPP_yr+NPP
                  Rh_yr =Rh_yr +Rhetero
                  Ra_yr=Ra_yr+Rauto
                  Rh4_yr=Rh4_yr+Rh_pools(1)
                  Rh5_yr=Rh5_yr+Rh_pools(2)
                  Rh6_yr=Rh6_yr+Rh_pools(3)
                  Rh7_yr=Rh7_yr+Rh_pools(4)
                  Rh8_yr=Rh8_yr+Rh_pools(5)
                  Pool1 = Pool1+QC(1)/8760.
                  Pool2 = Pool2+QC(2)/8760.
                  Pool3 = Pool3+QC(3)/8760.
                  Pool4 = Pool4+QC(4)/8760.
                  Pool5 = Pool5+QC(5)/8760.
                  Pool6 = Pool6+QC(6)/8760.
                  Pool7 = Pool7+QC(7)/8760.
                  Pool8 = Pool8+QC(8)/8760.
                  out1_yr=out1_yr+OutC(1)
                  out2_yr=out2_yr+OutC(2)
                  out3_yr=out3_yr+OutC(3)
                  out4_yr=out4_yr+OutC(4)
                  out5_yr=out5_yr+OutC(5)
                  out6_yr=out6_yr+OutC(6)
                  out7_yr=out7_yr+OutC(7)
                  out8_yr=out8_yr+OutC(8)
                  NEE_yr=NEE_yr+NEE
                  GL_yr=GL_yr+NPP*alpha_L
                  GW_yr=GW_yr+NPP*alpha_W
                  GR_yr=GR_yr+NPP*alpha_R
!                 numbering         
                  n=n+1
                  if((yr+first_year-1).eq.obs_soilwater(1,k1) .and.    &
                  &     days .eq. obs_soilwater(2,k1) .and.      &
                  &     (i-1).eq. obs_soilwater(3,k1))then
                  Simu_soilwater(1:10,k1)=wcl(1:10)
                  Simu_soiltemp(1:11,k1)=testout
                  Simu_watertable(1,k1)=zwt
                  
                  k1=k1+1
                  endif
                 ! write(*,*)yr,days,i,gpp,npp
                  if(isnan(gpp))then
                      write(*,*)'gpp is nan'
                  endif

              enddo              ! end of dtimes
              if((GDD5.gt.gddonset) .and. phenoset.eq.0) then
                pheno=days
                phenoset=1
              endif

!             +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++             
!             output results of canopy and soil models, daily
!             daily output
              INT_d=-9999
              DRAIN_d=-9999
              SH_d=-9999
              CSOIL_d=QC(6)+QC(7)+QC(8)
              LMA_d=QC(1)/0.48
              NSOIL_d=QN(6)+QN(7)+QN(8)+QNminer

              Rgrowth_d=-9999
              abvLitter=QC(4)  !-9999
              blLitter=QC(5) !-9999
!            write(*,*)yr,days,LAI,gpp_d,npp_d

            !if(yr.gt.yrs_eq)then  
                daily=daily+1
  
            Simu_dailyflux(1,daily)=GPP_d ! Leaf
            Simu_dailyflux(2,daily)=NEE_d	! Wood
            Simu_dailyflux(3,daily)=Reco_d	! Coarse roots
            Simu_dailyflux(4,daily)=QC(1)!*1.5
            Simu_dailyflux(5,daily)=GL_yr!*1.5
            Simu_dailyflux(6,daily)=QC(2)!*0.48
            Simu_dailyflux(7,daily)=GW_yr!*0.48
            Simu_dailyflux(8,daily)=QC(3)
            Simu_dailyflux(9,daily)=GR_yr
            Simu_dailyflux(10,daily)=(QC(6)+QC(7)+QC(8))!*13.8   ! Soil
            Simu_dailyflux(11,daily)=pheno   ! Soil
            Simu_dailyflux(12,daily)=LAI !QC(1)/(QC(1)+QC(4)) 
            
            Simu_dailyCH4(1,daily)=simuCH4_d
            Simu_dailyCH4(2,daily)=Pro_sum_d
            Simu_dailyCH4(3,daily)=Oxi_sum_d
            Simu_dailyCH4(4,daily)=Fdifu1_d
            Simu_dailyCH4(5,daily)=Ebu_sum_d
            Simu_dailyCH4(6,daily)=Pla_sum_d
            
            
            Simu_dailysoilt(1:11,daily)=soilt_d_simu(1:11)/24.
            
            Simu_dailyice(1:10,daily)=ice_d_simu(1:10)/24.
                                                        
            Simu_dailywatertable(1,daily)=zwt_d/24.
            Simu_snowdepth(1,daily)=snow_dsim  
            !endif            
        enddo                         ! end of idays
            storage=accumulation
            stor_use=Storage/times_storage_use
            if(yr.eq.yrs_eq+yr_length .and. MCMC.eq.1)then
            !write(*,*)yr,LAI,gpp_yr,NPP_yr,pheno   ! Yuanyuan
            write(61,601)year,LAI,gpp_yr,NPP_yr,real(pheno)
            endif
            if(MCMC.eq.0) then
                write(*,*)'year, LAI, GPP_YR,NPP_YR,PHENO:',year,LAI,gpp_yr,NPP_yr,pheno
                write(61,601)year,LAI,gpp_yr,NPP_yr,Ra_yr,Rh_yr, &
                &   Rh4_yr,Rh5_yr,Rh6_yr,Rh7_yr,Rh8_yr,GL_yr,    &
                &   GW_yr,GR_yr,Pool1,Pool2,Pool3,Pool4,Pool5,   &
                &   Pool6,Pool7,Pool8,out1_yr,out2_yr,out3_yr,   &
                &   out4_yr,out5_yr,out6_yr,out7_yr,out8_yr
            endif            
601         format(i7,",",29(f15.4,","))
            accumulation=0.0
            onset=0
         enddo            !end of simulations multiple years
         
         num=0
         if(MCMC.eq.0)then
             
         write(outfile,"(A20,I5.5,A4)") "output/Simu_dailyflux",num,".txt"
         outfile=trim(outfile)
         open(62,file=outfile)
         do i=1,daily
          write(62,602)i,(Simu_dailyflux(j,i),j=1,12)
         enddo
         close(62)          
         
         !******* added for MEMCMC output files, same format as above ********
         write(outfile,"(A20,I5.5,A4)") "output/Simu_dailyCH4",num,".txt"
         outfile=trim(outfile)
         open(664,file=outfile)
             do i=1,daily
                write (664,604) i, (Simu_dailyCH4(j,i),j=1,6)
             enddo         
         close (664)

!         open(63,file='Simu_soilwater001.txt')
!         do i=1,k1-1
!         write(63,603)i,(Simu_soilwater(j,i),j=1,10)
!         enddo   
         
         write(outfile,"(A26,I5.5,A4)") "output/Simu_soiltemp_daily",num,".txt"
         outfile=trim(outfile)
         open(63,file=outfile)
         do i=1,daily
             write(63,613)i,(Simu_dailysoilt(j,i),j=1,11)
         enddo
         close(63)
         
         write(outfile,"(A28,I5.5,A4)") "output/Simu_watertable_daily",num,".txt"
         outfile=trim(outfile)
         open(6333,file=outfile)
         do i=1,daily
             write(6333,604)i,(Simu_dailywatertable(j,i),j=1,1)
         enddo
         close(6333)
                 
!         open(64,file='test_yuan.txt')
!         do i=1,k1-1
!         write(64,613)i,(Simu_soiltemp(j,i),j=1,11)
!         enddo
                 
         endif
602     format((i7),",",12(f15.4,","))
603     format((i7),",",10(f15.4,","))
604     format((i7),",",10(f15.4,","))
613     format((i7),",",11(f15.4,","))
999      continue
    return
    end
      

!     ****************************************************************************
      subroutine consts(pi,tauL,rhoL,rhoS,emleaf,emsoil,                &
     &   Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat,    &
     &   wleaf,gsw0,Vcmx0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci,       &
     &   Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2) 
     
      real tauL(3), rhoL(3), rhoS(3)
      pi = 3.1415926
!     physical constants
      tauL(1)=0.1                  ! leaf transmittance for vis
      rhoL(1)=0.1                  ! leaf reflectance for vis
      rhoS(1)=0.1                  ! soil reflectance for vis
      tauL(2)=0.425                ! for NIR
      rhoL(2)=0.425                ! for NIR
      rhoS(2)=0.3                  ! for NIR - later function of soil water content
      tauL(3)=0.00                 ! for thermal
      rhoL(3)=0.00                 ! for thermal
      rhoS(3)=0.00                 ! for thermal
      emleaf=0.96
      emsoil=0.94
      Rconst=8.314                 ! universal gas constant (J/mol)
      sigma=5.67e-8                ! Steffan Boltzman constant (W/m2/K4)
      cpair=1010.                  ! heat capapcity of air (J/kg/K)
      Patm=101325. !1.e5           ! atmospheric pressure  (Pa)
      Trefk=293.2                  !reference temp K for Kc, Ko, Rd
      H2OLv0=2.501e6               !latent heat H2O (J/kg)
      AirMa=29.e-3                 !mol mass air (kg/mol)
      H2OMw=18.e-3                 !mol mass H2O (kg/mol)
      chi=0.93                     !gbH/gbw
      Dheat=21.5e-6                !molecular diffusivity for heat
!     plant parameters
      gsw0 = 1.0e-2                !g0 for H2O in BWB model
      eJmx0 = Vcmx0*2.7            !@20C Leuning 1996 from Wullschleger (1993)
      theta = 0.9
      wleaf=0.01                   !leaf width (m)

!     thermodynamic parameters for Kc and Ko (Leuning 1990)
      conKc0 = 302.e-6                !mol mol^-1
      conKo0 = 256.e-3                !mol mol^-1
      Ekc = 59430.                    !J mol^-1
      Eko = 36000.                    !J mol^-1
!     Erd = 53000.                    !J mol^-1
      o2ci= 210.e-3                   !mol mol^-1

!     thermodynamic parameters for Vcmax & Jmax (Eq 9, Harley et al, 1992; #1392)
      Eavm = 116300.               !J/mol  (activation energy)
      Edvm = 202900.               !J/mol  (deactivation energy)
      Eajm = 79500.                !J/mol  (activation energy) 
      Edjm = 201000.               !J/mol  (deactivation energy)
      Entrpy = 650.                !J/mol/K (entropy term, for Jmax & Vcmax)

!     parameters for temperature dependence of gamma* (revised from von Caemmerer et al 1993)
      gam0 = 28.0e-6               !mol mol^-1 @ 20C = 36.9 @ 25C
      gam1 = .0509
      gam2 = .0010
      return
      end
!****************************************************************************
!      a sub-model for calculating C flux and H2O flux of a canopy
!      adapted from a two-leaf canopy model developed by Wang Yingping
         
        subroutine canopy(gpp,evap,transp,Acanop,Hcanop,Rsoilabs, &  ! outputs
     &                    fwsoil,topfws,           &! from soil model
     &                    LAI,Sps,&
     &               doy,hour,radsol,tair,Dair,eairP,&! from climate data file,including 
     &               wind,rain,&
     &               Rnet,G,Esoil,Hcrop,Ecstot,Anet,&
     &               Tsoil,DepH2O,&
     &               wsmax,wsmin,&  !constants specific to soil and plant
     &               lat,co2ca,a1,Ds0,Vcmx0,extkU,xfang,alpha,&
     &               stom_n,pi,tauL,rhoL,rhoS,emleaf,emsoil,&
     &               Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,&
     &               H2OMw,chi,Dheat,wleaf,gsw0,eJmx0,theta,&
     &               conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,&
     &               Edjm,Entrpy,gam0,gam1,gam2,wcl,gddonset,sftmp,Tsoill,testout,ice,&
     &               water_table_depth,snow_depth,Tsnow,Twater,Tice,water_tw,ice_tw,diff_s,&
     &               diff_snow,albedo_snow,resht,thd_snow_depth,thksl,zwt,liq_water,shcap_snow,&
     &               condu_snow,condu_b,depth_ex,dcount_soil,do_soilphy)

      real lat,doy
      real gpp,evap,transp,LAI,Rsoilabs
      real tauL(3),rhoL(3),rhoS(3),reffbm(3),reffdf(3)
      real extkbm(3),extkdm(3)
      real Radabv(2),Qcan(3,2)
      real gddonset
!     extra variables used to run the model for the wagga data
      real topfws,Tsnow,Twater,Tice,water_tw,ice_tw        ! from siol subroutine      
      integer idoy,ihour,ileaf
      integer jrain,i,j,k
      real zwt
      logical do_soilphy

!     additional arrays to allow output of info for each layer
      real RnStL(5),QcanL(5),RcanL(5),AcanL(5),EcanL(5),HcanL(5)
      real GbwcL(5),GswcL(5),hG(5),hIL(5)
      real Gaussx(5),Gaussw(5),Gaussw_cum(5)
      real wcl(10),Tsoill(10),testout(11),ice(10),thksl(10),liq_water(10)
      character*80 commts
!     Normalised Gaussian points and weights (Goudriaan & van Laar, 1993, P98)
!     5-point
      data Gaussx/0.0469101,0.2307534,0.5,0.7692465,0.9530899/
      data Gaussw/0.1184635,0.2393144,0.2844444,0.2393144,0.1184635/
      data Gaussw_cum/0.11846,0.35777,0.64222,0.88153,1.0/
    
!     calculate beam fraction in incoming solar radiation
      call  yrday(doy,hour,lat,radsol,fbeam)
      idoy=int(doy)
      hours=idoy*1.0+hour/24.0
      coszen=sinbet(doy,lat,pi,hour)             !cos zenith angle of sun

!     set windspeed to the minimum speed to avoid zero Gb
      if(wind.lt.0.01) wind=0.01
!     calculate soil albedo for NIR as a function of soil water (Garratt pp292)
      if(topfws.gt.0.5) then
            rhoS(2)=0.18
      else
            rhoS(2)=0.52-0.68*topfws
      endif
!        assign plant biomass and leaf area index at time t
!        assume leaf biomass = root biomass
      FLAIT =LAI 
      eairP=esat(Tair)-Dair                !air water vapour pressure
      radabv(1)=0.5*radsol                 !(1) - solar radn
      radabv(2)=0.5*radsol                 !(2) - NIR
!     call multilayer model of Leuning - uses Gaussian integration but radiation scheme
!     is that of Goudriaan
      call xlayers(Sps,Tair,Dair,radabv,fbeam,eairP,&
     &           wind,co2ca,fwsoil,wcl,LAI,coszen,idoy,hours,&
     &           tauL,rhoL,rhoS,xfang,extkd,extkU,wleaf,&
     &           Rconst,sigma,emleaf,emsoil,theta,a1,Ds0,&
     &           cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,&
     &           gsw0,alpha,stom_n,wsmax,wsmin,&
     &           Vcmx0,eJmx0,conKc0,conKo0,Ekc,Eko,o2ci,&
     &           Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,&
     &           extKb,Rsoilabs,Acan1,Acan2,Ecan1,Ecan2,&
     &           RnStL,QcanL,RcanL,AcanL,EcanL,HcanL,GbwcL,GswcL,gddonset,testout,&
     &           Rsoilab1,Rsoilab2,QLleaf,QLair,raero,do_soilphy,G,Esoil,Hsoil)
     
     if (do_soilphy) then
       call Tsoil_simu(Rsoilab1,Rsoilab2,QLleaf,QLair,Tair,Dair,&
&         fbeam,FLAIT,sigma,emsoil,rhoS,Rconst,&
&         extkd,extkb,cpair,Patm,AirMa,H2OMw,&
&         H2OLv0,wcl,raero,wsmax,wsmin,wind,sftmp,Tsoill,testout,ht,ice,&
&         snow_depth,Tsnow,Twater,Tice,water_tw,ice_tw,diff_s,G,tsoil,&
&         diff_snow,albedo_snow,resht,thd_snow_depth,thksl,zwt,Esoil,Hsoil,liq_water,shcap_snow,&
&         condu_snow,condu_b,depth_ex,dcount_soil)
     endif
	 	
         Acanop=Acan1+Acan2
         Ecanop=Ecan1+Ecan2
         gpp=Acanop*3600.0*12.0                           ! every hour, g C m-2 h-1
         transp=AMAX1(Ecanop*3600.0/(1.0e6*(2.501-0.00236*Tair)),0.) ! mm H2O /hour
         evap=AMAX1(Esoil*3600.0/(1.0e6*(2.501-0.00236*Tair)),0.)
!        H2OLv0=2.501e6               !latent heat H2O (J/kg)

      return
      end      
!****************************************************************************
!   autotrophic respiration
    subroutine respiration(LAIMIN,GPP,Tair,Tsoil,DepH2O,&
     &                       Q10,Rl0,Rs0,Rr0,SNRauto,&
     &                       LAI,SLA,bmstem,bmroot,bmleaf,&
     &                       StemSap,RootSap,NSC,fnsc,&
     &                       RmLeaf,RmStem,RmRoot,Rmain)
!     calculate plant and soil respiration by the following equation:
!     RD=BM*Rd*Q10**((T-25)/10) (Sun et al. 2005. Acta Ecologica Sinica)
    implicit none
    real LAIMIN,LAI,GPP,SLA
    real Tair,Tsoil,DepH2O
    real bmstem,bmroot,bmleaf,StemSap,RootSap
    real NSC,fnsc
    real Q10
    real RmLeaf,RmStem,RmRoot,Rmain
    real Rl0,Rs0,Rr0,SNRauto
    real conv                  ! converter from "umol C /m2/s" to "gC/m2/hour"

    conv=3600.*12./1000000.    ! umol C /m2/s--> gC/m2/hour
    if(LAI.gt.LAIMIN) then
        RmLeaf=Rl0*SNRauto*bmleaf*0.48*SLA*0.1     &               
            &   *Q10**((Tair-10.)/10.)*fnsc*conv
        RmStem=Rs0*SNRauto*StemSap*0.001*Q10**((Tair-25.)/10.)*fnsc*conv
        RmRoot=Rr0*SNRauto*RootSap*0.001*Q10**((Tair-25.)/10.)*fnsc*conv
    else
        RmLeaf=0.3*GPP
        RmStem=0.3*GPP
        RmRoot=0.4*GPP
    endif
        Rmain=Rmleaf+Rmstem+Rmroot
    if(Rmain > 0.0015*NSC)then             ! If Total autotropic respiration greater than 0.15% of Nonstructure Carbon, rescale. 
        Rmleaf=Rmleaf/Rmain*0.0015*NSC
        Rmstem=Rmstem/Rmain*0.0015*NSC
        Rmroot=Rmstem/Rmain*0.0015*NSC
        Rmain=Rmleaf+Rmstem+Rmroot
    endif
    return
    end

!*******************************************************************
!     subroutine for soil moisture
    subroutine soilwater(wsmax,wsmin,rdepth,FRLEN,THKSL,    &   !constants specific to soil/plant
     &                rain,transp,evap,wcl,runoff,infilt,   &   !inputs
     &                fwsoil,topfws,omega,wsc,zwt,phi,liq_water,infilt_rate,melt,ta,day_mod, &
     &                do_soilphy,snow_depth,ice)               !outputs
!     All of inputs, the unit of water is 'mm', soil moisture or soil water content is a ratio
    implicit none
!   soil traits
    real wsmax,wsmin,wsmaxL(10),wsminL(10) !from input percent x%
    real FLDCAP,WILTPT ! ie. 0.xx
!   plant traits
    real rdepth,melt,ta,rain_new
    integer nfr,day_mod
!   climate conditions
    real rain,rain_t ! mm/hour
!   output from canopy model
    real evap,transp
!   output variables
    real fwsoil,topfws,omega
    real fw(10),ome(10)

    real thksl(10),depth(10),wsc(10),WUPL(10),EVAPL(10),SRDT(10),liq_water(10)
    real plantup(10)
    real Tsrdt
    real frlen(10) !fraction of root length in every layer
    real wcl(10) !volum ratio
 !      real fwcln(10) !  fraction of water in layers, like field capacity
    real DWCL(10),Tr_ratio(10)
    real wtadd,twtadd,infilt,runoff,tr_allo

    real exchangeL,supply,demand,omegaL(10)
    integer i,j,k
    real infilt_max
    logical do_soilphy
    real snow_depth
    real ice(10)
    
    !      characters annotation for water table module  -MS
    real vtot
    real phi
    real zmax
    real thetasmin
    real zthetasmin
    real az
    real zwt
    real zwt1
    real zwt2
    real zwt3
 !      water table characters annotation end here  -MS
    real infilt_rate

    infilt_max=25.
    WILTPT =wsmin/100.0
    FLDCAP =wsmax/100.0
    
    do i=1,10
        dwcl(i)=0.0
        evapl(i)=0.0
        WUPL(i)=0.0
        SRDT(i)=0.0
        DEPTH(i)=0.0
    enddo

!   Determine which layers are reached by the root system. 
!   Layer volume (cm3)
    DEPTH(1)=THKSL(1)
    DO i=2,10
        DEPTH(i)=DEPTH(i-1)+THKSL(i)
    enddo
    do i=1,10
        IF(rdepth.GT.DEPTH(i)) nfr=i+1
    enddo
    IF (nfr.GT.10) nfr=10
    
    
    if (do_soilphy) then 
       rain_new = rain
       if (ta .lt. -0.) rain_new =0.
       !if (day_mod .eq. 0) then
           rain_t =melt/24+rain_new
       !else
        !   rain_t=rain_new
       !endif
           infilt=infilt+rain_t
       if (ice(1) .gt. 0.0) then
           !infilt = 0.0
       endif
    else
        infilt=infilt+rain
    endif 
    
    
! *** water infiltration through layers
     !mm/hour
!   Loop over all soil layers.
    TWTADD=0
    IF(infilt.GT.0.0)THEN
!       Add water to this layer, pass extra water to the next.
        WTADD=AMIN1(INFILT,infilt_max,AMAX1((FLDCAP-wcl(1))*thksl(1)*10.0,0.0)) ! from cm to mm
!       change water content of this layer
        WCL(1)=(WCL(1)*(thksl(1)*10.0)+WTADD)/(thksl(1)*10.0)
!       FWCLN(I)=WCL(I)       !  /VOLUM(I)! update fwcln of this layer
        TWTADD=TWTADD+WTADD       !calculating total added water to soil layers (mm)
        INFILT=INFILT-WTADD !update infilt
        !write(*,*)"WTADD",WTADD
    ENDIF

!   runoff
    if (do_soilphy) then 
       if (wsc(1) .gt. 180.) then
!           if (snow_depth .gt. 0.0 .or. ice(1) .gt. 0.0) then
!               runoff=INFILT*0.95 
!
!           else
!               runoff= INFILT*infilt_rate
!  
!           endif
           
           runoff= INFILT*infilt_rate
           
                 
       endif
    else
          runoff=INFILT*0.019
    endif 
    infilt = infilt-runoff
    
!   water redistribution among soil layers
    do i=1,10
        wsc(i)=Amax1(0.00,(wcl(i)-wiltpt)*THKSL(i)*10.0)
        if (do_soilphy) then 
           omegaL(i)=Amax1(0.001,(liq_water(i)*100./thksl(i)-WILTPT)/(FLDCAP-WILTPT))
        else
           omegaL(i)=Amax1(0.001,(wcl(i)-WILTPT)/(FLDCAP-WILTPT))
        endif
            
    enddo
    supply=0.0
    demand=0.0
    do i=1,9
        if(omegaL(i).gt.0.3)then
            supply=wsc(i)*omegaL(i)
            demand=(FLDCAP-wcl(i+1))*THKSL(i+1)*10.0      &
                &               *(1.0-omegaL(i+1))
            exchangeL=AMIN1(supply,demand)
            wsc(i)=wsc(i)- exchangeL
            wsc(i+1)=wsc(i+1)+ exchangeL
            wcl(i)=wsc(i)/(THKSL(i)*10.0)+wiltpt
            wcl(i+1)=wsc(i+1)/(THKSL(i+1)*10.0)+wiltpt
        endif
    enddo
    wsc(10)=wsc(10)-wsc(10)*0.0000
    runoff = runoff+wsc(10)*0.0000
    wcl(10)=wsc(10)/(THKSL(10)*10.0)+wiltpt
!    end of water redistribution among soil layers

!   Redistribute evaporation among soil layers
    Tsrdt=0.0
    DO i=1,10
!   Fraction of SEVAP supplied by each soil layer
    SRDT(I)=EXP(-6.73*(DEPTH(I)-THKSL(I)/2.0)/100.0) !/1.987
!   SRDT(I)=AMAX1(0.0,SRDT(I)*(wcl(i)-wiltpt)) !*THKSL(I))
    Tsrdt=Tsrdt+SRDT(i)  ! to normalize SRDT(i)
    enddo
    do i=1,10
        EVAPL(I)=Amax1(AMIN1(evap*SRDT(i)/Tsrdt,wsc(i)),0.0)  !mm
        DWCL(I)=EVAPL(I)/(THKSL(I)*10.0) !ratio
        wcl(i)=wcl(i)-DWCL(i)
    enddo
    evap=0.0       
    do i=1,10
        evap=evap+EVAPL(I)
    enddo

!   Redistribute transpiration according to root biomass
!   and available water in each layer
    tr_allo=0.0
    do i=1,nfr
        tr_ratio(i)=FRLEN(i)*wsc(i) !*(wcl(i)-wiltpt)) !*THKSL(I))
        tr_allo=tr_allo+tr_ratio(i)
    enddo
    do i=1,nfr
        plantup(i)=AMIN1(transp*tr_ratio(i)/tr_allo, wsc(i)) !mm              
        wupl(i)=plantup(i)/(thksl(i)*10.0)
        wcl(i)=wcl(i)-wupl(i)
    enddo
    transp=0.0
    do i=1,nfr
        transp=transp+plantup(i)
    enddo
    

!******************************************************  
    
!******************************************************    
!   water table module starts here
    vtot = infilt+wsc(1)+wsc(2)+wsc(3)!+wsc(4)+wsc(5)   !total amount of water in top 300mm of soil  mm3/mm2
    !phi = 0.95  !soil porosity   mm3/mm3   the same unit with theta
    
    zmax = 600   !maximum water table depth   mm
    thetasmin = 0.25    !minimum volumetric water content at the soil surface   cm3/cm3
    zthetasmin = 100     !maximum depth where evaporation influences soil moisture   mm
    az = (phi-thetasmin)/zthetasmin     ! gradient in soil moisture resulting from evaporation at the soil surface    mm-1
    
    zwt1 = -sqrt(3.0*(phi*zmax-vtot)/(2.0*az))
    zwt2 = -(3.0*(phi*zmax-vtot)/(2.0*(phi-thetasmin)))
    zwt3 = vtot-phi*zmax                                   
    if ((zwt1 .ge. -100) .and. (zwt1 .le. 0))   zwt = zwt1  !the non-linear part of the water table changing line
    if (zwt2 .lt. -100)                         zwt = zwt2  !the linear part of the water table changing line

!    if ((zwt2 .lt. -100) .and. (zwt2 .ge. -300))zwt = zwt2 !the linear part of the water table changing line valid when Vtot>145mm
!    if (zwt2 .le. -300)                         zwt = -300
    if (phi*zmax .lt. vtot)                     zwt = zwt3  !the linear part when the water table is above the soil surface 
    
!     
!   water table module ends here
!******************************************************    
    
!   Output fwsoil, omega, and topfws
    do i=1,nfr       
        if (do_soilphy) then 
           ome(i)=(liq_water(i)*100./thksl(i)-WILTPT)/(FLDCAP-WILTPT)
        else 
           ome(i)=(wcl(i)-WILTPT)/(FLDCAP-WILTPT)
           ome(i)=AMIN1(1.0,AMAX1(0.0,ome(i)))
        endif 
        fw(i)=amin1(1.0,3.333*ome(i))
    enddo
    
        if (do_soilphy) then 
           topfws=amax1(0.0,topfws)
        else 
           topfws=amin1(1.0,(wcl(1)-WILTPT)/((FLDCAP-WILTPT)))
        endif     
    
    fwsoil=0.0
    omega=0.0
    do i=1,nfr
        fwsoil=fwsoil+fw(i)*frlen(i)
        omega=omega+ome(i)*frlen(i)
    enddo
    
!    write(81,181),rain_t,zwt,wsc(1),wsc(2),wsc(3),wsc(4),wsc(5),wsc(6),wsc(7),wsc(8),wsc(9),wsc(10)
!181	format(12(f11.4,","))   
    return
    end
    
    
!**********************************************************************
!     plant growth model
    subroutine plantgrowth(Tair,omega,GLmax,GRmax,GSmax,&
     &                       LAI,LAIMAX,LAIMIN,SLA,Tau_L,&
     &                       bmleaf,bmroot,bmstem,bmplant,&
     &                       Rootmax,Stemmax,SapS,SapR,&
     &                       StemSap,RootSap,Storage,GDD5,&
     &                       stor_use,onset,accumulation,gddonset,&
     &                       Sps,NSC,fnsc,NSCmin,NSCmax,&
     &                       NSN,CN,CN0,SNgrowth,N_deficit,&
     &                       store,add,L_fall,ht,&
     &                       NPP,alpha_L,alpha_W,alpha_R,&
     &                       RgLeaf,RgStem,RgRoot,Rgrowth)
      implicit none
      real NSC,NSCmin,NSCmax,fnsc,N_deficit
      real CN(8),CN0(8),NSN,nsCN
      real SnscnL,SnscnS,SnscnR
      real store,Storage,GDD5,stor_use,accumulation,gddonset
      integer onset
      real GLmax,GRmax,GSmax,TauLeaf
      real GrowthP,GrowthL,GrowthR,GrowthS
      real Tair,omega,LAI,LAIMAX,LAIMIN,SLA
!     biomass
      real bmleaf,bmroot,bmstem,bmplant,NPP
      real ht,hmax,hl0,CNP0
      REAL LAIMAX0,la0,GPmax,acP,c1,c2
      real Rootmax,Stemmax,SapS,SapR
      real bmL,bmR,bmP,bmS,StemSap,RootSap
      real Rgrowth,Rgroot,Rgleaf,Rgstem
!     scalars
      real St,Sw,Ss,Sn,SL_rs,SR_rs,Slai,Sps,SNgrowth,phiN
      real RS,RS0,RSw
      real gamma_W,gamma_Wmax,gamma_T,gamma_Tmax,gamma_N
      real beta_T,Tcold,Twarm,Topt
      real bW,bT,W
      real L_fall,L_add,add,NL_fall,NL_add,Tau_L
      real alpha_L,alpha_W,alpha_R,alpha_St
      integer i

    Twarm=35.0
    Tcold=5.0
!    Tcold=0.0       ! For SPRUCE
    Topt=30.
    phiN=0.33

    bmL=bmleaf*0.48   ! Carbon
    bmR=bmRoot*0.48
    bmS=bmStem*0.48

    if(bmL.lt.NSC/0.333)bmL=NSC/0.333
    if(bmR.lt.NSC/0.333)bmR=NSC/0.333
    if(bmS.lt.NSC/0.334)bmS=NSC/0.334
    StemSap=SapS*bmS  ! Weng 12/05/2008
    RootSap=SapR*bmR
    if(StemSap.lt.0.001)StemSap=0.001
    if(RootSap.lt.0.001)RootSap=0.001

    bmP=bmL+bmR+bmS					! Plant C biomass 
    acP=bmL+StemSap+bmS					! Plant available sapwood C  
    CNp0=bmP/(bmL/CN0(1)+bmR/CN0(3)+bmS/CN0(2))		! Plant CN ratio

    hmax=24.19   ! m
    hl0=0.00019  ! m2/kg C
    LAIMAX0=6.
    la0=0.2
    ht=hmax*(1.-exp(-hl0*bmP))				! Scaling plant C biomass to height
    LAIMAX=AMAX1(LAIMAX0*(1.-exp(-la0*ht)),LAIMIN+0.1)  ! Scaling plant height to maximum LAI

!   Phenology
    if((GDD5.gt.gddonset).and.onset.eq.0.and.storage.gt.stor_use) then
        onset=1
    endif
    if((onset.eq.1).and.(storage.gt.stor_use))then
        if(LAI.lt.LAIMAX)add=stor_use
        storage=storage-add
    else
        add=0.0
        onset=0
    endif
    if(accumulation.lt.(NSCmax+0.005*RootSap))then
        store=AMAX1(0.,0.005*NSC)			! 0.5% of nonstructure carbon is stored
    else
        store=0.0
    endif
    accumulation=accumulation+store

!   Scalars for plant growth
!      Sps=Amin1(1.0,3.33*AMAX1(0.0,1.0 - fnsc))
      Sps=Sps*(1.-exp(-phiN*NSN))							! Sps is not assigned previous, something is wrong. -JJJJJJJJJJJJJJJJJJJJJ
      Ss=AMIN1(1.0,2.*fnsc)
      RS0=1.0
      RS=bmR/bmL
      SL_rs=RS/(RS+RS0*(2.-W))
      SR_rs=(RS0*(2.-W))/(RS+RS0*(2.-W))
      Slai=amin1(1.0,2.333*(LAIMAX-LAI)/(LAIMAX-LAIMIN))
      St=AMAX1(0.0, 1.0-exp(-(Tair-gddonset/10.)/5.0))  !0.5 !
      Sw=AMAX1(0.333, 0.333+omega)
      W = AMIN1(1.0,3.333*omega)

!     Plant growth and allocation, based on LM3V
      GPmax=(GLmax*bmL+GSmax*StemSap+GRmax*bmR) !/acP					! Need to divide acP, g C day-1 -JJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
      GrowthP=AMIN1(GPmax*fnsc*St*(1.-exp(-NSN)),  & ! 
     &              0.004*NSC,&
     &              0.004*NSN*CNp0)
 
     !      c1=(bmR+200.)/bmL*CN(1)/CN0(1) !+N_deficit/NSN
!      c1=bmL/bmR*CN(1)/CN0(1) !+N_deficit/NSN
!      c2=0.5*250e3*SLA*0.00021*ht*2.
      !write(*,*)LAI
      !GrowthL=MAX(0.0,MIN(GrowthP*0.5,0.05*(LAIMAX-LAI)/SLA))    ! 1./(1.+c1+c2)
      GrowthL=MAX(0.0,GrowthP*0.5) 
      GrowthR=MIN(GrowthP*0.4,MAX(0.0,0.75/Sw*bmL-bmR))  ! *c1/(1.+c1+c2)
      GrowthS=MAX(0.0,GrowthP - (GrowthL+GrowthR) )         ! *c2/(1.+c1+c2)

      NPP = GrowthL + GrowthR + GrowthS + add       ! Modified by Jiang Jiang 2015/10/13

      if(NPP.eq.0.0)then
            alpha_L=0.333
            alpha_W=0.333
            alpha_R=0.333
      else
            alpha_L=(GrowthL+add)/NPP     
            alpha_W=GrowthS/NPP
            alpha_R=GrowthR/NPP
      endif
!     Carbon cost for growth
!     Rgrowth,Rgroot,Rgleaf,Rgstem, 0.5 is from IBIS and Amthor, 1984
      Rgleaf=0.5*GrowthL
      Rgstem=0.5*GrowthS
      Rgroot=0.5*GrowthR
      Rgrowth=Rgleaf+Rgstem+Rgroot

!     Leaf litter 

      gamma_Wmax=0.12/24. ! maxmum leaf fall rate per hour
      gamma_Tmax=0.12/24.

      bW=4.0
      bT=2.0

      if(Tair.gt.(Tcold+10.)) then
            beta_T=1.
      else 
            if(Tair.gt.Tcold)beta_T=(Tair-Tcold)/10.
            if(Tair.LE.Tcold)beta_T=0.0
      endif

      if (tau_L < 8760.)then
                gamma_W=(1. - W)     **bW  * gamma_Wmax
                gamma_T=(1. - beta_T)**bT * gamma_Tmax
      else
                gamma_W=0.
                gamma_T=0.
      endif
      gamma_N=1.0/Tau_L*Sw      ! Modify by Jiang Jiang 2015/10/20
      if(LAI < LAIMIN) then
            gamma_W=0.
            gamma_T=0.
            gamma_N=0.
      endif
    !  L_fall=bmleaf*0.48*AMIN1((gamma_T+gamma_N),0.99)
      L_fall=bmleaf*0.48*gamma_N

      return
      end

!************************************************************************
!     carbon transfer according to Xu et al. 2007 
      subroutine TCS_CN(Tair,Tsoil,omega,runoff,&
     &               NPP,alpha_L,alpha_W,alpha_R,L_fall,&
     &               tauC,QC,OutC,Rh_pools,Rnitrogen,NSC,&
     &               CNmin,CNmax,NSNmax,NSNmin,alphaN,    &        ! nitrogen
     &               NSN,N_uptake,N_miner,QN,QNminer,&
     &               CN,CN0,fnsc,rdepth,N_deficit,&
     &               N_leaf,N_wood,N_root,N_LF,N_WF,N_RF,&
     &               N_deposit,N_fixation,N_leach,N_vol,&
     &               SNvcmax,SNgrowth,SNRauto,SNrs,Q10,tsoill,testout,do_soilphy) 
         
      implicit none
      real NPP,NPP_L,NPP_W,NPP_R
      real L_fall,L_add,LAI,SLA,rdepth
      real Tair,Tsoil,omega,runoff
!     allocation ratios
      real alpha_L,alpha_W,alpha_R
!     pools
      real Q_plant,QC(8),TauC(8),OutC(8)
      real etaL,etaW,etaR                ! the percentage of fine litter of the litters from plant parts
      real f_F2M,f_C2M,f_C2S,f_M2S,f_M2P,f_S2P,f_S2M,f_P2M
      real Rh_pools(5),Q10h(5),Q10 ! Q10 of the litter and soil C pools
!     the fraction of C-flux which enters the atmosphere from the kth pool
      real f_CO2_fine,f_CO2_coarse,f_CO2_Micr,f_CO2_Slow,f_CO2_Pass
!     for nitrogen sub-model
      real CN0(8),CN(8),OutN(8),QN(8),QNminer,QNplant
      real CNmin,CNmax,NSNmax,NSNmin,NSN
      real CN_plant,CN_foliage
      real N_demand,N_deficit,N_immob,N_imm(5),N_fixation,Nfix0
      real N_transfer,N_miner,N_uptake,N_deposit,N_loss,N_leach,N_vol
      real alphaN,Qroot0,Cfix0
      real Scalar_N_flow,Scalar_N_T
      real N_leaf,N_wood,N_root
      real N_LF,N_WF,N_RF
      real NSC,fnsc,ksye
      real SNvcmax,SNgrowth,SNRauto,SNrs
      real kappaVcmax
      real SNfine,SNcoarse,SNmicr,SNslow,SNpass
      real Rnitrogen,costCuptake,costCfix,costCreuse
      real Creuse0,Nup0,N_deN0,LDON0
!     the variables relative to soil moisture calcualtion
      real S_omega    !  average values of the moisture scaling functions
      real S_t(5)     !  average values of temperature scaling functions
      real S_w_min    !  minimum decomposition rate at zero available water
!     For test
      real totalC1,totalN1,C_in,C_out,N_in,N_out,totalC2,totalN2
      real ScNloss
      real tsoill(10),frac_soc(10),testout(11),tsoil_layer(11)

      integer i,j,k,n,m
      integer day,week,month,year
      
      logical do_soilphy
      
      tsoil_layer = testout 

!     temperature sensitivity of Rh_pools
!      Q10h=(/2.0,2.0,2.0,2.0,2.0/)  ! for Oak Ridge
      Q10h=(/Q10,Q10,Q10,Q10,Q10/)
      frac_soc=(/0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1/)
      
      Qroot0=500.
      Nfix0=1./60.   ! maximum N fix ratio, N/C
      Nup0 =0.02     ! nitrogen uptake rate
      Cfix0=12.      ! C cost per N for fixation
      ksye=0.05      ! C cost per N for uptake
      Creuse0=2.     ! C cost per N for resorption
      ScNloss=1.
      N_deN0=1.E-3*ScNloss   ! 1.E-3, 5.E-3, 10.E-3, 20.E-3
      LDON0=1.E-3*ScNloss

      Rnitrogen=0.
!     for N scalars
      CNmin=40.0
      CNmax=200.0
!     Max and min NSN pool
      NSNmax = QN(1) + 0.2*QN(2) + QN(3)  ! 15.0
      NSNmin = 0.01
!     partitioning coefficients
      etaL=0.6          ! 60% of foliage litter is fine, didn't use
      etaW=0.15        ! 15% of woody litter is fine
      etaR=0.85         ! 85% of root litter is fine  , didn't use    
      f_F2M=0.55        ! *exp((CN0(4)-CN_fine)*0.1)
      f_C2M=0.275       ! *exp((CN0(5)-CN_coarse)*0.1)
      f_C2S=0.275       ! *exp((CN0(5)-CN_coarse)*0.1)
      f_M2S=0.3
      f_M2P=0.1
      f_S2P=0.2        !0.03 Change by Jiang Jiang 10/10/2015
      f_S2M=0.5
      f_P2M=0.45

!     calculating soil scaling factors, S_omega and S_tmperature
      S_w_min=0.08 !minimum decomposition rate at zero soil moisture
      S_omega=S_w_min + (1.-S_w_min) * Amin1(1.0,2.0*omega)

      if (do_soilphy) then 
          S_t=(/0.0,0.0,0.0,0.0,0.0/)
          do i=1,5
             if(i.lt.3) then
!               S_t(i)=Q10h(i)**((Tsoil-5.)/10.)  ! Oak
                S_t(i)=Q10h(i)**((tsoil_layer(1)-10.)/10.)  ! Duke
             else 
                do j=1,10
                    S_t(i)=S_t(i)+frac_soc(j)*Q10h(i)**((tsoil_layer(j+1)-10.)/10.)  ! Duke
                enddo
             endif
          enddo
      else
          do i=1,5
!            S_t(i)=Q10h(i)**((Tsoil-5.)/10.)  ! Oak
             S_t(i)=Q10h(i)**((Tsoil-10.)/10.)  ! Duke
          enddo  
      endif 


!     Calculating NPP allocation and changes of each C pool
      NPP_L = alpha_L * NPP           ! NPP allocation
      NPP_W = alpha_W * NPP
      NPP_R = alpha_R * NPP
!     N scalars on decomposition
      SNfine  =exp(-(CN0(4)-CN(4))/CN0(4)) 
      SNcoarse=exp(-(CN0(5)-CN(5))/CN0(5)) 
      SNmicr  =exp(-(CN0(6)-CN(6))/CN0(6)) 
      SNslow  =1. !exp(-(CN0(7)-CNC(7))/CN0(7)) 
      SNpass  =exp(-(CN0(8)-CN(8))/CN0(8)) 
      
!     the carbon leaving the pools
      OutC(1)=L_fall
      OutC(2)=QC(2)/tauC(2)*S_omega !*exp(CN(2)/CN0(2)-1.) 
      OutC(3)=QC(3)/tauC(3)*S_omega

      OutC(4)=QC(4)/tauC(4)*S_omega* S_T(1)*CN(4)/CN0(4)!*SNfine
      OutC(5)=QC(5)/tauC(5)*S_omega* S_T(2)*CN(5)/CN0(5)!*SNcoarse
      OutC(6)=QC(6)/tauC(6)*S_omega* S_T(3)!*SNmicr
      OutC(7)=QC(7)/tauC(7)*S_omega* S_T(4)!*SNslow
      OutC(8)=QC(8)/tauC(8)*S_omega* S_T(5)!*SNpass

!     heterotrophic respiration from each pool
      Rh_pools(1)=OutC(4)* (1. - f_F2M)
      Rh_pools(2)=OutC(5)* (1. - f_C2M - f_C2S)
      Rh_pools(3)=OutC(6)* (1. - f_M2S - f_M2P)
      Rh_pools(4)=OutC(7)* (1. - f_S2P - f_S2M)
      Rh_pools(5)=OutC(8)* (1. - f_P2M)
      
      !write(*,*)"out",OutC(4)+OutC(5)+OutC(6)+OutC(7)+OutC(8)
      
!      testout(1)=Rh_pools(1)
!      testout(2)=Rh_pools(2)
!      testout(3)=Rh_pools(3)
!      testout(4)=Rh_pools(4)
!      testout(5)=Rh_pools(5)
      
!========================================================================
!     Nitrogen part
!     nitrogen leaving the pools and resorption
      do i=1,8
         OutN(i) = OutC(i)/CN(i)
      enddo

!     nitrogen mineralization
      N_miner=OutN(4)* (1. - f_F2M)  &
     &       +OutN(5)* (1. - f_C2M - f_C2S) &
     &       +OutN(6)* (1. - f_M2S - f_M2P) &
     &       +OutN(7)* (1. - f_S2P - f_S2M) &
     &       +OutN(8)* (1. - f_P2M)

!     Nitrogen immobilization
      N_imm=0.
      N_immob=0.
      if(QNminer>0)then
          do i=4,8
            if(CN(i)>CN0(i))then
                N_imm(i-3)=Amin1(QC(i)/CN0(i)-QC(i)/CN(i)  &
     &             ,0.1*QNminer)
                N_immob=N_immob+N_imm(i-3)
            endif
          enddo
      endif

!     Let plant itself choose the strategy between using C to uptake
!     or fix N2 by comparing C invest.
!     N demand
      N_demand=NPP_L/CN0(1)+NPP_W/CN0(2)+NPP_R/CN0(3) !+N_deficit
!     Nitrogen input:
      N_transfer=0.
	N_uptake=0.
	N_fixation=0.
      costCuptake=0.
      costCfix=0.
      costCreuse=0.
!     1. Nitrogen resorption
      N_transfer=(OutN(1) + OutN(2) +OutN(3))*alphaN
      costCreuse= Creuse0*N_transfer
      N_demand=N_demand-N_transfer

      If(N_demand>0.0)then

!     2.  N uptake
          if(ksye/QNminer<Cfix0)then
              N_uptake=AMIN1(N_demand+N_deficit,      &
     &                       QNminer*QC(3)/(QC(3)+Qroot0), &
     &                       Nup0*NSC/(ksye/QNminer))
              costCuptake=N_uptake*(ksye/QNminer)
              N_demand=N_demand-N_uptake
          elseif(NSN<24.*30.*N_demand)then
!     3.  Nitrogen fixation
              N_fixation=Amin1(N_demand,fnsc*Nfix0*NSC)
              costCfix=Cfix0*N_fixation
              N_demand=N_demand-N_fixation
          endif
      endif
      N_deficit=N_deficit+N_demand

!     update NSN
      NSN=NSN+N_transfer+N_uptake+N_fixation
!     Total C cost for nitrogen
      Rnitrogen=costCuptake+costCfix+costCreuse

!      Oak Ridge N fixation rate: 
!      asymbiotic: 2 mg N/m2/yr ;  symbiotic: 65 mg/m2/yr, Oak Ridge
!      N_fixation=0.067/8760. ! Oak Ridge
!      N_fixation=0.23/8760.  ! Duke

!     Nitrogen using, non-structural nitrogen pool, NSN
      N_leaf =AMIN1(NPP*alpha_L/CN(1)+QC(1)/CN0(1)-QC(1)/CN(1),0.2*NSN)
      N_wood =AMIN1(NPP*alpha_W/CN(2)                         ,0.1*NSN)
      N_root =AMIN1(NPP*alpha_R/CN(3)+QC(3)/CN0(3)-QC(3)/CN(3),0.2*NSN)
      NSN=NSN-(N_leaf+N_wood+N_root)

      N_LF=OutN(1)*(1.-alphaN)
      N_WF=OutN(2)*(1.-alphaN)
      N_RF=OutN(3)*(1.-alphaN)

!     update QNminer
      QNminer=QNminer+N_miner+N_deposit  &
     &              -(N_uptake+N_immob)

!     Loss of mineralized N and dissolved organic N
      Scalar_N_flow=0.5*runoff/rdepth
!      Scalar_N_T=0.005*(Tsoil+273.)/(Tsoil+273+333.)
      if (do_soilphy) then 
          Scalar_N_T = 0.0 
          do j=1,10
              Scalar_N_T = Scalar_N_T + frac_soc(j)*N_deN0*exp((tsoil_layer(j+1)-25.)/10.)  
          enddo

      else
          Scalar_N_T=N_deN0*exp((Tsoil-25.)/10.)
      endif 
      
      N_leach=Scalar_N_flow*QNminer+Scalar_N_flow*QN(6)*LDON0
      N_vol  =Scalar_N_T*QNminer
      N_loss =N_leach + N_vol

!     update QNminer
      QNminer=QNminer-N_loss
!     update plant carbon pools, ! daily change of each pool size
      QC(1)=QC(1) - OutC(1) + NPP_L    
      QC(2)=QC(2) - OutC(2) + NPP_W
      QC(3)=QC(3) - OutC(3) + NPP_R
      QC(4)=QC(4) - OutC(4) + OutC(1)+etaW*OutC(2)+OutC(3)
      QC(5)=QC(5) - OutC(5) + (1.-etaW)*OutC(2)
      QC(6)=QC(6) - OutC(6) + f_F2M*OutC(4)+f_C2M*OutC(5)     &         
     &                      + f_S2M*OutC(7)+f_P2M * OutC(8)
      QC(7)=QC(7) - OutC(7)+f_C2S*OutC(5)+f_M2S*OutC(6)
      QC(8)=QC(8) - OutC(8)+f_M2P*OutC(6)+f_S2P*OutC(7)

      Q_plant =QC(1) + QC(2) + QC(3)
!     update nitrogen pools
      QN(1)=QN(1) - OutN(1) + N_leaf
      QN(2)=QN(2) - OutN(2) + N_wood
      QN(3)=QN(3) - OutN(3) + N_root
      QN(4)=QN(4) - OutN(4)+ N_imm(1)     &
     &            + (OutN(1) + etaW*OutN(2) + OutN(3))*(1.-alphaN)
      QN(5)=QN(5) - OutN(5) + N_imm(2)   &
     &            + (1.-etaW)*OutN(2)*(1.-alphaN)

      QN(6)=QN(6) - OutN(6) + N_imm(3) - Scalar_N_flow*QN(6)*LDON0  &
     &            + f_F2M*OutN(4)+f_C2M*OutN(5)  &
     &            + f_S2M*OutN(7)+f_P2M*OutN(8)
      QN(7)= QN(7) - OutN(7) + N_imm(4)  &
     &                         + f_C2S*OutN(5) &
     &                         + f_M2S*OutN(6)
      QN(8)= QN(8) - OutN(8) + N_imm(5) &
     &         + f_M2P*OutN(6) + f_S2P*OutN(7)
      QNplant = QN(1) + QN(2)+ QN(3)

!     update C/N ratio
      CN=QC/QN
      CN_foliage=(QC(1)+QC(3))/(QN(1)+QN(3))

!     calculate N related scalars for Oak Ridge
!      SNvcmax =exp(-(CN(1)-CN0(1))) ! /CN0(1) ! Oak Ridge
!      SNgrowth=exp(-(CN(1)-CN0(1))/CN0(1)) !  AMAX1((CNmax-CN_foliage)/(CNmax-CNmin),0.0)+0.25
!      SNRauto =AMAX1((CNmax-CN_foliage)/(CNmax-CNmin),0.0)+0.5
!      SNrs=1.

!     calculate N related scalars for Duke FACE
      kappaVcmax=CN0(1)/1.
      SNvcmax =exp(-kappaVcmax*(CN(1)-CN0(1))/CN0(1)) ! /CN0(1) ! Duke
      SNgrowth=exp(-(CN(1)-CN0(1))/CN0(1)) !  AMAX1((CNmax-CN_foliage)/(CNmax-CNmin),0.0)+0.25
      SNRauto =exp(-(CN(1)-CN0(1))/CN0(1)) !  AMAX1((CNmax-CN_foliage)/(CNmax-CNmin),0.0)+0.5
      SNrs=1.

      return
      end
      
      subroutine snow_d(rain,lat,days,ta,snow_dsim,fa,fsub,rho_snow,melt,dcount,decay_m)
          
       real lat,tr,daylength,dec,melt,fa,sublim,dsnow,snow_in,decay_m,fsub
       real rain,snow_dsim,rho_snow,d_count
       integer days
       real snow_dsim_pre
       
!       rho_snow =100.
!       fa=0.1
!       fsub=0.1
       tr=0.0174532925
       
       dec=sin(((real(days)-70.)/365.)*360.*tr)*23.44
       daylength=acos(-tan(lat*tr)*tan(dec*tr))/7.5 
       daylength=daylength/tr/24.
            
       if (snow_dsim .ge. 0.) then
           dcount = dcount +1.
       else 
           dcount =0.
       endif
       
       sublim=0.
       if (ta .gt. 0. .and. snow_dsim .gt. 0.) sublim=fsub*715.5*daylength*esat(ta)/(ta+273.2)*0.001   ! 0.001 from Pa to kPa
       !if (sublim .lt. 0.1) sublim=0.
       !sublim=0.
       !sublim=AMIN1(sublim,0.2)
       !melt=fa*(2.63+2.55*ta+0.0912*ta*rain)
       
       !if (snow_dsim .gt. 0.7) sublim = 10.
       melt=0.
       if (ta .gt. 0. .and. snow_dsim .gt. 0.) melt=fa*(2.63+2.55*ta+0.0912*ta*rain)
       !if (ta .gt. 0. .and. snow_dsim .gt. 0.) melt=fa*(0.55*ta)
       
       if (dcount .gt.0. .and. ta .lt.5.) then
           !write(*,*)'melt_befor',melt
           melt=melt*EXP(-decay_m*dcount/365.)
           !write(*,*)'melt_after',melt
       endif

       !write(*,*),EXP(-3.*dcount/365.)
       
       !melt=AMIN1(melt, snow_dsim*rho_snow-sublim)
       !if (melt .lt. 2.) melt=0.
       
       if (ta .le. 0.) snow_in =rain
       
       dsnow=snow_in-sublim-melt 
       
       snow_dsim_pre = snow_dsim
       snow_dsim =snow_dsim + dsnow/rho_snow 
       
       if (snow_dsim .le. 0.0) then 
          snow_dsim=0.0 
          melt = snow_dsim_pre*rho_snow +snow_in-sublim    !! for water part
       endif 
       melt=AMAX1(melt, 0.)
       
 !      write(*,*)'snow_dsim',snow_dsim
!       write(*,*)'snow_dsim',snow_dsim
      ! write(*,*)'doy',pi   
       
       write(88,188),melt,snow_dsim,snow_in,ta
188    format(12(f11.4,","))  
    
    return
    end
      
!    !    ========================================================================================
      subroutine Tsoil_simu(Rsoilab1,Rsoilab2,QLleaf,QLair,Tair,Dair,&
&         fbeam,FLAIT,sigma,emsoil,rhoS,Rconst,&
&         extkd,extkb,cpair,Patm,AirMa,H2OMw,&
&         H2OLv0,wcl,raero,wsmax,wsmin,wind,sftmp,Tsoill,testout,ht,ice,&
&         snow_depth,Tsnow,Twater,Tice,water_tw,ice_tw,diff_s,G,tsoil,&
&         diff_snow,albedo_snow,resht,thd_snow_depth,thksl,zwt,Esoil,Hsoil,liq_water,&
&         shcap_snow,condu_snow,condu_b,depth_ex,dcount_soil)          
      implicit none 
      integer i
      real tsoil
      real Rsoilab1,Rsoilab2,qlleaf,qlair,tair,Dair,fbeam,flait,sigma,emsoil
      real rhoS(3),rconst,extkd,extkb,cpair,patm,airma,h2omw,h2olv0,raero,wsmax,wsmin
      real esoil,G,hsoil,wind,ht,esat,theta_sat_min,Rsoilabs,Rsoil,difsv2,difsv1
      real delta
      !real thksl(10)
      real TairK,H2OLv
      real ice(10)
      real wcl(10),thksl(10)
      real,dimension(10):: Tsoill,liq_water
      real,dimension(11)::testout
      real WILTPT,FILDCP,temph1,temph2     
      real sftmp,hitmax,rflen,zopnd,thkns1,thkns2
      real Twater, flux_water,Tsnow,flux_snow,Tice
      real condu_water,shcap_water,shcap_ice,shcap_snow,condu_snow,depth_ex
      real albedo_snow, albedo_water,ice_incr,heat_excess,heat_adjust,ice_tw,water_tw
      real inter_var,latent_heat_fusion,QLsoil,Rsoilab3
      real rhocp,slope,psyc,Cmolar,fw1,resoil,rLAI,resht,resdh,dnr,dsh,dgh,dle,drsdh
      real f_om,theta_sat_om,b_om,b_min,phi_om,phi_min,theta_sat,b_tot,phi_sat,gravi
      real water_table_depth,snow_depth,temph_water,temph_snow
      real condu_air,shcap_air,condu(10), shcap(10), condu_ice,tsoill_pre, thd_t
      real ice_density,condu_soil,shcap_soil 
      real thd_snow_depth,resht_lai,zwt,snow_depth_t
      real diff_s, diff_snow,condu_s,tsoill_0,diff_air,d_cor,condu_b,dcount,dcount_soil
      real sftmp_pre
      integer n_layers
      
      real, allocatable ::depth_z(:) 
      n_layers=10
      allocate(depth_z(n_layers))
         
      
      !write(*,*),thd_t
      ! soil thermal conductivity W m-2 K-1
      ice_density=916.!916.
      thkns1=thksl(1)/2.
      shcap_ice=2117.27*ice_density
      condu_ice=2.29
      condu_water=0.56!0.56
      shcap_water=4188000.
      condu_soil=0.25
      shcap_soil=2600000.
      condu_s=0.25
      thd_t=0.0
      
      diff_snow=3600.*condu_snow/shcap_snow*10000.
      diff_s=3600.*condu_b/shcap_soil*10000.
           
      latent_heat_fusion = 333700.   ! j kg-1
      condu_air=0.023
      shcap_air=1255.8
      
      diff_air=3600.*condu_air/shcap_air*10000. 
           
      water_tw=zwt*0.001-ice_tw
      water_table_depth=zwt*0.1
      
      snow_depth_t = snow_depth - 0.46*0.0     ! warming in Tair impact on snow_depth
                                               ! in unit cm 0.46 based on snow_depth vs. tair regression     
      if (snow_depth_t .lt. thd_snow_depth) snow_depth_t =0.0
      
       if (snow_depth_t .gt. 0.) then
           dcount_soil = dcount_soil +1./24.
       else 
           dcount_soil =0.
       endif

      if (water_table_depth .lt. 4. .and. water_table_depth .gt. 0.0) water_table_depth =0.    ! avoid numerical issues when 
      if (water_table_depth .lt. -99.) water_table_depth =-30.    ! temporary for NaN    
        
      albedo_water =0.1      
      ! soil water conditions
      WILTPT=wsmin/100.
      FILDCP=wsmax/100.
      TairK=Tair+273.2   
         
      flux_snow = 0.0 
      
      depth_z=(/0., 0., 0., 0., 0., 0., 0.,0.,0.,0./) 
   
!      if (snow_depth_t .gt. 0.0) then 
!          emsoil =0.98
!      elseif (water_table_depth .gt. 0.0) then
!          emsoil =0.99
!      endif
      
      QLsoil=emsoil*sigma*((sftmp+273.2)**4)
      Rsoilab3=(QLair+QLleaf)*(1.0-rhoS(3))-QLsoil       

   ! Total radiation absorbed by soil
      if (snow_depth_t .gt. 0.0) then 
         Rsoilabs=(Rsoilab1+Rsoilab2)*(1-albedo_snow)/(1-0.1)+Rsoilab3  
      elseif (water_table_depth .gt. 0.0) then 
         Rsoilabs=(Rsoilab1+Rsoilab2)*(1-albedo_water)/(1-0.1)+Rsoilab3  
      else
      Rsoilabs=Rsoilab1+Rsoilab2+Rsoilab3
      endif
      
!    thermodynamic parameters for air
      rhocp=cpair*Patm*AirMa/(Rconst*TairK)      
      H2OLv=H2oLv0-2.365e3*Tair
      slope=(esat(Tair+0.01)-esat(Tair))/0.01   
   
      psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar=Patm/(Rconst*TairK)
      fw1=AMIN1(AMAX1((FILDCP-wcl(1))/(FILDCP-WILTPT),0.3),1.0)
!      
      if (water_table_depth .gt. 0.0) then 
         Rsoil = 0. 
      else 
         Rsoil=30.*exp(0.2/fw1)
      endif 
      !Rsoil=40.
      !Rsoil=5.
      rLAI=exp(FLAIT)     
!     latent heat flux into air from soil
!           Eleaf(ileaf)=1.0*
!     &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    !2* Weng 0215
!     &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
!
      Esoil=(slope*(Rsoilabs-G)+rhocp*Dair/(raero+rLAI))/       &
     &      (slope+psyc*(Rsoil/(raero+rLAI)+1.))
     
     !if (snow_depth_t .gt. 0.) Esoil=0.
          
!      endif
     !Esoil=0.
     resht_lai=resht*FLAIT
     !resht_lai= resht*exp(FLAIT)/15. ! need improvement, should be a function of LAI 
      !if (water_table_depth .gt. 0.0) resht_lai=resht/FLAIT*0.2

      !resht_lai=200.      
      Hsoil=rhocp*(sftmp-Tair)/resht_lai
      
      !Hsoil=rhocp*(sftmp-Tair)/resht_lai
      !Hsoil=1010.*1.17*(sftmp-Tair)/resht_lai
      i=1;
      condu(i)=(FILDCP-wcl(i))*condu_air+liq_water(i)/(thksl(i)*0.01)*condu_water+ &
         &  ice(i)/(thksl(i)*0.01)*condu_ice +(1-FILDCP)*condu_soil
      shcap(i)=(FILDCP-wcl(i))*shcap_air+liq_water(i)/(thksl(i)*0.01)*shcap_water+ &
            ice(i)/(thksl(i)*0.01)*shcap_ice +(1-FILDCP)*shcap_soil
      difsv1=3600.*condu(i)/shcap(i)*10000.
          
      G=condu(1)*(sftmp-tsoill(1))/(thksl(1)/2.*0.01)
      if (snow_depth_t .gt. 0.0) then 
          G=condu_snow*(sftmp-Tsnow)/(snow_depth_t/2.*0.01)
      endif
      
      ! thksl(1)
      !G=0.   
      ! Residual heat energy.
      RESDH=Rsoilabs-Hsoil-Esoil-G
      !G=RESDH
      
      ! First derivative of net radiation; sensible heat; ground heat;
      DNR=4.*emsoil*sigma*(sftmp+273.2)**3
      DSH=rhocp/resht_lai 
      DGH=condu_s/(thksl(1)/2.*0.01)
      DLE=(DNR+DGH)*slope/(slope+psyc*(Rsoil/(raero+rLAI)+1.))      
      drsdh=-DNR-DSH-DGH-DLE
    ! Calculate increment DELTA.
      DELTA=resdh/drsdh
      sftmp_pre=sftmp
      sftmp=sftmp-DELTA
      if (ABS(sftmp_pre -sftmp) .gt. 20. ) sftmp=sftmp_pre
      
      tsoill_0=sftmp
    ! Temperature dynamics along soil profile
     !difsv1=diff_s
             
!     if (snow_depth_t .gt. 0.) then
!        d_cor=20.
!        diff_snow=d_cor*diff_snow/snow_depth_t
!        write(*,*),snow_depth_t
!     endif

     do i=1,10
        Tsoill_pre=tsoill(i) 

             
        if (water_table_depth .lt. 0.0 .and. -water_table_depth .lt. depth_z(i)) then
            liq_water(i)=FILDCP*thksl(i)*0.01-ice(i)
        else
            liq_water(i)=wcl(i)*thksl(i)*0.01-ice(i)
        endif

                       
        if (i .eq. 1) then 
            depth_z(1)=thksl(1)
        else 
            depth_z(i)=depth_z(i-1)+thksl(i)
        endif
        
         thkns2=(thksl(i)+thksl(i+1))/2.
                 
        if (i .eq. 10) then
         difsv2=3600.*condu(i)/shcap(i)*10000. 
        else
         condu(i+1)=(FILDCP-wcl(i+1))*condu_air+liq_water(i+1)/(thksl(i+1)*0.01)*condu_water+ &
         &  ice(i+1)/(thksl(i+1)*0.01)*condu_ice +(1-FILDCP)*condu_soil
         shcap(i+1)=(FILDCP-wcl(i+1))*shcap_air+liq_water(i+1)/(thksl(i+1)*0.01)*shcap_water+ &
            ice(i+1)/(thksl(i+1)*0.01)*shcap_ice +(1-FILDCP)*shcap_soil 

         difsv2=3600.*condu(i+1)/shcap(i+1)*10000.
        endif
        
        temph2=(difsv1+difsv2)*(Tsoill(i)-Tsoill(i+1))/thkns2 

         
        !(dcount_soil/150.)**3.*
      !!!!!!!!!!!!!!!!!!!! start first layer !!!!!!!!!!!!!!!!!!!!!!
        !!!!!! adjust if there are snow or water layer above !!!!!!!!!!!!!!!!!!!!
        if(i.eq.1) then
           if (snow_depth_t .gt. 0.) then   
               temph_snow = Amin1(diff_snow,difsv1)*(Tsnow-Tsoill(1))/((snow_depth_t+thksl(1))/2.)
               Tsnow=Tsnow+(exp(-depth_ex*snow_depth_t)*diff_snow*(sftmp-Tsnow)/(snow_depth_t/2.) &
        &              -temph_snow)/(snow_depth_t/2.+(snow_depth_t+thksl(1))/2.) 
            
               Tsoill(1)=Tsoill(1)+(temph_snow &
        &              -temph2)/((snow_depth_t+thksl(1))/2.+thkns2) 
        
               if (Tsnow .gt.0.0) then 
                   Tsnow =0.0   
                   Tsoill(1)=0.
               endif
                
                !write(*,*),temph2
               drsdh =0.0    ! temporarily set drsdh =0 for heat adjustment of soil when  
               tsoill_0= (Tsoill(1)+Tsnow)/2.
           elseif (water_table_depth .gt. 0.) then  
               temph_water = (3600.*condu_water/shcap_water*10000.+difsv1)*(Twater-Tsoill(1))/((water_table_depth+thksl(1))/2.)! there is snow layer 
               Twater=Twater+(2.*3600.*condu_water/shcap_water*10000.*(sftmp-Twater)/(water_table_depth/2.) &
        &              -temph_water)/(water_table_depth/2.+(water_table_depth+thksl(1))/2.) 
    
         !!!!!!!!!!!!!!!!!!  Phase change surface water !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               if (Twater .lt. 0.0 .and. water_tw .gt. 0.0) then  ! freeze 
                   heat_excess=-(shcap_water/360000.*water_tw*100.-drsdh)*Twater
                   ice_incr=heat_excess*3600./latent_heat_fusion/ice_density
                   !write(*,*)'water_tw',water_tw
                   if (ice_incr .lt. water_tw) then
                     ice_tw=ice_tw +ice_incr
                     water_tw=water_tw-ice_incr
                     Twater=0.0
                     Tice=0.0
                   else
                     ice_tw=ice_tw +water_tw
                     water_tw=0.0
                     Tice = Tice - latent_heat_fusion*(ice_incr-water_tw)*ice_density/(shcap_ice*ice_tw)
                   endif     
               elseif (Twater .gt. 0.0 .and. ice_tw .gt. 0.0) then    ! thraw              
                   heat_excess=(shcap_water/360000.*ice_tw*100.-drsdh)*Twater
                   ice_incr=heat_excess*3600./latent_heat_fusion/ice_density
                   if (ice_incr .lt. ice_tw) then
                     ice_tw=ice_tw -ice_incr
                     water_tw=water_tw+ice_incr
                     Twater=0.0
                     Tice=0.0
                   else
                     water_tw=water_tw +ice_tw
                     ice_tw=0.0
                     Twater = Twater + latent_heat_fusion*(ice_incr-ice_tw)*ice_density/(shcap_water*water_tw)
                   endif
              !write(*,*)'heat_excess',ice_incr-ice_tw
               endif                       
   !!!!!!!!!!!!!!!!!!!!!!!!! end of phase change for surface layer !!!!!!!!!!!!!!!!!!!  

               temph2=(difsv1+3600.*condu_water/shcap_water*10000.)*(Tsoill(i)-Tsoill(i+1))/thkns2 
               if (water_tw .eq. 0.0 .and. ice_tw .gt. 0.0) then 
                   Tsoill(1)=Tsoill(1)+(2.*3600.*condu_ice/shcap_ice*10000.*(Tice-Tsoill(1))/thkns1 &
        &              -temph2)/(thkns1+thkns2) 
               else 
                   Tsoill(1)=Tsoill(1)+(2.*3600.*condu_water/shcap_water*10000.*(Twater-Tsoill(1))/thkns1 &
        &              -temph2)/(thkns1+thkns2) 
               endif
               drsdh =0.0    ! temporarily set drsdh =0 for heat adjustment of soil       
          else   
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              
               Tsoill(1)=Tsoill(1)+(diff_s*(sftmp-Tsoill(1))/thkns1 &
        &              -temph2)/(thkns1+thkns2)
          endif
        !!!!!  phase change in top soil       
          heat_excess=drsdh*(thd_t-Tsoill(i))+shcap(i)*thksl(i)*(Tsoill(i)-thd_t)/360000.         
          ice_incr=heat_excess*3600./latent_heat_fusion/ice_density        
          inter_var = ice(i)   
          if (ice_incr .lt. 0.) then     ! freeze             
              ice(i)=Amin1(liq_water(i)+inter_var,ice(i)-ice_incr)            
          else 
              ice(i) = Amax1(ice(i)-ice_incr,0.0)              
          endif
         !! readjust energy and temp 
          heat_adjust=heat_excess-latent_heat_fusion*(inter_var-ice(i))*ice_density/3600.
          Tsoill(i)=thd_t+heat_adjust/(shcap(i)*thksl(i)/360000.-drsdh)      
        else
            if ( i .gt. 9) then 
                temph2=0
                thkns2=500  ! boundary conditions, rethink
            endif
            Tsoill(i)=Tsoill(i)+(temph1-temph2)/(thkns1+thkns2)    
            heat_excess=shcap(i)*thksl(i)*(Tsoill(i)-thd_t)/360000.        
            ice_incr=heat_excess*3600./latent_heat_fusion/ice_density         
            inter_var = ice(i) 
            if (ice_incr .lt. 0.) then     ! freeze             
               ice(i)=Amin1(liq_water(i)+inter_var,ice(i)-ice_incr)             
            else 
               ice(i) = Amax1(ice(i)-ice_incr,0.0)              
            endif         
               !! readjust energy and temp 
            heat_adjust=heat_excess-latent_heat_fusion*(inter_var-ice(i))*ice_density/3600.
            Tsoill(i)=thd_t+heat_adjust/(shcap(i)/360000.*thksl(i))
        endif
       
        if (ABS(tsoill_pre -tsoill(i)) .gt. 5. ) Tsoill(i)=tsoill_pre
          TEMPH1=TEMPH2
          THKNS1=THKNS2
          DIFSV1=DIFSV2        
     enddo
     testout(1)=tsoill_0
     !testout(1)=tsnow
     testout(2:11)=tsoill(1:10)
     !testout(1:10)=tsoill(2)
     !testout(1)=dcount_soil
     !testout(1:10)=liq_water(1:10)
     !testout(1:10)=ice
     !testout(1)=Hsoil
     
     deallocate(depth_z)
     return 
     end 

 subroutine methane (Rh_pools,Tsoil,zwt,wsc,      &      !update single value in a hourly loop when MEMCMC=0
    &             phi,LAIMIN,LAIMAX,           &
    &             ProCH4,Pro_sum,OxiCH4,Oxi_sum,Fdifu,Ebu_sum,Pla_sum,simuCH4,CH4,   &
    &             r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi, &
    &             testout, do_soilphy)  !update single value of Rh_pools,Tsoil,zwt,wsc 
    !                                                           in a hourly loop when MEMCMC=1
!******************************************************************************************************************************
!****!introduce variables and constants used in this subroutine
!****************************************************************************************************************************** 
!     set soil layers
!****************************************************************************************************************************** 
      implicit none
      integer i,MEMCMC
      integer,parameter :: nlayers=10       !use this statement to set the parameter value
      real zwt    
      real consum
!****************************************************************************************************************************** 
!     set values for MEMCMC
!******************************************************************************************************************************       
      integer,parameter :: miterms=17
      integer,parameter :: ilines=9000      
!******************************************************************************************************************************
!     CH4 Production      
!******************************************************************************************************************************      
!      real Rhetero
      real Rh(nlayers),Rh_pools(5),Rh_h,ProCH4(nlayers),Pro_sum
      real r_me         !release ratio of CH4 to CO2
      real Q10pro
      real fSTP         !CH4 production factor of soil temperature
      real vt,xt
      real Tmax_me,Tpro_me
      real fpH          !CH4 production factor of soil pH
      real fEhP         !CH4 production factor of soil redox potential
      real FRLEN(nlayers)        !fraction of root in each layer
      real Tsoil
!****************************************************************************************************************************** 
!     CH4 Oxidation
!******************************************************************************************************************************      
      real CH4(nlayers)          !CH4 concentration
      real wsc(10)      
      real OxiCH4(nlayers),Oxi_sum       !CH4 oxidation
      real Omax_layers(nlayers),Omax       !maximum oxidation rate
      real kCH4_layers(nlayers),kCH4       !system specific half saturation constant
      real Q10oxi
      real fCH4(nlayers)         !CH4 oxidation factor of CH4 concentration
      real fSTO         !CH4 oxidation factor of soil temperature
      real fEhO         !CH4 oxidation factor of soil redox potential
      real Toxi
!****************************************************************************************************************************** 
!     CH4 Diffusion
!******************************************************************************************************************************      
      real Deff(nlayers)     !CH4 effective diffusivity  
      real D_CH4_a           !CH4 diffusion coefficient in air  !unit cm2 s-1   diffusivity of CH4 in air
      real D_CH4_w           !CH4 diffusion coefficient in water  !unit cm2 s-1   diffusivity of CH4 in water
      real phi          !soil porosity  also used in water table part
      real fwater,fair
      real D_CH4_soil(nlayers)
      real fcoarse      !relative volume of coarse pores depending on soil texture  Zhuang 2004
      real ftort        !tortuosity coefficient with a value of 0.66    Walter and Heimann 2000
      !suggesting that the distance covered by diffusion is about two thirds of the length of the real average path
      real SAND         !relative contents of sand (%) in the soil
      real PVSAND       !relative volume of coarse pores in sandy soils     set to 0.45     value from Walter 2001
      real SILT         !relative contents of silt (%) in the soil
      real PVSILT       !relative volume of coarse pores in silty soils     set to 0.20     value from Walter 2001
      real CLAY         !relative contents of clay (%) in the soil
      real PVCLAY       !relative volume of coarse pores in clayish soils     set to 0.14   value from Walter 2001
      real DEPTH(10)        !depth in soil  will define it inside this subroutine again      resolution 100mm 200mm
      real THKSL(10)        !will define it inside this subroutine again  
      real Fdifu(nlayers)
      real CH4_atm      !concentration of CH4 in atmosphere     seen as 0 cause the value is too low someone use 0.076
      real simuCH4      !simulated CH4 emission

!****************************************************************************************************************************** 
!     Ebullition 
!******************************************************************************************************************************      
      real CH4_thre,CH4_thre_layers(nlayers),EbuCH4(nlayers),Kebu
      real Ebu_sum_unsat,Ebu_sum_sat,Ebu_sum          !sum value one dimension is enough 
      integer wtlevelindex
!****************************************************************************************************************************** 
!     Plant transport
!******************************************************************************************************************************      
      real PlaCH4(nlayers),Pla_sum
      real LAIMIN,LAIMAX
      real Tveg,Tgr,Tmat,fgrow,Pox,Kpla
      ! Yuan added for soil temp  
      logical do_soilphy
      real testout(11), tsoil_layer(11)
!****************************************************************************************************************************** 
!******************************************************************************************************************************
      MEMCMC=1   ! note here, any changes here result unexpected bug 
      Rh_h=Rh_pools(1)+Rh_pools(2)+Rh_pools(3)+Rh_pools(4)+Rh_pools(5)  !hourly Rh_f + Rh_c + Rh_Micr + Rh_Slow + Rh_Pass
      
      tsoil_layer = testout
!      if (zwt .le. -199.99) then
!          zwt = -101.0
!      endif
!      Tsoil = 10.
!      wsc(1)=59.14
!      wsc(2)=94.9
!      wsc(3)=94.9
!      wsc(4)=94.9
!      wsc(5)=94.9
!      wsc(6)=189.8
!      wsc(7)=189.8
!      wsc(8)=189.8
!      wsc(9)=189.8
!      wsc(10)=189.8                    these are used for debugging
            
      data FRLEN /0.1,0.25,0.25,0.2,0.1,0.05,0.025,0.015,0.005,0.005/
      data thksl /20.,20.,20.,20.,20.,20.,20.,20.,20.,20./
    
      do i = 1, nlayers
                                        !!!!!!!put it out of the subroutine
         !****************************************************
         !* Rh weighed according to the distribution of root *
         !****************************************************
              if (i .LE. 3) then                                 ! the empirical method used here is from CLM4.5
                 Rh(i)= 0.5*Rh_h*FRLEN(i)+((0.5*Rh_h)/0.3)*0.1   
                 ! Rh(h,i)Rh produced by each layer per hour  unit should be g C m-2 h-1 
              else                                               ! i*10: depth of ith soil layers
                 Rh(i)= 0.5*Rh_h*FRLEN(i)
              endif
!              Rh(i) = Rh(i) + OxiCH4(i)*(11/4)
      enddo   
       
!****************************************************************************************************************************** 
!******************************************************************************************************************************            

     !****************************************************          
     !A. methane production     hourly  gC m-2 hour-1
     !Methane production is modeled as an anaerobic process that occurs in the saturated zone of the soil profile ZHUANG
     !****************************************************
     !Rh_h=Rh_pools(1)+Rh_pools(2)+Rh_pools(3)+Rh_pools(4)+Rh_pools(5)  !hourly Rh_f + Rh_c + Rh_Micr + Rh_Slow + Rh_Pass
     !r assignment
!      r_me=0.3      !find in parafile
      Tmax_me=45.0
!      Tpro_me=10.0
!      Q10pro=3.0    !find in parafile
          
      if (do_soilphy) then           
          if (Tsoil .lt. 0.0) then
              fSTP = 0.0
          else if (Tsoil .gt. Tmax_me) then
              fSTP = 0.0
          else if (Tsoil .ge. 0.0 .and. Tsoil .le. Tmax_me) then
              fSTP = Q10pro**((Tsoil-Tpro_me)/10)        !Tsoil is the only variable
          endif
      else 
          if (Tsoil .lt. 0.0) then
              fSTP = 0.0
          else if (Tsoil .gt. Tmax_me) then
              fSTP = 0.0
          else if (Tsoil .ge. 0.0 .and. Tsoil .le. Tmax_me) then
              fSTP = Q10pro**((Tsoil-Tpro_me)/10)        !Tsoil is the only variable
          endif
      endif
      
     !fpH assignment
      fpH=1.0
     !fEhP assignment
      fEhP=1.0
    
      depth(1)=10.0                                  !calculate soil depth unit cm
      do i=2,nlayers
          depth(i)=depth(i-1)+THKSL(i)
      enddo
      
      Pro_sum=0.0
      do i = 1,nlayers
!          (depth(i)*10)                   !convert unit from cm to mm
!          (THKSL(i)*10)                   !convert unit from cm to mm convert the unit in each of the equations
          if ((depth(i)*10) .le. -zwt) then
                  ProCH4(i)=0.0
          else
              if (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then
                  ProCH4(i)=Rh(i)*r_me*fSTP*fpH*fEhP*(((depth(i)*10.0)-(-zwt))/(THKSL(i)*10.0))     ! *percent
              elseif (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then
                  ProCH4(i)=Rh(i)*r_me*fSTP*fpH*fEhP
              endif
          endif
          Pro_sum=Pro_sum+ProCH4(i)
      enddo

     !**************************************************
     !Add CH4 production to CH4 pool    (gC layer -1)=(gC m-2)
     !**************************************************

      do i=1,nlayers
          CH4(i) = CH4(i) + ProCH4(i)
      enddo
!      CH4=0.0
!      ProCH4=0.0
!     END OF METHANE PRODUCTION
      
!     ********************************************************************************************************************
!     B. methane oxidation      hourly  unit gC m-2 h-1     !!!!!!!!!!!method of CLM and Zhuang!!!!!!!!!!!!
!     Methane oxidation is modeled as an aerobic process that occurs in the unsaturated zone of the soil profile ZHUANG
!     ********************************************************************************************************************
!     fSTO assignment
!     ***************          
      Q10oxi=2.0      !Zhu 2014 results from previous studies  unit 1  also used by zhang
!      Toxi=10.0       !Zhuang 2004 table1 Boreal Forest Wetland
      fSTO=Q10oxi**((Tsoil-Toxi)/10.0)

!     fEhO assignment
      fEhO=1.0        !Walter 2000  did not consider it, equal to value of 1

!     Omax assignment
!     ***************
      Oxi_sum=0.0
      do i = 1,nlayers
!      Omax=1.5
!      Omax=15.0  !!find in parafile Zhuang 2004 table1 Boreal Forest Wetland µmol L-1 h-1 system specific maximum oxidation coefficient
!     convert the unit of Omax from µmol L-1 h-1 to gC m-2 h-1
     !/1000,000 to get mol
     !*12 cmass to get gC
     !*1000 to get from dm-3(L) to m-3
     !*(wsc*0.001) to get unit of omax_layers from m-3 to m-2     !caution that wsc unit is mm
     !** w  /   (w/t)           CLM used 
      Omax_layers(i)=(Omax/(1000000))*12*1000*(wsc(i)*0.001)     !convert the unit of Omax from µmol L-1 h-1 to gC m-2 h-1
     !in unsaturated part of oxidation in CLM, they used the Omax/10 but did not expained why   they also /water volume
      
     !fCH4 assignment
!     ***************          
!      kCH4=5.0     !!find in parafile Zhuang 2004 range between 1 and 66.2µmol L-1 system specific half saturation constant  1.0e
!     convert the unit of kCH4 from µmol L-1 to gC m-2
      kCH4_layers(i)=(kCH4/(1000000))*12*1000*(wsc(i)*0.001)    !convert the unit of kCH4 from µmol L-1 to gC m-2
!     then calculate fCH4 with CH4(i) and kCH4_layers(i) 
      fCH4(i)=CH4(i)/(kCH4_layers(i)+CH4(i))   !??????????????  CH4 concentration factor

          if ((depth(i)*10.0) .le. -zwt) then                !unit of Omax: gC m-2 h-1
                  OxiCH4(i)=Omax_layers(i)*fCH4(i)*fSTO*fEhO!*0.1      !wrong:*(THKSL(i)/1000)!mm to m account for the thickness
!                  OxiCH4(i)=CH4(i)*0.001
          else
              if (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then
                  if (i .eq. 1) then
                  OxiCH4(i)=Omax_layers(i)*fCH4(i)*fSTO*fEhO*((-zwt)/(THKSL(i)*10.0))
                  else
                  OxiCH4(i)=Omax_layers(i)*fCH4(i)*fSTO*fEhO*(((-zwt)-(depth(i-1)*10.0))/(THKSL(i)*10.0))      !  *percent
                  endif
!                  OxiCH4(i)=CH4(i)*0.001*(((-zwt)-(depth(i-1)*10.0))/(THKSL(i)*10.0))
              else if (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then
                  OxiCH4(i)= 0.0
              endif
          endif
          
                    
          if (OxiCH4(i) .gt. CH4(i)) then
              OxiCH4(i)=CH4(i)
          endif
          
      Oxi_sum=Oxi_sum+OxiCH4(i)  
      
      enddo 

     !*******************************************************************
     !minus CH4 oxidation from CH4 pool     
     !*******************************************************************
      do i=1,nlayers
          CH4(i) = CH4(i) - OxiCH4(i)               !minus CH4 oxidation from CH4 pool
      enddo
!      CH4=0.0
!      OxiCH4=0.0
!     END OF METHANE OXIDATION
      
!     ****************************************************
!     C. methane diffusion
!     ****************************************************
!     Parameters assignment 
      D_CH4_a=0.2            !unit cm2 s-1   D_CH4_a is the molecular diffusion coefficient of methane in air
      D_CH4_a=(D_CH4_a/10000.0)*3600.0        !unit m2 h-1
     
      D_CH4_w=0.00002        !unit cm2 s-1   D_CH4_a is the molecular diffusion coefficient of methane in water
      D_CH4_w=(D_CH4_w/10000.0)*3600.0        !unit m2 h-1          
!      write(*,*)'D_ch4', D_CH4_w
      ftort=0.66        !tortuosity coefficient with a value of 0.66    Walter and Heimann 2000
!      phi=0.95          !might need to be divided into different layers
!     WFPS=         !water-filled pore space      
!     parameters for fcoarse algorithm      
      SAND=0.4             !   %   SPRUCE site value
      SILT=0.4             !   %   SPRUCE site value
      CLAY=0.2             !   %   SPRUCE site value
      PVSAND=0.45       !relative volume of coarse pores in sandy soils       set to 0.45     value from Walter 2001 zhuang
      PVSILT=0.20       !relative volume of coarse pores in silty soils       set to 0.20     value from Walter 2001 zhuang
      PVCLAY=0.14       !relative volume of coarse pores in clayish soils     set to 0.14     value from Walter 2001 zhuang  
      fcoarse=SAND*PVSAND+SILT*PVSILT+CLAY*PVCLAY
      CH4_atm=0.0
      
!!      ***************************************************
!!      * Old way for Deff: abandoned due to oscillation  *
!!      ***************************************************      
!!      do i = 1,nlayers
!!          if ((depth(i)*10.0) .le. -zwt) then
!!                Deff(i)=ftort*D_CH4_a*fcoarse       !unit m2 h-1
!!          else
!!              if (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then
!!                  Deff(i)=ftort*D_CH4_a*fcoarse*(((-zwt)-(depth(i-1)*10.0))/(THKSL(i)*10.0))    &
!!                  &     +ftort*D_CH4_w*fcoarse*(((depth(i)*10.0)-(-zwt))/(THKSL(i)*10.0)) 
!!              else if (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then
!!                  Deff(i)=ftort*D_CH4_w*fcoarse     !unit m2 h-1
!!              endif
!!          endif
!!      enddo      

      
!       ***************************************************
!       * New solution for Deff: Equations for D_CH4_soil *
!       ***************************************************
!!!        phi = 0.95;
!!!        zthetasmin = 0.1*1000;  %m to mm
!!!        thetasmin = 0.25;
!!!        depth = 100; % mm depth i-1
!!!        % wsc = 0:0.01:94.9;
!!!        zwt = 101:1:200;  % mm
!!!        m = length(zwt);
!!!        az = (phi-thetasmin)/zthetasmin;            %constant
!!!          thetas = max(thetasmin,phi-(az*zwt*0.01))                  !volumetric water content in mires quoted by Granberg
!!!          z = (zwt-depth)/2+depth                                    !use the middle position to represent depth for whole layer 
!!!        thetaus = min(phi,(thetas + (phi-thetas)*(z/zwt)**2))        !volumetric water content in the unsaturated zone
!!!          fwater = thetaus/(zwt-depth)
!!!          fair = phi-fwater
!!
!!!        D_CH4_soil_1 = (((fair)**(10/3))/((phi)**2))*D_CH4_a
!!!        D_CH4_soil_2 = D_CH4_W
!!!        
!!!        Dsoil = ge(fair,0.05).*Dsoil_1 + lt(fair,0.05)*Dsoil_2

!      do i = 1,nlayers
!          fwater = wsc(i)/(THKSL(i)*10)
!          fair = phi-fwater
!          D_CH4_soil(i)=(((fair)**(10/3))/((phi)**2))*D_CH4_a
!          
!          if ((depth(i)*10.0) .le. -zwt) then
!              Deff(i)=ftort*D_CH4_soil(i)*fcoarse                          !unit m2 h-1
!          else
!              if (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then
!                  if (fair .le. 0.05) then
!                      Deff(i)=D_CH4_w
!                  else
!                      Deff(i)=D_CH4_soil(i)
!                  endif
!              else if (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then
!                  Deff(i)=D_CH4_w                                                    !unit m2 h-1
!              endif
!          endif
!      enddo 
      
      
      
      do i = 1,nlayers
          fwater = wsc(i)/(THKSL(i)*10)
!          fair = phi-fwater
!          D_CH4_soil(i)=(((fair)**(10/3))/((phi)**2))*D_CH4_a
                    
          Deff(i) = D_CH4_a*fcoarse*ftort*phi*(phi-fwater)+D_CH4_w*fwater
                   
!          if ((depth(i)*10.0) .le. -zwt) then
!              Deff(i)=D_CH4_soil(i)                                                  !unit m2 h-1
!          else
!              if (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then
!                  if (fair .le. 0.05) then
!                      Deff(i)=D_CH4_w
!                  else
!                      Deff(i)=D_CH4_soil(i)
!                  endif
!              else if (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then
!                  Deff(i)=D_CH4_w                                                    !unit m2 h-1
!              endif
!          endif
      enddo 
!             
      Fdifu(1) = Deff(1)*(CH4(1)-CH4_atm)/(THKSL(1)*0.01)         !refer to flux from layer 1 to atmosphere  cm to m
      do i = 2,nlayers                                  !refer to flux from layer ii to ii-1 
          Fdifu(i)= Deff(i)*(CH4(i)-CH4(i-1))/(THKSL(i)*0.01)      !###**** is THKSL(ii) correct?
      enddo                                             !need the value of CH4 concentration in atm    assume to be 0
      
      do i=1,nlayers
          if (Fdifu(i) .gt. 0.0 .and. (Fdifu(i))/(THKSL(i)*0.01) .gt. CH4(i)) then
              Fdifu(i)=CH4(i)*(THKSL(i)*0.01)
          else if (Fdifu(i) .lt. 0.0 .and. (abs(Fdifu(i))/(THKSL(i)*0.01)) .gt. CH4(i-1)) then
              Fdifu(i)=-CH4(i-1)*(THKSL(i-1)*0.01)    
          endif
      enddo
      
!      CH4(1) = CH4(1) + (0.0+Fdifu(1))/(THKSL(1)*0.01)
      do i = 1,nlayers-1                                  !loop of time
          CH4(i) = CH4(i) + (Fdifu(i+1)-Fdifu(i))/(THKSL(i)*0.01)       !###****PLUS OR MINUS??  
!          CH4(i) = CH4(i) - 0.1*CH4(i)
          if (CH4(i) .lt. 0.0) then
              CH4(i) = 0.0
          endif
      enddo    
        CH4(10) = CH4(10) - Fdifu(10)/(THKSL(10)*0.01)                  !MODIFIED ON 07/25/2016
          if (CH4(10) .lt. 0.0) then                                    !defined the Fdifu(11) to be 0.0
              CH4(10)= 0.0
          endif
    ! measured CH4flux is Fdifu(1)
!      simuCH4 = simuCH4+Fdifu(1)
      simuCH4 = (Fdifu(1)-0.0)/(THKSL(1)*0.01)            !starting to calculate the CH4 emission       0.0 atm
!      CH4=0.0
!      Fdifu=0.0 

!     ********************************************************************************************************************      
    ! D. methane ebullition     !assume bubbles can reach the water table within 1 h&
                                !& the bubbles is added to the methane concentration in the soil layer just above the wt
                                !& and then diffused through layers
!     ********************************************************************************************************************
      Kebu=1.0                    !unit  h-1   rate constant
    ! here in other papers the depth resolution is 1 centimeter
    ! ##might need to be multiplied by 10??????
      
      do i=1,nlayers
!          CH4_thre=500.0  !!find in parafile  !unit  umol L-1 according to Walter's 500-1000
          CH4_thre_layers(i)=(CH4_thre*1.0e-6)*12*1000*(wsc(i)*0.001)    !convert the unit of CH4_thre from µmol L-1 to gC m-2
      enddo
    
      if (zwt .ge. 0.0) then                                  !when water table is above the soil surface
          Ebu_sum_unsat=0.0
          Ebu_sum_sat=0.0                                      !initial value
          do i=1,nlayers
              if (CH4(i) .gt. CH4_thre_layers(i)) then                  
                      EbuCH4(i)=Kebu*(CH4(i)-CH4_thre_layers(i))     !only if the concentration is larger than threshold
              else !if (CH4(i) .le. CH4_thre_layers(i)) then
                      EbuCH4(i)=0.0
              endif
              Ebu_sum_sat=Ebu_sum_sat+EbuCH4(i)               !& the bubbles are directly added into CH4 efflux into atmosphere
              CH4(i)=CH4(i)- EbuCH4(i)                        !& update the concentration at the end of this hour in each layers
          enddo
          Ebu_sum=Ebu_sum_unsat+Ebu_sum_sat
          simuCH4=simuCH4+Ebu_sum_sat                         !& the bubbles are directly added into CH4 efflux into atmosphere

      endif
      
      if (zwt .lt. 0.0) then                                  !when water table is below the soil surface
        Ebu_sum_sat=0.0
        Ebu_sum_unsat=0.0                                     !initial value
        do i=1,nlayers
            if ((depth(i)*10.0) .le. -zwt) then               !acrotelm layers
                EbuCH4(i)=0.0
                Ebu_sum_unsat=Ebu_sum_unsat+EbuCH4(i)         
                CH4(i)=CH4(i)- EbuCH4(i) 
            else
                if (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then       !partly acrotelm layer
                    wtlevelindex = i
                    if (CH4(i) .gt. CH4_thre_layers(i)) then                  
                      EbuCH4(i)=Kebu*(CH4(i)-CH4_thre_layers(i))!*(((depth(i)*10.0)-(-zwt))/(THKSL(i)*10.0))        ! * percent
                    else !if (CH4(i) .le. CH4_thre_layers(i)) then                     ??????????,??????????????????
                      EbuCH4(i)=0.0
                    endif 
                  CH4(i)=CH4(i)- EbuCH4(i)
                  Ebu_sum_unsat=Ebu_sum_unsat+EbuCH4(i)
                else if (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then   !catotelm layers
                    if (CH4(i) .gt. CH4_thre_layers(i)) then                  
                      EbuCH4(i)=Kebu*(CH4(i)-CH4_thre_layers(i))
                    else !if (CH4(i) .le. CH4_thre_layers(i)) then
                      EbuCH4(i)=0.0
                    endif 
                  CH4(i)=CH4(i)- EbuCH4(i)
                  Ebu_sum_unsat=Ebu_sum_unsat+EbuCH4(i)
                endif
            endif
        enddo
!        CH4(wtlevelindex)=CH4(wtlevelindex)+ Ebu_sum_unsat           ! assume that bubbles emission into air within 1 hour
        Ebu_sum=Ebu_sum_unsat+Ebu_sum_sat
        simuCH4=simuCH4+Ebu_sum_unsat                                 ! assume that bubbles emission into air within 1 hour
      endif
      !Ebu_sum_sat=0.0 
!      Ebu_sum_unsat=0.0
!      EbuCH4=0.0
            
!     ******************************************************************************************************
    ! E. plant mediated methane transportation      totoally used Walter's model also used by Zhuang et. al
!     ******************************************************************************************************
      Kpla=0.01         !unit h-1
!      Tveg=0.3 ! a factor describing the quality of plant-mediated transport depend on the density of plant stands and plant types 
!      0 for boreal forest and 0.5 for tundra
      !find in parafile !
    ! the Tsoil used here would be better if refer to the 20cm soil temperature after &
    ! & the accomplishment of soil heat dynamics module. according to Zhuang. however Walter used 50cm soil temp.
      Tgr=2.0               !unit degree Celsius if annual mean temp is below 5 (otherwise 7)
      Tmat=Tgr+10.0         !unit degree Celsius
      Pox=0.5               !50% of mediated methane are oxidised 
    ! define fgrow
      if (Tsoil .lt. Tgr) then
          fgrow=LAIMIN
      else if (Tsoil .ge. Tgr .and. Tsoil .le. Tmat) then
          fgrow=LAIMIN+LAIMAX*(1-((Tmat-Tsoil)/(Tmat-Tgr))**2)
      else if (Tsoil .gt. Tmat) then
          fgrow=LAIMAX
      endif
      
      Pla_sum=0.0
      do i=1,nlayers
          PlaCH4(i)=Kpla*Tveg*FRLEN(i)*fgrow*CH4(i)*(1-Pox)         !not sensitive at all to this change, but better
!          PlaCH4(i)=Kpla*Tveg*FRLEN(i)*fgrow*CH4(i)
          Pla_sum=Pla_sum+PlaCH4(i)
          CH4(i)=CH4(i)-PlaCH4(i)
      enddo
      
      simuCH4=simuCH4+Pla_sum
      consum=simuCH4+OxiCH4(1)+OxiCH4(2)+OxiCH4(3)+OxiCH4(4)+OxiCH4(5)+OxiCH4(6)+OxiCH4(7)+OxiCH4(8)+OxiCH4(9)+OxiCH4(10)
      
      if (MEMCMC .eq. 0) then
      !write(*,*) 'zwt',zwt,'simuCH4',simuCH4         !show on screen          
      !write(82,182)zwt,Ebu_sum_sat,Ebu_sum_unsat,Pla_sum,simuCH4, &
      !        & Rh(1),Rh(2),Rh(3),Rh(4),Rh(5),Rh(6),Rh(7),Rh(8),Rh(9),Rh(10),   &
      !        & consum,Pro_sum, &
      !        & ProCH4(1),ProCH4(2),ProCH4(3),ProCH4(4),ProCH4(5),ProCH4(6),ProCH4(7),ProCH4(8),ProCH4(9),ProCH4(10),   &
      !        & CH4(1),CH4(2),CH4(3),CH4(4),CH4(5),CH4(6),CH4(7),CH4(8),CH4(9),CH4(10), &
      !        & Fdifu(1),Fdifu(2),Fdifu(3),Fdifu(4),Fdifu(5),Fdifu(6),Fdifu(7),Fdifu(8),Fdifu(9),Fdifu(10), &
      !        & OxiCH4(1),OxiCH4(2),OxiCH4(3),OxiCH4(4),OxiCH4(5),OxiCH4(6),OxiCH4(7),OxiCH4(8),OxiCH4(9),OxiCH4(10),   &
      !        & wsc(1),wsc(2),wsc(3),wsc(4),wsc(5),wsc(6),wsc(7),wsc(8),wsc(9),wsc(10), &
      !        & Rh_pools(1),Rh_pools(2),Rh_pools(3),Rh_pools(4),Tsoil

182   format(72(f15.9,","))  

      !write(83,183)zwt,Tsoil,Rh_pools(1),Rh_pools(2),Rh_pools(3),Rh_pools(4),Rh_pools(5), &
       !       & wsc(1),wsc(2),wsc(3),wsc(4),wsc(5),wsc(6),wsc(7),wsc(8),wsc(9),wsc(10)

183   format(17(f15.9,","))      
      endif    
      
      return
      end    
       
!     ========================================================================================
!     subroutines used by canopy submodel
      subroutine xlayers(Sps,Tair,Dair,radabv,fbeam,eairP,&
     &           wind,co2ca,fwsoil,wcl,FLAIT,coszen,idoy,hours,&
     &           tauL,rhoL,rhoS,xfang,extkd,extkU,wleaf,&
     &           Rconst,sigma,emleaf,emsoil,theta,a1,Ds0,&
     &           cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,&
     &           gsw0,alpha,stom_n,wsmax,wsmin,&
     &           Vcmx0,eJmx0,conKc0,conKo0,Ekc,Eko,o2ci,&
     &           Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,&
     &           extKb,Rsoilabs,Acan1,Acan2,Ecan1,Ecan2,&
     &           RnStL,QcanL,RcanL,AcanL,EcanL,HcanL,GbwcL,GswcL,gddonset,testout,&
     &           Rsoilab1,Rsoilab2,QLleaf,QLair,raero,do_soilphy,G,Esoil,Hsoil)


     
!    the multi-layered canopy model developed by 
!    Ray Leuning with the new radiative transfer scheme   
!    implemented by Y.P. Wang (from Sellers 1986)
!    12/Sept/96 (YPW) correction for mean surface temperature of sunlit
!    and shaded leaves
!    Tleaf,i=sum{Tleaf,i(n)*fslt*Gaussw(n)}/sum{fslt*Gaussw(n)} 
!    
      real Gaussx(5),Gaussw(5)
      real layer1(5),layer2(5)
      real tauL(3),rhoL(3),rhoS(3),Qabs(3,2),Radabv(2),Rnstar(2)
      real Aleaf(2),Eleaf(2),Hleaf(2),Tleaf(2),co2ci(2)
      real gbleaf(2),gsleaf(2),QSabs(3,2),Qasoil(2)
      integer ng,nw
      real rhoc(3,2),reff(3,2),kpr(3,2),scatt(2)       !Goudriaan
      real testout(11)

      real rsoil,rlai,raero
      real wsmax,wsmin,WILTPT,FILDCP,wcl(10)
      real gddonset
!    additional arrays to allow output of info for each Layer
      real RnStL(5),QcanL(5),RcanL(5),AcanL(5),EcanL(5),HcanL(5)
      real GbwcL(5),GswcL(5)
      logical do_soilphy
      
! Normalised Gaussian points and weights (Goudriaan & van Laar, 1993, P98)
!* 5-point
      data Gaussx/0.0469101,0.2307534,0.5,0.7692465,0.9530899/
      data Gaussw/0.1184635,0.2393144,0.2844444,0.2393144,0.1184635/

!     soil water conditions
      WILTPT=wsmin/100.
      FILDCP=wsmax/100.
!     reset the vairables
      Rnst1=0.0        !net rad, sunlit
      Rnst2=0.0        !net rad, shaded
      Qcan1=0.0        !vis rad
      Qcan2=0.0
      Rcan1=0.0        !NIR rad
      Rcan2=0.0
      Acan1=0.0        !CO2
      Acan2=0.0
      Ecan1=0.0        !Evap
      Ecan2=0.0
      Hcan1=0.0        !Sens heat
      Hcan2=0.0
      Gbwc1=0.0        !Boundary layer conductance
      Gbwc2=0.0
      Gswc1=0.0        !Canopy conductance
      Gswc2=0.0
      Tleaf1=0.0       !Leaf Temp
      Tleaf2=0.0  
  
!     aerodynamic resistance                                                
      raero=50./wind                           

!    Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
      xphi1 = 0.5 - 0.633*xfang -0.33*xfang*xfang
      xphi2 = 0.877 * (1.0 - 2.0*xphi1)
      funG=xphi1 + xphi2*coszen                             !G-function: Projection of unit leaf area in direction of beam
      
      if(coszen.gt.0) then                                  !check if day or night
        extKb=funG/coszen                                   !beam extinction coeff - black leaves
      else
        extKb=100.
      end if

!     Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
!     Effective extinction coefficient for diffuse radiation Goudriaan & van Laar Eq 6.6)
      pi180=3.1416/180.
      cozen15=cos(pi180*15)
      cozen45=cos(pi180*45)
      cozen75=cos(pi180*75)
      xK15=xphi1/cozen15+xphi2
      xK45=xphi1/cozen45+xphi2
      xK75=xphi1/cozen75+xphi2
      transd=0.308*exp(-xK15*FLAIT)+0.514*exp(-xK45*FLAIT)+     &
     &       0.178*exp(-xK75*FLAIT)
      extkd=(-1./FLAIT)*alog(transd)
      extkn=extkd                        !N distribution coeff 

!canopy reflection coefficients (Array indices: first;  1=VIS,  2=NIR
!                                               second; 1=beam, 2=diffuse
      do nw=1,2                                                      !nw:1=VIS, 2=NIR
       scatt(nw)=tauL(nw)+rhoL(nw)                      !scattering coeff
       if((1.-scatt(nw))<0.0)scatt(nw)=0.9999           ! Weng 10/31/2008
       kpr(nw,1)=extKb*sqrt(1.-scatt(nw))               !modified k beam scattered (6.20)
       kpr(nw,2)=extkd*sqrt(1.-scatt(nw))             !modified k diffuse (6.20)
       rhoch=(1.-sqrt(1.-scatt(nw)))/(1.+sqrt(1.-scatt(nw)))            !canopy reflection black horizontal leaves (6.19)
       rhoc15=2.*xK15*rhoch/(xK15+extkd)                                !canopy reflection (6.21) diffuse
       rhoc45=2.*xK45*rhoch/(xK45+extkd)
       rhoc75=2.*xK75*rhoch/(xK75+extkd)
       rhoc(nw,2)=0.308*rhoc15+0.514*rhoc45+0.178*rhoc75
       rhoc(nw,1)=2.*extKb/(extKb+extkd)*rhoch                          !canopy reflection (6.21) beam 
       reff(nw,1)=rhoc(nw,1)+(rhoS(nw)-rhoc(nw,1))   &                   !effective canopy-soil reflection coeff - beam (6.27)
     &            *exp(-2.*kpr(nw,1)*FLAIT) 
       reff(nw,2)=rhoc(nw,2)+(rhoS(nw)-rhoc(nw,2))   &                   !effective canopy-soil reflection coeff - diffuse (6.27)
     &            *exp(-2.*kpr(nw,2)*FLAIT)  
      enddo


!     isothermal net radiation & radiation conductance at canopy top - needed to calc emair
      call Radiso(flait,flait,Qabs,extkd,Tair,eairP,cpair,Patm, &
     &            fbeam,airMa,Rconst,sigma,emleaf,emsoil,       &
     &            emair,Rnstar,grdn)
      TairK=Tair+273.2

!     below      
      do ng=1,5
         flai=gaussx(ng)*FLAIT
!        radiation absorption for visible and near infra-red
         call goudriaan(FLAI,coszen,radabv,fbeam,reff,kpr,      &
     &                  scatt,xfang,Qabs) 
!        isothermal net radiation & radiation conductance at canopy top
         call Radiso(flai,flait,Qabs,extkd,Tair,eairP,cpair,Patm,   &
     &               fbeam,airMa,Rconst,sigma,emleaf,emsoil,        &
     &               emair,Rnstar,grdn)
         windUx=wind*exp(-extkU*flai)             !windspeed at depth xi
         scalex=exp(-extkn*flai)                    !scale Vcmx0 & Jmax0
         Vcmxx=Vcmx0*scalex
         eJmxx=eJmx0*scalex
         if(radabv(1).ge.10.0) then                          !check solar Radiation > 10 W/m2
!           leaf stomata-photosynthesis-transpiration model - daytime
            call agsean_day(Sps,Qabs,Rnstar,grdn,windUx,Tair,Dair,      &
     &               co2ca,wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,  &
     &               Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,  &
     &               gsw0,alpha,stom_n,                                 &
     &               Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,            &
     &               Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,         &
     &               Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci,gddonset)       
         else
            call agsean_ngt(Sps,Qabs,Rnstar,grdn,windUx,Tair,Dair,      &
     &               co2ca,wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,  &
     &               Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,  &
     &               gsw0,alpha,stom_n,                                 &
     &               Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,            &
     &               Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,         &
     &               Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci)
         endif  
         fslt=exp(-extKb*flai)                        !fraction of sunlit leaves
         fshd=1.0-fslt                                !fraction of shaded leaves
         Rnst1=Rnst1+fslt*Rnstar(1)*Gaussw(ng)*FLAIT  !Isothermal net rad`
         Rnst2=Rnst2+fshd*Rnstar(2)*Gaussw(ng)*FLAIT
         RnstL(ng)=Rnst1+Rnst2
!
         Qcan1=Qcan1+fslt*Qabs(1,1)*Gaussw(ng)*FLAIT  !visible
         Qcan2=Qcan2+fshd*Qabs(1,2)*Gaussw(ng)*FLAIT
         QcanL(ng)=Qcan1+Qcan2
!
         Rcan1=Rcan1+fslt*Qabs(2,1)*Gaussw(ng)*FLAIT  !NIR
         Rcan2=Rcan2+fshd*Qabs(2,2)*Gaussw(ng)*FLAIT
         RcanL(ng)=Rcan1+Rcan2
!
         if(Aleaf(1).lt.0.0)Aleaf(1)=0.0      !Weng 2/16/2006
         if(Aleaf(2).lt.0.0)Aleaf(2)=0.0      !Weng 2/16/2006

         Acan1=Acan1+fslt*Aleaf(1)*Gaussw(ng)*FLAIT*stom_n    !amphi/hypostomatous
         Acan2=Acan2+fshd*Aleaf(2)*Gaussw(ng)*FLAIT*stom_n
         AcanL(ng)=Acan1+Acan2

         layer1(ng)=Aleaf(1)
         layer2(ng)=Aleaf(2)

         Ecan1=Ecan1+fslt*Eleaf(1)*Gaussw(ng)*FLAIT
         Ecan2=Ecan2+fshd*Eleaf(2)*Gaussw(ng)*FLAIT
         EcanL(ng)=Ecan1+Ecan2
!
         Hcan1=Hcan1+fslt*Hleaf(1)*Gaussw(ng)*FLAIT
         Hcan2=Hcan2+fshd*Hleaf(2)*Gaussw(ng)*FLAIT
         HcanL(ng)=Hcan1+Hcan2
!
         Gbwc1=Gbwc1+fslt*gbleaf(1)*Gaussw(ng)*FLAIT*stom_n
         Gbwc2=Gbwc2+fshd*gbleaf(2)*Gaussw(ng)*FLAIT*stom_n
!
         Gswc1=Gswc1+fslt*gsleaf(1)*Gaussw(ng)*FLAIT*stom_n
         Gswc2=Gswc2+fshd*gsleaf(2)*Gaussw(ng)*FLAIT*stom_n
!
         Tleaf1=Tleaf1+fslt*Tleaf(1)*Gaussw(ng)*FLAIT
         Tleaf2=Tleaf2+fshd*Tleaf(2)*Gaussw(ng)*FLAIT
      enddo  ! 5 layers

      FLAIT1=(1.0-exp(-extKb*FLAIT))/extkb
      Tleaf1=Tleaf1/FLAIT1
      Tleaf2=Tleaf2/(FLAIT-FLAIT1)

!     Soil surface energy and water fluxes
!    Radiation absorbed by soil
      Rsoilab1=fbeam*(1.-reff(1,1))*exp(-kpr(1,1)*FLAIT)        &
     &         +(1.-fbeam)*(1.-reff(1,2))*exp(-kpr(1,2)*FLAIT)          !visible
      Rsoilab2=fbeam*(1.-reff(2,1))*exp(-kpr(2,1)*FLAIT)        &
     &         +(1.-fbeam)*(1.-reff(2,2))*exp(-kpr(2,2)*FLAIT)          !NIR
      Rsoilab1=Rsoilab1*Radabv(1)
      Rsoilab2=Rsoilab2*Radabv(2)
!  
      Tlk1=Tleaf1+273.2
      Tlk2=Tleaf2+273.2
!      temp1=-extkd*FLAIT
      QLair=emair*sigma*(TairK**4)*exp(-extkd*FLAIT)
      QLleaf=emleaf*sigma*(Tlk1**4)*exp(-extkb*FLAIT)           &
     &      +emleaf*sigma*(Tlk2**4)*(1.0-exp(-extkb*FLAIT))
      QLleaf=QLleaf*(1.0-exp(-extkd*FLAIT)) 
      QLsoil=emsoil*sigma*(TairK**4)
      Rsoilab3=(QLair+QLleaf)*(1.0-rhoS(3))-QLsoil

!    Net radiation absorbed by soil
!    the old version of net long-wave radiation absorbed by soils 
!    (with isothermal assumption)
!     Rsoil3=(sigma*TairK**4)*(emair-emleaf)*exp(-extkd*FLAIT)         !Longwave
!     Rsoilab3=(1-rhoS(3))*Rsoil3

!    Total radiation absorbed by soil    
      Rsoilabs=Rsoilab1+Rsoilab2+Rsoilab3 

!    thermodynamic parameters for air
      TairK=Tair+273.2
      rhocp=cpair*Patm*AirMa/(Rconst*TairK)
      H2OLv=H2oLv0-2.365e3*Tair
      slope=(esat(Tair+0.1)-esat(Tair))/0.1
      psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar=Patm/(Rconst*TairK)
      fw1=AMIN1(AMAX1((FILDCP-wcl(1))/(FILDCP-WILTPT),0.05),1.0)
      Rsoil=30.*exp(0.2/fw1)
      rLAI=exp(FLAIT)
!     latent heat flux into air from soil
!           Eleaf(ileaf)=1.0*
!     &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    !2* Weng 0215
!     &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
      
      if (do_soilphy) then
          Esoil=(slope*(Rsoilabs-G)+rhocp*Dair/(raero+rLAI))/       &
          &      (slope+psyc*(rsoil/(raero+rLAI)+1.))
!     sensible heat flux into air from soil
          Hsoil=Rsoilabs-Esoil-G
      endif 
      !testout(1)=Rsoilabs
      return
      end 

!     ****************************************************************************
      subroutine goudriaan(FLAI,coszen,radabv,fbeam,reff,kpr,   &
     &                  scatt,xfang,Qabs)
     
!    for spheric leaf angle distribution only
!    compute within canopy radiation (PAR and near infra-red bands)
!    using two-stream approximation (Goudriaan & vanLaar 1994)
!    tauL: leaf transmittance
!    rhoL: leaf reflectance
!    rhoS: soil reflectance
!    sfang XiL function of Ross (1975) - allows for departure from spherical LAD
!         (-1 vertical, +1 horizontal leaves, 0 spherical)
!    FLAI: canopy leaf area index
!    funG: Ross' G function
!    scatB: upscatter parameter for direct beam
!    scatD: upscatter parameter for diffuse
!    albedo: single scattering albedo
!    output:
!    Qabs(nwave,type), nwave=1 for visible; =2 for NIR,
!                       type=1 for sunlit;   =2 for shaded (W/m2)

      real radabv(2)
      real Qabs(3,2),reff(3,2),kpr(3,2),scatt(2)
      xu=coszen                                         !cos zenith angle
      
!     Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
      xphi1 = 0.5 - 0.633*xfang -0.33*xfang*xfang
      xphi2 = 0.877 * (1.0 - 2.0*xphi1)
      funG=xphi1 + xphi2*xu                             !G-function: Projection of unit leaf area in direction of beam
      
      if(coszen.gt.0) then                                  !check if day or night
        extKb=funG/coszen                                   !beam extinction coeff - black leaves
      else
        extKb=100.
      end if
                       
! Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
      do nw=1,2
       Qd0=(1.-fbeam)*radabv(nw)                                          !diffuse incident radiation
       Qb0=fbeam*radabv(nw)                                               !beam incident radiation
       Qabs(nw,2)=Qd0*(kpr(nw,2)*(1.-reff(nw,2))*exp(-kpr(nw,2)*FLAI))+  & !absorbed radiation - shaded leaves, diffuse
     &            Qb0*(kpr(nw,1)*(1.-reff(nw,1))*exp(-kpr(nw,1)*FLAI)-   & !beam scattered
     &            extKb*(1.-scatt(nw))*exp(-extKb*FLAI))
       Qabs(nw,1)=Qabs(nw,2)+extKb*Qb0*(1.-scatt(nw))                     !absorbed radiation - sunlit leaves 
      end do
      return
      end

!****************************************************************************
      subroutine Radiso(flai,flait,Qabs,extkd,Tair,eairP,cpair,Patm,    &
     &                  fbeam,airMa,Rconst,sigma,emleaf,emsoil,         &
     &                  emair,Rnstar,grdn)
!     output
!     Rnstar(type): type=1 for sunlit; =2 for shaded leaves (W/m2)
!     23 Dec 1994
!     calculates isothermal net radiation for sunlit and shaded leaves under clear skies
!     implicit real (a-z)
      real Rnstar(2)
      real Qabs(3,2)
      TairK=Tair+273.2

! thermodynamic properties of air
      rhocp=cpair*Patm*airMa/(Rconst*TairK)   !volumetric heat capacity (J/m3/K)

! apparent atmospheric emissivity for clear skies (Brutsaert, 1975)
      emsky=0.642*(eairP/Tairk)**(1./7)       !note eair in Pa
     
! apparent emissivity from clouds (Kimball et al 1982)
      ep8z=0.24+2.98e-12*eairP*eairP*exp(3000/TairK)
      tau8=amin1(1.0,1.0-ep8z*(1.4-0.4*ep8z))            !ensure tau8<1
      emcloud=0.36*tau8*(1.-fbeam)*(1-10./TairK)**4      !10 from Tcloud = Tair-10

! apparent emissivity from sky plus clouds      
!      emair=emsky+emcloud
! 20/06/96
      emair=emsky

      if(emair.gt.1.0) emair=1.0
      
! net isothermal outgoing longwave radiation per unit leaf area at canopy
! top & thin layer at flai (Note Rn* = Sn + Bn is used rather than Rn* = Sn - Bn in Leuning et al 1985)
      Bn0=sigma*(TairK**4.)
      Bnxi=Bn0*extkd*(exp(-extkd*flai)*(emair-emleaf)       &
     &    + exp(-extkd*(flait-flai))*(emsoil-emleaf))
!     isothermal net radiation per unit leaf area for thin layer of sunlit and
!     shaded leaves
      Rnstar(1)=Qabs(1,1)+Qabs(2,1)+Bnxi
      Rnstar(2)=Qabs(1,2)+Qabs(2,2)+Bnxi
!     radiation conductance (m/s) @ flai
      grdn=4.*sigma*(TairK**3.)*extkd*emleaf*               &       ! corrected by Jiang Jiang 2015/9/29
     &    (exp(-extkd*flai)+exp(-extkd*(flait-flai)))       &
     &    /rhocp
      return
      end
!     ****************************************************************************
      subroutine agsean_day(Sps,Qabs,Rnstar,grdn,windUx,Tair,Dair,      &
     &               co2ca,wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,  &
     &               Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,  &
     &               gsw0,alpha,stom_n,                                 &
     &               Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,            &
     &               Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,         &
     &               Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci,gddonset)

!    implicit real (a-z)
      integer kr1,ileaf
      real Aleaf(2),Eleaf(2),Hleaf(2),Tleaf(2),co2ci(2)
      real gbleaf(2), gsleaf(2)
      real Qabs(3,2),Rnstar(2)
!    thermodynamic parameters for air
      TairK=Tair+273.2
      rhocp=cpair*Patm*AirMa/(Rconst*TairK)
      H2OLv=H2oLv0-2.365e3*Tair
      slope=(esat(Tair+0.1)-esat(Tair))/0.1
      psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar=Patm/(Rconst*TairK)
      weighJ=1.0
!    boundary layer conductance for heat - single sided, forced convection
!    (Monteith 1973, P106 & notes dated 23/12/94)
      if(windUx/wleaf>=0.0)then
          gbHu=0.003*sqrt(windUx/wleaf)    !m/s
      else
          gbHu=0.003 !*sqrt(-windUx/wleaf)
      endif         ! Weng 10/31/2008
!     raero=0.0                        !aerodynamic resistance s/m
      do ileaf=1,2              ! loop over sunlit and shaded leaves
!        first estimate of leaf temperature - assume air temp
         Tleaf(ileaf)=Tair
         Tlk=Tleaf(ileaf)+273.2    !Tleaf to deg K
!        first estimate of deficit at leaf surface - assume Da
         Dleaf=Dair                !Pa
!        first estimate for co2cs
         co2cs=co2ca               !mol/mol
         Qapar = (4.6e-6)*Qabs(1,ileaf)
!    ********************************************************************
         kr1=0                     !iteration counter for LE
!        return point for evaporation iteration
         do               !iteration for leaf temperature
!          single-sided boundary layer conductance - free convection (see notes 23/12/94)
           Gras=1.595e8*ABS(Tleaf(ileaf)-Tair)*(wleaf**3.)     !Grashof
           gbHf=0.5*Dheat*(Gras**0.25)/wleaf
           gbH=gbHu+gbHf                         !m/s
           rbH=1./gbH                            !b/l resistance to heat transfer
           rbw=0.93*rbH                          !b/l resistance to water vapour
!          Y factor for leaf: stom_n = 1.0 for hypostomatous leaf;  stom_n = 2.0 for amphistomatous leaf
           rbH_L=rbH*stom_n/2.                   !final b/l resistance for heat  
           rrdn=1./grdn
           Y=1./(1.+ (rbH_L+raero)/rrdn)
!          boundary layer conductance for CO2 - single side only (mol/m2/s)
           gbc=Cmolar*gbH/1.32            !mol/m2/s
           gsc0=gsw0/1.57                 !convert conductance for H2O to that for CO2
           varQc=0.0
           weighR=1.0
           call photosyn(Sps,CO2Ca,CO2Cs,Dleaf,Tlk,Qapar,Gbc,   &   !Qaparx<-Qapar,Gbcx<-Gsc0
     &         theta,a1,Ds0,fwsoil,varQc,weighR,                &
     &         gsc0,alpha,Vcmxx,eJmxx,weighJ,                   &
     &         conKc0,conKo0,Ekc,Eko,o2ci,Rconst,Trefk,         &
     &         Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,       &
     &         Aleafx,Gscx,gddonset)  !outputs
!          choose smaller of Ac, Aq
           Aleaf(ileaf) = Aleafx      !0.7 Weng 3/22/2006          !mol CO2/m2/s
!          calculate new values for gsc, cs (Lohammer model)
           co2cs = co2ca-Aleaf(ileaf)/gbc
           co2Ci(ileaf) = co2cs-Aleaf(ileaf)/gscx
!          scale variables
!           gsw=gscx*1.56      !gsw in mol/m2/s, oreginal:gsw=gsc0*1.56,Weng20060215
           gsw=gscx*1.56       !gsw in mol/m2/s, oreginal:gsw=gscx*1.56,Weng20090226
           gswv=gsw/Cmolar                           !gsw in m/s
           rswv=1./gswv
!          calculate evap'n using combination equation with current estimate of gsw
           Eleaf(ileaf)=1.0*                                    &
     &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    &   !2* Weng 0215
     &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))        

!          calculate sensible heat flux
           Hleaf(ileaf)=Y*(Rnstar(ileaf)-Eleaf(ileaf))
!          calculate new leaf temperature (K)
           Tlk1=273.2+Tair+Hleaf(ileaf)*(rbH/2.+raero)/rhocp
!          calculate Dleaf use LE=(rhocp/psyc)*gsw*Ds
           Dleaf=psyc*Eleaf(ileaf)/(rhocp*gswv)
           gbleaf(ileaf)=gbc*1.32*1.075
           gsleaf(ileaf)=gsw
!          compare current and previous leaf temperatures
           if(abs(Tlk1-Tlk).le.0.1) exit ! original is 0.05 C Weng 10/31/2008
!          update leaf temperature  ! leaf temperature calculation has many problems! Weng 10/31/2008
           Tlk=Tlk1
           Tleaf(ileaf)=Tlk1-273.2
           kr1=kr1+1
           if(kr1 > 500)then
               Tlk=TairK
               exit
           endif
           if(Tlk < 200.)then
                Tlk=TairK
                exit 
           endif                     ! Weng 10/31/2008
!        goto 100                          !solution not found yet
         enddo
! 10  continue
      enddo
      return
      end
!     ****************************************************************************
      subroutine agsean_ngt(Sps,Qabs,Rnstar,grdn,windUx,Tair,Dair,co2ca,    &
     &               wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,            &
     &               Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,      &
     &               gsw0,alpha,stom_n,                                     &
     &               Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,                &
     &               Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,             &
     &               Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci)
!    implicit real (a-z)
      integer kr1,ileaf
      real Aleaf(2),Eleaf(2),Hleaf(2),Tleaf(2),co2ci(2)
      real gbleaf(2), gsleaf(2)
      real Qabs(3,2),Rnstar(2)
!    thermodynamic parameters for air
      TairK=Tair+273.2
      rhocp=cpair*Patm*AirMa/(Rconst*TairK)
      H2OLv=H2oLv0-2.365e3*Tair
      slope=(esat(Tair+0.1)-esat(Tair))/0.1
      psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar=Patm/(Rconst*TairK)
      weighJ=1.0

!     boundary layer conductance for heat - single sided, forced convection
!     (Monteith 1973, P106 & notes dated 23/12/94)
      gbHu=0.003*sqrt(windUx/wleaf)    !m/s
!     raero=0.0                        !aerodynamic resistance s/m

      do ileaf=1,2                  ! loop over sunlit and shaded leaves
!        first estimate of leaf temperature - assume air temp
         Tleaf(ileaf)=Tair
         Tlk=Tleaf(ileaf)+273.2    !Tleaf to deg K
!        first estimate of deficit at leaf surface - assume Da
         Dleaf=Dair                !Pa
!        first estimate for co2cs
         co2cs=co2ca               !mol/mol
         Qapar = (4.6e-6)*Qabs(1,ileaf)
!        ********************************************************************
         kr1=0                     !iteration counter for LE
         do
!100        continue !    return point for evaporation iteration
!           single-sided boundary layer conductance - free convection (see notes 23/12/94)
            Gras=1.595e8*abs(Tleaf(ileaf)-Tair)*(wleaf**3)     !Grashof
            gbHf=0.5*Dheat*(Gras**0.25)/wleaf
            gbH=gbHu+gbHf                         !m/s
            rbH=1./gbH                            !b/l resistance to heat transfer
            rbw=0.93*rbH                          !b/l resistance to water vapour
!           Y factor for leaf: stom_n = 1.0 for hypostomatous leaf;  stom_n = 2.0 for amphistomatous leaf
            rbH_L=rbH*stom_n/2.                   !final b/l resistance for heat  
            rrdn=1./grdn
            Y=1./(1.+ (rbH_L+raero)/rrdn)
!           boundary layer conductance for CO2 - single side only (mol/m2/s)
            gbc=Cmolar*gbH/1.32            !mol/m2/s
            gsc0=gsw0/1.57                        !convert conductance for H2O to that for CO2
            varQc=0.0                  
            weighR=1.0
!           respiration      
            Aleafx=-0.0089*Vcmxx*exp(0.069*(Tlk-293.2))
            gsc=gsc0
!           choose smaller of Ac, Aq
            Aleaf(ileaf) = Aleafx                     !mol CO2/m2/s
!           calculate new values for gsc, cs (Lohammer model)
            co2cs = co2ca-Aleaf(ileaf)/gbc
            co2Ci(ileaf) = co2cs-Aleaf(ileaf)/gsc
!           scale variables
            gsw=gsc*1.56                              !gsw in mol/m2/s
            gswv=gsw/Cmolar                           !gsw in m/s
            rswv=1./gswv
!           calculate evap'n using combination equation with current estimate of gsw
            Eleaf(ileaf)=                                       &
     &      (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/   &
     &      (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
!           calculate sensible heat flux
            Hleaf(ileaf)=Y*(Rnstar(ileaf)-Eleaf(ileaf))
!           calculate new leaf temperature (K)
            Tlk1=273.2+Tair+Hleaf(ileaf)*(rbH/2.+raero)/rhocp
!           calculate Dleaf use LE=(rhocp/psyc)*gsw*Ds
            Dleaf=psyc*Eleaf(ileaf)/(rhocp*gswv)
            gbleaf(ileaf)=gbc*1.32*1.075
            gsleaf(ileaf)=gsw

!          compare current and previous leaf temperatures
            if(abs(Tlk1-Tlk).le.0.1)exit
            if(kr1.gt.500)exit
!           update leaf temperature
            Tlk=Tlk1 
            Tleaf(ileaf)=Tlk1-273.2
            kr1=kr1+1
         enddo                          !solution not found yet
10    continue
      enddo
      return
      end
!     ****************************************************************************
      subroutine ciandA(Gma,Bta,g0,X,Rd,co2Cs,gammas,ciquad,Aquad)
!     calculate coefficients for quadratic equation for ci
      b2 = g0+X*(Gma-Rd)
      b1 = (1.-co2cs*X)*(Gma-Rd)+g0*(Bta-co2cs)-X*(Gma*gammas+Bta*Rd)
      b0 = -(1.-co2cs*X)*(Gma*gammas+Bta*Rd)-g0*Bta*co2cs

      bx=b1*b1-4.*b2*b0
      if(bx.gt.0.0) then 
!       calculate larger root of quadratic
        ciquad = (-b1+sqrt(bx))/(2.*b2)
      endif

      IF(ciquad.lt.0.or.bx.lt.0.) THEN
        Aquad = 0.0
        ciquad = 0.7 * co2Cs
      ELSE
        Aquad = Gma*(ciquad-gammas)/(ciquad+Bta)
      ENDIF
      return
      end

!****************************************************************************
      subroutine goud1(FLAIT,coszen,radabv,fbeam,               &
     &                  Tair,eairP,emair,emsoil,emleaf,sigma,   &
     &                  tauL,rhoL,rhoS,xfang,extkb,extkd,       &
     &                  reffbm,reffdf,extkbm,extkdm,Qcan)
!    use the radiation scheme developed by
!    Goudriaan (1977, Goudriaan and van Larr 1995)
!=================================================================
!    Variable      unit      defintion
!    FLAIT         m2/m2     canopy leaf area index       
!    coszen                  cosine of the zenith angle of the sun
!    radabv(nW)    W/m2      incoming radiation above the canopy
!    fbeam                   beam fraction
!    fdiff                   diffuse fraction
!    funG(=0.5)              Ross's G function
!    extkb                   extinction coefficient for beam PAR
!    extkd                   extinction coefficient for diffuse PAR
!    albedo                  single scattering albedo
!    scatB                   upscattering parameter for beam
!    scatD                   upscattering parameter for diffuse
! ==================================================================
!    all intermediate variables in the calculation correspond
!    to the variables in the Appendix of of Seller (1985) with
!    a prefix of "x".
      integer nW
      real radabv(3)
      real rhocbm(3),rhocdf(3)
      real reffbm(3),reffdf(3),extkbm(3),extkdm(3)
      real tauL(3),rhoL(3),rhoS(3),scatL(3)
      real Qcan(3,2)
!
!     for PAR: using Goudriann approximation to account for scattering
      fdiff=1.0-fbeam
      xu=coszen
      xphi1 = 0.5 -0.633*xfang - 0.33*xfang*xfang
      xphi2 = 0.877 * (1.0 - 2.0*xphi1)
      funG = xphi1 + xphi2*xu
      extkb=funG/xu
                       
!     Effective extinction coefficient for diffuse radiation Goudriaan & van Laar Eq 6.6)
      pi180=3.1416/180.
      cozen15=cos(pi180*15)
      cozen45=cos(pi180*45)
      cozen75=cos(pi180*75)
      xK15=xphi1/cozen15+xphi2
      xK45=xphi1/cozen45+xphi2
      xK75=xphi1/cozen75+xphi2
      transd=0.308*exp(-xK15*FLAIT)+0.514*exp(-xK45*FLAIT)+     &
     &       0.178*exp(-xK75*FLAIT)
      extkd=(-1./FLAIT)*alog(transd)

!     canopy reflection coefficients (Array indices: 1=VIS,  2=NIR
      do nw=1,2                                                         !nw:1=VIS, 2=NIR
         scatL(nw)=tauL(nw)+rhoL(nw)                                    !scattering coeff
         if((1.-scatL(nw))<0.0) scatL(nw)=0.9999                        !Weng 10/31/2008
         extkbm(nw)=extkb*sqrt(1.-scatL(nw))                            !modified k beam scattered (6.20)
         extkdm(nw)=extkd*sqrt(1.-scatL(nw))                            !modified k diffuse (6.20)
         rhoch=(1.-sqrt(1.-scatL(nw)))/(1.+sqrt(1.-scatL(nw)))          !canopy reflection black horizontal leaves (6.19)
         rhoc15=2.*xK15*rhoch/(xK15+extkd)                              !canopy reflection (6.21) diffuse
         rhoc45=2.*xK45*rhoch/(xK45+extkd)
         rhoc75=2.*xK75*rhoch/(xK75+extkd)   
       
         rhocbm(nw)=2.*extkb/(extkb+extkd)*rhoch                        !canopy reflection (6.21) beam 
         rhocdf(nw)=0.308*rhoc15+0.514*rhoc45+0.178*rhoc75
         reffbm(nw)=rhocbm(nw)+(rhoS(nw)-rhocbm(nw))        &               !effective canopy-soil reflection coeff - beam (6.27)
     &             *exp(-2.*extkbm(nw)*FLAIT)                              
         reffdf(nw)=rhocdf(nw)+(rhoS(nw)-rhocdf(nw))        &            !effective canopy-soil reflection coeff - diffuse (6.27)
     &             *exp(-2.*extkdm(nw)*FLAIT)  

!        by the shaded leaves
         abshdn=fdiff*(1.0-reffdf(nw))*extkdm(nw)                       &           !absorbed NIR by shaded
     &      *(funE(extkdm(nw),FLAIT)-funE((extkb+extkdm(nw)),FLAIT))    &
     &      +fbeam*(1.0-reffbm(nw))*extkbm(nw)                          &
!    &      *(funE(extkbm(nw),FLAIT)-funE((extkb+extkdm(nw)),FLAIT))    ! error found by De Pury
     &      *(funE(extkbm(nw),FLAIT)-funE((extkb+extkbm(nw)),FLAIT))    &
     &      -fbeam*(1.0-scatL(nw))*extkb                                &
     &      *(funE(extkb,FLAIT)-funE(2.0*extkb,FLAIT))
!        by the sunlit leaves
         absltn=fdiff*(1.0-reffdf(nw))*extkdm(nw)                       &  !absorbed NIR by sunlit
     &      *funE((extkb+extkdm(nw)),FLAIT)                             &
     &      +fbeam*(1.0-reffbm(nw))*extkbm(nw)                          &
!    &      *funE((extkb+extkdm(nw)),FLAIT)                         ! error found by De Pury
     &      *funE((extkb+extkbm(nw)),FLAIT)                             &
     &      +fbeam*(1.0-scatL(nw))*extkb                                &
     &      *(funE(extkb,FLAIT)-funE(2.0*extkb,FLAIT))

!        scale to real flux 
!        sunlit    
          Qcan(nw,1)=absltn*radabv(nw)
!        shaded
          Qcan(nw,2)=abshdn*radabv(nw)
      enddo
!     
!    calculate the absorbed (iso)thermal radiation
      TairK=Tair+273.2
      
!     apparent atmospheric emissivity for clear skies (Brutsaert, 1975)
      emsky=0.642*(eairP/Tairk)**(1./7)      !note eair in Pa

!     apparent emissivity from clouds (Kimball et al 1982)
      ep8z=0.24+2.98e-12*eairP*eairP*exp(3000.0/TairK)
      tau8=amin1(1.0,1-ep8z*(1.4-0.4*ep8z))                !ensure tau8<1
      emcloud=0.36*tau8*(1.-fbeam)*(1-10./TairK)**4        !10 from Tcloud = Tair-10 

!     apparent emissivity from sky plus clouds      
!     emair=emsky+emcloud
! 20/06/96
      emair=emsky
      if(emair.gt.1.0) emair=1.0                             

      Bn0=sigma*(TairK**4)
      QLW1=-extkd*emleaf*(1.0-emair)*funE((extkd+extkb),FLAIT)      &
     &     -extkd*(1.0-emsoil)*(emleaf-emair)*exp(-2.0*extkd*FLAIT) &
     &     *funE((extkb-extkd),FLAIT)
      QLW2=-extkd*emleaf*(1.0-emair)*funE(extkd,FLAIT)              &
     &     -extkd*(1.0-emsoil)*(emleaf-emair)                       &
     &     *(exp(-extkd*FLAIT)-exp(-2.0*extkd*FLAIT))/extkd         &
     &     -QLW1
      Qcan(3,1)=QLW1*Bn0
      Qcan(3,2)=QLW2*Bn0
      return
      end

!****************************************************************************
      subroutine photosyn(Sps,CO2Ca,CO2Csx,Dleafx,Tlkx,Qaparx,Gbcx, &
     &         theta,a1,Ds0,fwsoil,varQc,weighR,                    &
     &         g0,alpha,                                            &
     &         Vcmx1,eJmx1,weighJ,conKc0,conKo0,Ekc,Eko,o2ci,       &
     &         Rconst,Trefk,Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,  &
     &         Aleafx,Gscx,gddonset)

!     calculate Vcmax, Jmax at leaf temp (Eq 9, Harley et al 1992)
!     turned on by Weng, 2012-03-13
!     VcmxT = Vjmax(Tlkx,Trefk,Vcmx1,Eavm,Edvm,Rconst,Entrpy)
!     eJmxT = Vjmax(Tlkx,Trefk,eJmx1,Eajm,Edjm,Rconst,Entrpy)
      CO2Csx=AMAX1(CO2Csx,0.6*CO2Ca)
!    check if it is dark - if so calculate respiration and g0 to assign conductance 
      if(Qaparx.le.0.) then                            !night, umol quanta/m2/s
        Aleafx=-0.0089*Vcmx1*exp(0.069*(Tlkx-293.2))   ! original: 0.0089 Weng 3/22/2006
        Gscx=g0
      endif
!     calculate  Vcmax, Jmax at leaf temp using Reed et al (1976) function J appl Ecol 13:925
      TminV=gddonset/10.  ! original -5.        !-Jiang Jiang 2015/10/13
      TmaxV=50.
      ToptV=35.
      
      TminJ=TminV
      TmaxJ=TmaxV
      ToptJ=ToptV 
      
      Tlf=Tlkx-273.2
      VcmxT=VJtemp(Tlf,TminV,TmaxV,ToptV,Vcmx1)
      eJmxT=VJtemp(Tlf,TminJ,TmaxJ,ToptJ,eJmx1)      
!     calculate J, the asymptote for RuBP regeneration rate at given Q
      eJ = weighJ*fJQres(eJmxT,alpha,Qaparx,theta)
!     calculate Kc, Ko, Rd gamma*  & gamma at leaf temp
      conKcT = EnzK(Tlkx,Trefk,conKc0,Rconst,Ekc)
      conKoT = EnzK(Tlkx,Trefk,conKo0,Rconst,Eko)
!     following de Pury 1994, eq 7, make light respiration a fixed proportion of
!     Vcmax
      Rd = 0.0089*VcmxT*weighR                              !de Pury 1994, Eq7
      Tdiff=Tlkx-Trefk
      gammas = gam0*(1.+gam1*Tdiff+gam2*Tdiff*Tdiff)       !gamma*
!     gamma = (gammas+conKcT*(1.+O2ci/conKoT)*Rd/VcmxT)/(1.-Rd/VcmxT)
      gamma = 0.0
!     ***********************************************************************
!     Analytical solution for ci. This is the ci which satisfies supply and demand
!     functions simultaneously
!     calculate X using Lohammer model, and scale for soil moisture
      a1= 1./(1.-0.7)
      X = a1*fwsoil/((co2csx - gamma)*(1.0 + Dleafx/Ds0))
!     calculate solution for ci when Rubisco activity limits A
      Gma = VcmxT  
      Bta = conKcT*(1.0+ o2ci/conKoT)
      call ciandA(Gma,Bta,g0,X,Rd,co2Csx,gammas,co2ci2,Acx)
!     calculate +ve root for ci when RuBP regeneration limits A
      Gma = eJ/4.
      Bta = 2.*gammas
!    calculate coefficients for quadratic equation for ci
      call ciandA(Gma,Bta,g0,X,Rd,co2Csx,gammas,co2ci4,Aqx)
!     choose smaller of Ac, Aq
      sps=AMAX1(0.001,sps)                  !Weng, 3/30/2006
      Aleafx = (amin1(Acx,Aqx) - Rd) !*sps     ! Weng 4/4/2006
!      if(Aleafx.lt.0.0) Aleafx=0.0    ! by Weng 3/21/2006
!    calculate new values for gsc, cs (Lohammer model)
      CO2csx = co2ca-Aleafx/Gbcx
      Gscx=g0 + X*Aleafx  ! revised by Weng
      return
      end
!***********************************************************************
      function funeJ(alpha,eJmxT,Qaparx)
      funeJ=alpha*Qaparx*eJmxT/(alpha*Qaparx+2.1*eJmxT)
      return
      end
!****************************************************************************
      real function esat(T)
!     returns saturation vapour pressure in Pa
      esat=610.78*exp(17.27*T/(T+237.3))
      return
      end

!****************************************************************************
      real function evapor(Td,Tw,Patm)
!* returns vapour pressure in Pa from wet & dry bulb temperatures
      gamma = (64.6 + 0.0625*Td)/1.e5
      evapor = esat(Tw)- gamma*(Td-Tw)*Patm
      return
      end

!****************************************************************************
      real function Vjmax(Tk,Trefk,Vjmax0,Eactiv,Edeact,Rconst,Entrop)
      anum = Vjmax0*EXP((Eactiv/(Rconst*Trefk))*(1.-Trefk/Tk))
      aden = 1. + EXP((Entrop*Tk-Edeact)/(Rconst*Tk))
      Vjmax = anum/aden
      return
      end
!****************************************************************************
      real function funE(extkbd,FLAIT)
      funE=(1.0-exp(-extkbd*FLAIT))/extkbd
      return
      end

!     ****************************************************************************
!     Reed et al (1976, J appl Ecol 13:925) equation for temperature response
!     used for Vcmax and Jmax
      real function VJtemp(Tlf,TminVJ,TmaxVJ,ToptVJ,VJmax0)
      if(Tlf.lt.TminVJ) Tlf=TminVJ   !constrain leaf temperatures between min and max
      if(Tlf.gt.TmaxVJ) Tlf=TmaxVJ
      pwr=(TmaxVJ-ToptVJ)/(ToptVj-TminVj)
      VJtemp=VJmax0*((Tlf-TminVJ)/(ToptVJ-TminVJ))*     &
     &       ((TmaxVJ-Tlf)/(TmaxVJ-ToptVJ))**pwr 
      return
      end

!     ****************************************************************************
      real function fJQres(eJmx,alpha,Q,theta)
      AX = theta                                 !a term in J fn
      BX = alpha*Q+eJmx                          !b term in J fn
      CX = alpha*Q*eJmx                          !c term in J fn
      if((BX*BX-4.*AX*CX)>=0.0)then
          fJQres = (BX-SQRT(BX*BX-4.*AX*CX))/(2*AX)
      else
          fJQres = (BX)/(2*AX)                   !Weng 10/31/2008
      endif

      return
      end

!     *************************************************************************
      real function EnzK(Tk,Trefk,EnzK0,Rconst,Eactiv)

      temp1=(Eactiv/(Rconst* Trefk))*(1.-Trefk/Tk)
!      if (temp1<50.)then
      EnzK = EnzK0*EXP((Eactiv/(Rconst* Trefk))*(1.-Trefk/Tk))
!      else
!      EnzK = EnzK0*EXP(50.)                                          ! Weng 10/31/2008
!      endif

      return
      end

!     *************************************************************************
      real function sinbet(doy,lat,pi,timeh)
      real lat
!     sin(bet), bet = elevation angle of sun
!     calculations according to Goudriaan & van Laar 1994 P30
      rad = pi/180.
!     sine and cosine of latitude
      sinlat = sin(rad*lat)
      coslat = cos(rad*lat)
!     sine of maximum declination
      sindec=-sin(23.45*rad)*cos(2.0*pi*(doy+10.0)/365.0)
      cosdec=sqrt(1.-sindec*sindec)
!     terms A & B in Eq 3.3
      A = sinlat*sindec
      B = coslat*cosdec
      sinbet = A+B*cos(pi*(timeh-12.)/12.)
      return
      end

!     *************************************************************************
      subroutine yrday(doy,hour,lat,radsol,fbeam)
      real lat
      pi=3.14159256
      pidiv=pi/180.0
      slatx=lat*pidiv
      sindec=-sin(23.4*pidiv)*cos(2.0*pi*(doy+10.0)/365.0)
      cosdec=sqrt(1.-sindec*sindec)
      a=sin(slatx)*sindec
      b=cos(slatx)*cosdec
      sinbet=a+b*cos(2*pi*(hour-12.)/24.)
      solext=1370.0*(1.0+0.033*cos(2.0*pi*(doy-10.)/365.0))*sinbet
      
      tmprat=radsol/solext

      tmpR=0.847-1.61*sinbet+1.04*sinbet*sinbet
      tmpK=(1.47-tmpR)/1.66
      if(tmprat.le.0.22) fdiff=1.0
      if(tmprat.gt.0.22.and.tmprat.le.0.35) then
        fdiff=1.0-6.4*(tmprat-0.22)*(tmprat-0.22)
      endif
      if(tmprat.gt.0.35.and.tmprat.le.tmpK) then
        fdiff=1.47-1.66*tmprat
      endif
      if(tmprat.ge.tmpK) then
        fdiff=tmpR
      endif
      fbeam=1.0-fdiff
      if(fbeam.lt.0.0) fbeam=0.0
      return
      end


!===============================================================================
!========================================================
! Subroutine 1. Read parameters from file
    subroutine Getparameters(lat,longi,wsmax,wsmin,     &              
    &   LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax,           &
    &   SapR,SapS,SLA,GLmax,GRmax,Gsmax,stom_n,         &
    &   a1,Ds0,Vcmx0,extkU,xfang,alpha,                 &
    &   Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C,         &
    &   Tau_Micro,Tau_slowSOM,Tau_Passive,              &
    &   gddonset,Q10,RL0,Rs0,Rr0,parafile,              &
    &   r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi)                  !added for MEMCMC
    
    implicit none
    real lat,longi,wsmax,wsmin                  
    real LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax
    real SapR,SapS,SLA,GLmax,GRmax,Gsmax,stom_n
    real a1,Ds0,Vcmx0,extkU,xfang,alpha
    real Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C
    real Tau_Micro,Tau_slowSOM,Tau_Passive
    real gddonset, Q10,Rl0,Rs0,Rr0
    real r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi
    character(len=50) parafile,commts
    
    parafile=TRIM(parafile)
!   open and read input file for getting climate data
    open(10,file=parafile,status='old')
    read(10,11)commts
    read(10,*)lat,longi,wsmax,wsmin
    read(10,11)commts
    read(10,*)LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax    
    read(10,11)commts
    read(10,*)SapR,SapS,SLA,GLmax,GRmax,Gsmax,stom_n
    read(10,11)commts
    read(10,*)a1,Ds0,Vcmx0,extkU,xfang,alpha
    read(10,11)commts
    read(10,*)Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C,Tau_Micro,Tau_slowSOM,Tau_Passive
    read(10,11)commts
    read(10,*)gddonset,Q10,Rl0,Rs0,Rr0
    read(10,11)commts
    read(10,*)r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi                      !this line is for MCMEME
    
11  format(a132)
    close(10)
    return
    end

! ====================================================================
! Subroutine 2. Read climatic forcing from file   
    subroutine Getclimate(year_seq,doy_seq,hour_seq,          &
    &   forcing_data,climatefile,lines,yr_length)
    
    implicit none
    integer, parameter :: ilines=90000
    integer, parameter :: iiterms=7
    integer,dimension(ilines):: year_seq,doy_seq,hour_seq
    real forcing_data(iiterms,ilines)
    character(len=50) climatefile,commts
    integer m,n,istat1,lines,yr_length

    open(11,file=climatefile,status='old',ACTION='read',     &
    &     IOSTAT=istat1)
!     skip 2 lines of input met data file
      read(11,'(a160)') commts
    m=0  ! to record the lines in a file
    yr_length=0 ! to record years of a dataset
    do    ! read forcing files
        m=m+1
        read(11,*,IOSTAT=istat1)year_seq(m),      &
        &       doy_seq(m),hour_seq(m),           &
        &       (forcing_data(n,m),n=1,iiterms)
        if(istat1<0)exit
    enddo ! end of reading the forcing file
    lines=m-1
    yr_length=(year_seq(lines)-year_seq(1))+1
    close(11)    ! close forcing file
    return
    end
    
    
  !Subroutine  Read water table from file   Yuanyuan   
    subroutine Getwatertable(year_seq,doy_seq,hour_seq,          &
    &   water_table,watertablefile,lines,yr_length)
    
    implicit none
    integer, parameter :: ilines=90000
    integer,dimension(ilines):: year_seq,doy_seq,hour_seq
    real water_table(ilines)
    character(len=50) watertablefile,commts
    integer m,n,istat1,lines,yr_length
    real water_table_read
    integer year,doy,hour,istat111

    open(111,file=watertablefile,status='old',ACTION='read',     &
    &     IOSTAT=istat1)
!     skip 2 lines of input met data file
    read(111,'(a160)') commts
    m=0  ! to record the lines in a file
! ERROR At line 4371 of file /source/frc_01.f90 (unit = 111, file = 'fort.111')   
    do
        m=m+1
        read (111,*,IOSTAT=istat111)year,doy,hour,water_table_read
        if(istat111<0)exit
        water_table(m)=water_table_read   
         
    enddo
  
    
    close(111)    ! close watertable file
    return
    end 
    
 !Subroutine  Read snow_depth from file   Yuanyuan   
    subroutine Getsnowdepth(year_seq,doy_seq,hour_seq,          &
    &   snow_in,snowdepthfile,lines,yr_length)
    
    implicit none
    integer, parameter :: ilines=90000
    integer,dimension(ilines):: year_seq,doy_seq,hour_seq
    real snow_in(ilines)
    character(len=50) snowdepthfile,commts
    integer m,n,istat1,lines,yr_length
    real snow_depth_read
    integer year,doy,hour,istat111

    open(1111,file=snowdepthfile,status='old',ACTION='read',     &
    &     IOSTAT=istat1)
!     skip 2 lines of input met data file
    read(1111,'(a160)') commts
    m=0  ! to record the lines in a file
   
    do
        !write(*,*)'snow_inaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'
        m=m+1
        !write(*,*)'mmm',m
        read (1111,*,IOSTAT=istat111)year,doy,hour,snow_depth_read
        if(istat111<0)exit
        snow_in(m)=snow_depth_read   
        
    enddo
  
    
    close(1111)    ! close snow depth file
    return
    end    

    subroutine GetObsCH4Data(obs_CH4_MEMCMC,std_CH4,len3,obsfile3)
    implicit none
    real obs_CH4_MEMCMC(2,1000),std_CH4(2,1000)
    real CH4_flux,CH4_flux_sd
    integer m,istat4,len3,days
    character(len=80) commts,obsfile3

    open(15,file=obsfile3,status='old')

!	read observed data   daily
    read (15,901) commts		!format need to be understood
    m=0
    do
        read (15,*,IOSTAT=istat4)days,CH4_flux,CH4_flux_sd		!recble variables
        if(istat4<0) exit
                m=m+1
                obs_CH4_MEMCMC(1,m)=real(days)          !QUESTION what's real() for?
                obs_CH4_MEMCMC(2,m)=CH4_flux
                std_CH4(1,m)=real(days)
                std_CH4(2,m)=CH4_flux_sd
                if(days.lt.365*4)then                   !         Why (days.lt.365*4)
                    len3=m                              !   How many days of obs data
                endif
    enddo

901 format(a80)    

    close(15)
    return
    end	
	  
				  
    
!===============================================================================
! Subroutine 3. read observation data from files
    subroutine GetObsData(obs_spruce,std,obs_soilwater,     &
            &   len1,len2,obsfile1,obsfile2,obs_soilt,obsfile3,obs_watertable_d,obs_soilt_d,&
            &   obsfile4,obsfile5,len3)
         
    implicit none
    real obs_spruce(12,1000),std(12,1000),obs_soilwater(5,40000),obs_soilt(11,40000)
    real obs_watertable_d(2,4000),obs_soilt_d(8,4000)
    real foliage,fnpp,wood,wnpp,root,rnpp,soilc,phenology
    real foliage_sd,fnpp_sd,wood_sd,wnpp_sd,root_sd,rnpp_sd
    real soilc_sd,pheno_sd,gpp_sd,nee_sd,er_sd
    real gpp,nee,er,sw
    integer days,m,istat2,istat3,len1,len2,istat4,len3
    integer year,doy,hour,month
    character(len=80) commts,obsfile1,obsfile2,obsfile3,obsfile4,obsfile5
    real soiltsurf,soilt5,soilt20,soilt40,soilt80,soilt160,soilt200
    real wt

    open(12,file=obsfile1,status='old')
    open(13,file=obsfile2,status='old')
    open(133,file=obsfile3,status='old')
    open(134,file=obsfile4,status='old')
    open(135,file=obsfile5,status='old')
!   Observed hourly data
    read(12,901) commts
    m=0
    do
        read (12,*,IOSTAT=istat2)days,gpp,gpp_sd,nee,nee_sd,er,er_sd,   &
        &   foliage,foliage_sd,fnpp,fnpp_sd,wood,wood_sd,wnpp,wnpp_sd,  &
        &   root,root_sd,rnpp,rnpp_sd,soilc,soilc_sd,phenology,pheno_sd
        if(istat2<0)exit
            m=m+1
            obs_spruce(1,m)=real(days)
            obs_spruce(2,m)=gpp
            obs_spruce(3,m)=nee
            obs_spruce(4,m)=er
            obs_spruce(5,m)=foliage
            obs_spruce(6,m)=fnpp
            obs_spruce(7,m)=wood
            obs_spruce(8,m)=wnpp
            obs_spruce(9,m)=root
            obs_spruce(10,m)=rnpp
            obs_spruce(11,m)=soilc
            obs_spruce(12,m)=real(phenology)
            std(1,m)=real(days)
            std(2,m)=gpp_sd
            std(3,m)=nee_sd
            std(4,m)=er_sd
            std(5,m)=foliage_sd
            std(6,m)=fnpp_sd
            std(7,m)=wood_sd
            std(8,m)=wnpp_sd
            std(9,m)=root_sd
            std(10,m)=rnpp_sd
            std(11,m)=soilc_sd
            std(12,m)=real(pheno_sd)
            
            if(days.lt.365*4)then
                len1=m
            endif
    enddo
 
!   Observed daily soil water data
    read(13,901) commts
    m=0
    do
        m=m+1
        read (13,*,IOSTAT=istat3)year,doy,hour,month,sw
        if(istat3<0)exit
        obs_soilwater(1,m)=real(year)
        obs_soilwater(2,m)=real(doy)
        obs_soilwater(3,m)=real(hour)
        obs_soilwater(4,m)=real(month)
        obs_soilwater(5,m)=sw        
    enddo
    len2=m-1
    
    read(133,901) commts
    m=0
    do
        m=m+1
        read (133,*,IOSTAT=istat4)year,doy,hour,month,soiltsurf,soilt5,soilt20,soilt40,soilt80,soilt160,soilt200
        if(istat4<0)exit
        obs_soilt(1,m)=real(year)
        obs_soilt(2,m)=real(doy)
        obs_soilt(3,m)=real(hour)
        obs_soilt(4,m)=real(month)
        obs_soilt(5,m)=soiltsurf 
        obs_soilt(6,m)=soilt5 
        obs_soilt(7,m)=soilt20 
        obs_soilt(8,m)=soilt40 
        obs_soilt(9,m)=soilt80 
        obs_soilt(10,m)=soilt160 
        obs_soilt(11,m)=soilt200
        
        !write(*,*)'soilt',obs_soilt(11,m)
    enddo
    len2=m-1
    
    !   Observed daily soil water data
    read(134,901) commts
    m=0
    do
        m=m+1
        read (134,*,IOSTAT=istat3)doy,soiltsurf,soilt5,soilt20,soilt40,soilt80,soilt160,soilt200
        if(istat3<0)exit
        obs_soilt_d(1,m)=real(doy)
        obs_soilt_d(2,m)=soiltsurf
        obs_soilt_d(3,m)=soilt5
        obs_soilt_d(4,m)=soilt20
        obs_soilt_d(5,m)=soilt40    
        obs_soilt_d(6,m)=soilt80 
        obs_soilt_d(7,m)=soilt160
        obs_soilt_d(8,m)=soilt200
    enddo
    len3=m-1
    
    read(135,901) commts
    m=0
    do
        m=m+1
        read (135,*,IOSTAT=istat3)doy,wt
        if(istat3<0)exit
        obs_soilt_d(1,m)=real(doy)
        obs_watertable_d(2,m)=wt
        
    enddo
    
      
901 format(a80)    

    close(12)
    close(13)
    close(133)
    close(134)
    close(135)
    return
    end
      
      
! **********************************************************    
    subroutine getCov(gamma,covfile,npara)
    implicit none
    integer npara,i,k
    real gamma(npara,npara)    
    character(len=80) covfile
    
    open(14,file=covfile,status='old')

    do i=1,npara
        read (14,*)(gamma(i,k),k=1,npara)
    enddo  
    return
    end        
     
!========================================================================
!    	cost function for observed data
    
    subroutine costFObsNee(Simu_dailyflux,Simu_soilwater,   &
            &   obs_spruce,std,obs_soilwater,len1,len2,   &
            &   J_last,upgraded,Simu_soiltemp,obs_soilt,Simu_watertable,water_table,&
            &   Simu_dailysoilt,Simu_dailywatertable,obs_soilt_d,obs_watertable_d,len3,&
            &   Simu_dailyice,Simu_dailyCH4,obs_CH4_MEMCMC,std_CH4,len6,do_soilt_da,do_snow_da,do_watertable_da, &
            &   do_methane_da,do_co2_da)             
    implicit none
    real Simu_dailyflux(12,4000),Simu_soilwater(10,40000),Simu_soiltemp(11,40000)
    real Simu_watertable(1,40000),water_table(40000)
    real obs_spruce(12,1000),std(12,1000)
    real obs_soilwater(5,40000),obs_soilt(11,40000)
    real obs_soilt_d(8,4000),obs_watertable_d(2,4000)
    integer day
    real J_new,J_last,delta_J,var2
    real J_gpp,J_nee,J_er,J_sw,J_foliage,J_fnpp
    real J_wood,J_wnpp,J_root,J_rnpp,J_soilc,J_pheno
    real J_soilt1,J_soilt2,J_soilt3,J_soilt4,J_soilt5,J_soilt6,J_soilt7
    real tmp1,dObsSim,random_harvest
    integer i,upgraded,len1,len2,len3,len6
    integer j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,js1,js2,js3,js4,js5,js6,js7,jwt1,j_dwt1
    real r_num
    real J_watertable,J_dwatertable
    real Simu_dailysoilt(11,4000),Simu_dailywatertable(1,4000),Simu_dailyice(10,4000)
    
    real J_dsoilt1,J_dsoilt2,J_dsoilt3,J_dsoilt4,J_dsoilt5,J_dsoilt6,J_dsoilt7
    real j_ds1,j_ds2,j_ds3,j_ds4,j_ds5,j_ds6,j_ds7,J_CH4flux
    
    real Simu_dailyCH4(6,4000)
    real obs_CH4_MEMCMC(2,1000),std_CH4(2,1000)
    integer jch4
    logical do_soilt_da, do_snow_da, do_watertable_da,do_methane_da,do_co2_da  
                               
!    compute J_obs
   if (do_soilt_da) then 
    J_gpp=0.0
    J_nee=0.0
    J_er=0.0
    J_foliage=0.0
    J_fnpp=0.0
    J_wood=0.0
    J_wnpp=0.0
    J_root=0.0
    J_rnpp=0.0
    J_soilc=0.0
    J_pheno=0.0
    j1=0
    j2=0
    j3=0
    j4=0
    j5=0
    j6=0
    j7=0
    j8=0
    j9=0
    j10=0
    j11=0
    do i=1,len1
        day=int(obs_spruce(1,i))
        if(obs_spruce(2,i).gt.-999)then
            j1=j1+1
            dObsSim=Simu_dailyflux(1,day)-obs_spruce(2,i)
            J_gpp=J_gpp+(dObsSim*dObsSim)/(2*std(2,i)*std(2,i))
        endif
        if(obs_spruce(3,i).gt.-999)then
            j2=j2+1
            dObsSim=Simu_dailyflux(2,day)-obs_spruce(3,i)
            J_nee=J_nee+(dObsSim*dObsSim)/(2*std(3,i)*std(3,i))
        endif
        if(obs_spruce(4,i).gt.-999)then
            j3=j3+1
            dObsSim=Simu_dailyflux(3,day)-obs_spruce(4,i)
            J_er=J_er+(dObsSim*dObsSim)/(2*std(4,i)*std(4,i))
        endif
        if(obs_spruce(5,i).gt.-999)then
            j4=j4+1
            dObsSim=Simu_dailyflux(4,day)-obs_spruce(5,i)
            J_foliage=J_foliage+(dObsSim*dObsSim)/(2*std(5,i)*std(5,i))
        endif
        if(obs_spruce(6,i).gt.-999)then
            j5=j5+1
            dObsSim=Simu_dailyflux(5,day)-obs_spruce(6,i)
            J_fnpp=J_fnpp+(dObsSim*dObsSim)/(2*std(6,i)*std(6,i))
        endif
        if(obs_spruce(7,i).gt.-999)then
            j6=j6+1
            dObsSim=Simu_dailyflux(6,day)-obs_spruce(7,i)
            J_wood=J_wood+(dObsSim*dObsSim)/(2*std(7,i)*std(7,i))
        endif
        if(obs_spruce(8,i).gt.-999)then
            j7=j7+1
            dObsSim=Simu_dailyflux(7,day)-obs_spruce(8,i)
            J_wnpp=J_wnpp+(dObsSim*dObsSim)/(2*std(8,i)*std(8,i))
        endif
        if(obs_spruce(9,i).gt.-999)then
            j8=j8+1
            dObsSim=Simu_dailyflux(8,day)-obs_spruce(9,i)
            J_root=J_root+(dObsSim*dObsSim)/(2*std(9,i)*std(9,i))
        endif
        if(obs_spruce(10,i).gt.-999)then
            j9=j9+1
            dObsSim=Simu_dailyflux(9,day)-obs_spruce(10,i)
            J_rnpp=J_rnpp+(dObsSim*dObsSim)/(2*std(10,i)*std(10,i))
        endif
        if(obs_spruce(11,i).gt.-999)then
            j10=j10+1
            dObsSim=Simu_dailyflux(10,day)-obs_spruce(11,i)
            J_soilc=J_soilc+(dObsSim*dObsSim)/(2*std(11,i)*std(11,i))
        endif
        if(obs_spruce(12,i).gt.-999)then
            j11=j11+1
            dObsSim=Simu_dailyflux(11,day)-obs_spruce(12,i)
            J_pheno=J_pheno+(dObsSim*dObsSim)/(2*std(12,i)*std(12,i))
        endif
    enddo
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    J_gpp=J_gpp/real(j1)
    J_nee=J_nee/real(j2)
    J_er=J_er/real(j3)
    J_foliage=J_foliage/real(j4)
    J_fnpp=J_fnpp/real(j5)
    J_wood=J_wood/real(j6)
    J_wnpp=J_wnpp/real(j7)
    J_root=J_root/real(j8)
    J_rnpp=J_rnpp/real(j9)
    J_soilc=J_soilc/real(j10)
    J_pheno=J_pheno/real(j11)
        
    !J_new=J_gpp+J_er+J_foliage+J_wood   
    !J_new = J_gpp+J_er+J_foliage+J_wood+J_fnpp+J_wnpp
    !J_new = J_gpp+J_er+J_foliage+J_wood+J_root+J_soilc
    J_new=J_gpp+J_er+J_foliage+J_wood+J_fnpp+J_wnpp+J_root+J_rnpp+J_soilc+J_pheno  
   endif  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! soilwater cost fun!!!!!!!    
    J_sw=0.0
    do i=1,len2
        if(obs_soilwater(5,i).gt.0.)then
        dObsSim=Simu_soilwater(1,i)-obs_soilwater(5,i)
        J_sw=J_sw+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif      
    enddo
!    
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  !! soil temp cost fun!!!!!!!!!!!!!!!!!!1 
  if (do_soilt_da) then 
    J_soilt1=0.0
    J_soilt2=0.0
    J_soilt3=0.0
    J_soilt4=0.0
    J_soilt5=0.0
    J_soilt6=0.0
    J_soilt7=0.0
    
    js1=0
    js2=0
    js3=0
    js4=0
    js5=0
    js6=0
    js7=0
        
    do i=1,len2
        if(obs_soilt(5,i).gt.-999.)then
           js1=js1+1 
           dObsSim=Simu_soiltemp(1,i)-obs_soilt(5,i)
           J_soilt1=J_soilt1+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif
        if(obs_soilt(6,i).gt.-999.)then
           js2=js2+1 
           dObsSim=Simu_soiltemp(2,i)-obs_soilt(6,i)
           J_soilt2=J_soilt2+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif
        if(obs_soilt(7,i).gt.-999.)then
           js3=js3+1
           dObsSim=Simu_soiltemp(3,i)-obs_soilt(7,i)
           J_soilt3=J_soilt3+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif
        if(obs_soilt(8,i).gt.-999.)then
           js4=js4+1 
           dObsSim=Simu_soiltemp(4,i)-obs_soilt(8,i)
           J_soilt4=J_soilt4+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif
        if(obs_soilt(9,i).gt.-999.)then
            js5=js5+1 
           dObsSim=Simu_soiltemp(5,i)-obs_soilt(9,i)
           J_soilt5=J_soilt5+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif
        if(obs_soilt(10,i).gt.-999.)then
            js6=js6+1 
           dObsSim=Simu_soiltemp(6,i)-obs_soilt(10,i)
           J_soilt6=J_soilt6+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif
        if(obs_soilt(11,i).gt.-999.)then
           js7=js7+1 
           dObsSim=Simu_soiltemp(7,i)-obs_soilt(11,i)
           J_soilt7=J_soilt7+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif
    enddo  
    J_soilt1=J_soilt1/real(js1)
    J_soilt2=J_soilt2/real(js2)
    J_soilt3=J_soilt3/real(js3)
    J_soilt4=J_soilt4/real(js4)
    J_soilt5=J_soilt5/real(js5)
    J_soilt6=J_soilt6/real(js6)
    J_soilt7=J_soilt7/real(js7)   
    !J_new=J_soilt1+J_soilt2+J_soilt3+J_soilt4+J_soilt5+J_soilt6+J_soilt7
  endif 
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  !! soil temp cost fun!!!!!!!!!!!!!!!!!!1  
  if (do_soilt_da) then 
    J_dsoilt1=0.0
    J_dsoilt2=0.0
    J_dsoilt3=0.0
    J_dsoilt4=0.0
    J_dsoilt5=0.0
    J_dsoilt6=0.0
    J_dsoilt7=0.0
    
    j_ds1=0
    j_ds2=0
    j_ds3=0
    j_ds4=0
    j_ds5=0
    j_ds6=0
    j_ds7=0
        
    do i=1,len3
        if(obs_soilt_d(2,i).gt.-999.)then
           j_ds1=j_ds1+1 
           dObsSim=Simu_dailysoilt(1,i)-obs_soilt_d(2,i)
           J_dsoilt1=J_dsoilt1+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif
        if(obs_soilt_d(3,i).gt.-999.)then
           j_ds2=j_ds2+1 
           dObsSim=Simu_dailysoilt(2,i)-obs_soilt_d(3,i)
           J_dsoilt2=J_dsoilt2+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif
        if(obs_soilt_d(4,i).gt.-999.)then
           j_ds3=j_ds3+1
           dObsSim=Simu_dailysoilt(3,i)-obs_soilt_d(4,i)
           J_dsoilt3=J_dsoilt3+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif
        if(obs_soilt_d(5,i).gt.-999.)then
           j_ds4=j_ds4+1 
           dObsSim=Simu_dailysoilt(4,i)-obs_soilt_d(5,i)
           J_dsoilt4=J_dsoilt4+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif
        if(obs_soilt_d(6,i).gt.-999.)then
            j_ds5=j_ds5+1 
           dObsSim=Simu_dailysoilt(5,i)-obs_soilt_d(6,i)
           J_dsoilt5=J_dsoilt5+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif
        if(obs_soilt_d(7,i).gt.-999.)then
            j_ds6=j_ds6+1 
           dObsSim=Simu_dailysoilt(6,i)-obs_soilt_d(7,i)
           J_dsoilt6=J_dsoilt6+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif
        if(obs_soilt_d(8,i).gt.-999.)then
           j_ds7=j_ds7+1 
           dObsSim=Simu_dailysoilt(7,i)-obs_soilt_d(8,i)
           J_dsoilt7=J_dsoilt7+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif
    enddo  
    J_dsoilt1=J_dsoilt1/real(j_ds1)
    J_dsoilt2=J_dsoilt2/real(j_ds2)
    J_dsoilt3=J_dsoilt3/real(j_ds3)
    J_dsoilt4=J_dsoilt4/real(j_ds4)
    J_dsoilt5=J_dsoilt5/real(j_ds5)
    J_dsoilt6=J_dsoilt6/real(j_ds6)
    J_dsoilt7=J_dsoilt7/real(j_ds7) 
    J_new=J_dsoilt2+J_dsoilt3+J_dsoilt4+J_dsoilt5+J_dsoilt6+J_dsoilt7!+J_dwatertable
  endif
!        
           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! water table cost function 
  if (do_watertable_da) then
    jwt1=0.
    do i=1,len2
        if(water_table(i).gt.-999.)then
           jwt1=jwt1+1 
           dObsSim=Simu_watertable(1,i)-water_table(i)
           J_watertable=J_watertable+(dObsSim*dObsSim)/(2*1.94*1.94)
        endif
    enddo
     J_watertable=J_watertable/real(jwt1)
    !J_new=J_watertable
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! daily water table  cost function 
  if (do_watertable_da) then
    j_dwt1=0.
    do i=1,len3
        if(obs_watertable_d(2,i).gt.-999.)then
           j_dwt1=j_dwt1+1 
           dObsSim=Simu_dailywatertable(1,i)-obs_watertable_d(2,i)
           J_dwatertable=J_dwatertable+(dObsSim*dObsSim)/(2*1.94*1.94*9.)
        endif
    enddo
   J_dwatertable=J_dwatertable/real(j_dwt1)
   J_new=J_dwatertable
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    
!!!!!!!!!!!!!! Methane cost function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (do_methane_da) then
    J_CH4flux=0.0
    jch4=0	
    do i=1,len6			!len1=m undefined here      m is the length of simu_dailyflux
        day=int(obs_CH4_MEMCMC(1,i))        !index for simudailyCH4 data  QUESTION int()  real()
        jch4=jch4+1
        dObsSim=Simu_dailyCH4(1,day)-obs_CH4_MEMCMC(2,i)			!need to be modified here Simu_dailyCH4
        J_CH4flux=J_CH4flux+(dObsSim*dObsSim)/(2*std_CH4(2,i)*std_CH4(2,i))
    enddo	
    J_CH4flux=J_CH4flux/real(jch4)
    J_new=J_CH4flux
   endif
     
    delta_J=J_new-J_last         !    delta_J=(J_new-J_last)/J_last   
    CALL random_number(r_num)

!     delta_J=-1     !accept all samples, No data
    if(ISNAN(J_new))then
        write(*,*)'NaN return, upgraded', upgraded 
        return
    endif
    !if(AMIN1(1.0,exp(-delta_J)).gt.r_num)then
    if(-delta_J.gt.0.0)then
        upgraded=upgraded+1
        J_last=J_new    
    endif
    return
    end


!******************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Square root of a matrix							  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine racine_mat(M, Mrac,npara)

    integer npara,i
    real M(npara,npara),Mrac(npara,npara)
    real valpr(npara),vectpr(npara,npara)
    Mrac=0.
    call jacobi(M,npara,npara,valpr,vectpr,nrot)
    do i=1,npara
	if(valpr(i).ge.0.) then
            Mrac(i,i)=sqrt(valpr(i))
	else
            print*, 'WARNING!!! Square root of the matrix is undefined.'
            print*, ' A negative eigenvalue has been set to zero - results may be wrong'
            Mrac=M
            return
	endif
    enddo
    Mrac=matmul(matmul(vectpr, Mrac),transpose(vectpr))

end subroutine racine_mat      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Extraction of the eigenvalues and the eigenvectors !!
!! of a matrix (Numerical Recipes)					  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE jacobi(a,n,np,d,v,nrot)
INTEGER :: n,np,nrot
REAL :: a(np,np),d(np),v(np,np)
INTEGER, PARAMETER :: NMAX=500
INTEGER :: i,ip,iq,j
REAL :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      
do ip=1,n
	do iq=1,n
		v(ip,iq)=0.
	end do
	v(ip,ip)=1.
end do

do ip=1,n
	b(ip)=a(ip,ip)
	d(ip)=b(ip)
	z(ip)=0.
end do

nrot=0
do i=1,50
	sm=0.
	do ip=1,n-1
		do iq=ip+1,n
			sm=sm+abs(a(ip,iq))
		end do
	end do
	if(sm.eq.0.)return
	if(i.lt.4)then
		tresh=0.2*sm/n**2
	else
		tresh=0.
	endif
	do ip=1,n-1
		do iq=ip+1,n
			g=100.*abs(a(ip,iq))
			if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
				a(ip,iq)=0.
			else if(abs(a(ip,iq)).gt.tresh)then
				h=d(iq)-d(ip)
				if(abs(h)+g.eq.abs(h))then
					t=a(ip,iq)/h
				else
					theta=0.5*h/a(ip,iq)
					t=1./(abs(theta)+sqrt(1.+theta**2))
					if(theta.lt.0.) then
						t=-t
					endif
				endif
				c=1./sqrt(1+t**2)
				s=t*c
				tau=s/(1.+c)
				h=t*a(ip,iq)
				z(ip)=z(ip)-h
				z(iq)=z(iq)+h
				d(ip)=d(ip)-h
				d(iq)=d(iq)+h
				a(ip,iq)=0.
				do j=1,ip-1
					g=a(j,ip)
					h=a(j,iq)
					a(j,ip)=g-s*(h+g*tau)
					a(j,iq)=h+s*(g-h*tau)
				end do
				do j=ip+1,iq-1
					g=a(ip,j)
					h=a(j,iq)
					a(ip,j)=g-s*(h+g*tau)
					a(j,iq)=h+s*(g-h*tau)
				end do
				do j=iq+1,n
					g=a(ip,j)
					h=a(iq,j)
					a(ip,j)=g-s*(h+g*tau)
					a(iq,j)=h+s*(g-h*tau)
				end do
				do j=1,n
					g=v(j,ip)
					h=v(j,iq)
					v(j,ip)=g-s*(h+g*tau)
					v(j,iq)=h+s*(g-h*tau)
				end do
				nrot=nrot+1
			endif
		end do
	end do
	do ip=1,n
		b(ip)=b(ip)+z(ip)
		d(ip)=b(ip)
		z(ip)=0.
	end do
end do
print*, 'too many iterations in jacobi' 
return
END subroutine jacobi


!===================================================
!       generate new coefficents
        subroutine coefgenerate(coefac,coefmax,coefmin,coef,search_length,npara)
        
        integer npara
        real coefac(npara),coefmax(npara),coefmin(npara),coef(npara)
        real r,coefmid,random_harvest
        integer i
        real search_length
        do i=1,npara
999         continue
            CALL random_number(random_harvest)
            r=random_harvest-0.5
            coef(i)=coefac(i)+r*(coefmax(i)-coefmin(i))*search_length
            if(coef(i).gt.coefmax(i).or.coef(i).lt.coefmin(i))goto 999
        enddo
        return
        end
!============
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Generation of a random vector from a multivariate  !!
!! normal distribution with mean zero and covariance  !!
!! matrix gamma.									  !!
!! Beware!!! In order to improve the speed of the	  !!
!! algorithms, the subroutine use the Square root	  !!
!! matrix of gamma									  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gengaussvect(gamma_racine,xold,xnew,npara)

integer npara
real gamma_racine(npara,npara)
real x(npara),xold(npara),xnew(npara)

do i=1,npara
    x(i)=rangauss(25)
enddo

x = matmul(gamma_racine, x)
xnew = xold + x
end subroutine gengaussvect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Generation of a random number from a standard	  !!
!! normal distribution. (Numerical Recipes)           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function rangauss(idum)


integer idum
real v1, v2, r, fac, gset
real r_num

data iset/0/
if(iset==0) then
1	CALL random_number(r_num)
        v1=2.*r_num-1
        CALL random_number(r_num)
	v2=2.*r_num-1
	r=(v1)**2+(v2)**2
	if(r>=1) go to 1
	fac=sqrt(-2.*log(r)/r)
	gset=v1*fac
	rangauss=v2*fac
	iset=1
else
	rangauss=gset
	iset=0
end if

return
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! variance matrix of a matrix of data				  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine varcov(tab,varcovar,npara,ncov)

integer npara,ncov
real tab(ncov,npara),tab2(ncov,npara)
real varcovar(npara,npara)

call centre(tab,tab2,npara,ncov)

varcovar = matmul(transpose(tab2), tab2)*(1./real(ncov))

end subroutine varcov

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute the centered matrix, ie. the matrix minus  !!
!! the column means									  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine centre(mat,mat_out,npara,ncov)

    integer npara,ncov
    real mat(ncov,npara),mat_out(ncov,npara)
    real mean

do i=1,npara
    mat_out(:,i) = mat(:,i) - mean(mat(:,i),ncov)
enddo

end subroutine centre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mean of a vector									  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Function mean(tab,ncov)
    integer ncov
real tab(ncov)
real mean,mean_tt
mean_tt=0.
do i=1,ncov	
mean_tt=mean_tt+tab(i)/real(ncov)
enddo
mean=mean_tt
End Function


!        
        







