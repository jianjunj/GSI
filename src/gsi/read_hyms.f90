subroutine read_hyms(mype,val_tovs,ithin,isfcalc,&
     rmesh,jsatid,gstime,infile,lunout,obstype,&
     nread,ndata,nodata,twind,sis, &
     mype_root,mype_sub,npe_sub,mpi_comm_sub,nobs, &
     nrec_start,nrec_start_ears,nrec_start_db,dval_use,radmod)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_hyms                 read hyms 1b data
!   prgmmr:         Jin          org: ncep/emc                date: 2026-11-05
!
! abstract:  This routine reads BUFR format HYMS (hyperspectral microwave
!            sounder) brightness temperature files (BUFR message type
!            NC021000).  It is derived from read_mws.f90.
!
!            Each BUFR subset is one field of view (FOVN 0-39, 40 per scan)
!            with a replicated block (HYMSCH) of 1759 channels
!            (about 49.36-57.91 GHz).  Only the channels listed in the
!            satinfo file for this sensor (sis) are kept; the BUFR
!            channel number (CHNM) is matched to the satinfo channel number.
!
!            The data are thinned to a specified resolution using simple
!            quality control checks.
!
!            When running the gsi in regional mode, the code only
!            retains those observations that fall within the regional
!            domain
!
! program history log:
!  2026-11-05  Modified from read_mws.f90
!
!   input argument list:
!     mype     - mpi task id
!     val_tovs - weighting factor applied to super obs
!     ithin    - flag to thin data
!     isfcalc  - method to calculate surface fields within FOV.
!                FOV-based calculation is not available for HYMS;
!                it is reset to zero (bilinear interpolation).
!     rmesh    - thinning mesh size (km)
!     jsatid   - satellite to read
!     gstime   - analysis time in minutes from reference date
!     infile   - unit from which to read BUFR data
!     lunout   - unit to which to write data for further processing
!     obstype  - observation type to process
!     twind    - input group time window(hours)
!     sis      - sensor/instrument/satellite indicator
!     mype_root - "root" task for sub-communicator
!     mype_sub - mpi task id within sub-communicator
!     npe_sub  - number of data read tasks
!     mpi_comm_sub - sub-communicator for data read
!     nrec_start - first subset with useful information
!     nrec_start_ears - first ears subset with useful information
!     nrec_start_db - first db subset with useful information
!
!   output argument list:
!     nread    - number of BUFR HYMS 1b observations read
!     ndata    - number of BUFR HYMS 1b profiles retained for further processing
!     nodata   - number of BUFR HYMS 1b observations retained for further processing
!     nobs     - array of observations on each subdomain for each processor
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,r_double,i_kind
  use satthin, only: super_val,itxmax,makegrids,destroygrids,checkob, &
      finalcheck,map2tgrid,score_crit
  use satthin, only: radthin_time_info,tdiff2crit
  use obsmod,  only: time_window_max
  use radinfo, only: iuse_rad,nuchan,nusis,jpch_rad, &
      use_edges,radedge1,radedge2,radstart,radstep,maxscan
  use gridmod, only: diagnostic_reg,regional,nlat,nlon,tll2xy,txy2ll,rlats,rlons
  use constants, only: deg2rad,zero,one,two,rad2deg,r60inv,r100,rearth_equator
  use crtm_module, only : max_sensor_zenith_angle
  use gsi_4dvar, only: l4dvar,l4densvar,iwinbgn,winlen
  use deter_sfc_mod, only: deter_sfc
  use gsi_nstcouplermod, only: nst_gsi,nstinfo
  use gsi_nstcouplermod, only: gsi_nstcoupler_skindepth,gsi_nstcoupler_deter
  use mpimod, only: npe
  use radiance_mod, only: rad_obs_type

  implicit none

! Declare passed variables
  character(len=*),intent(in   ) :: infile,obstype,jsatid
  character(len=20),intent(in  ) :: sis
  integer(i_kind) ,intent(in   ) :: mype,lunout,ithin
  integer(i_kind) ,intent(in   ) :: nrec_start,nrec_start_ears,nrec_start_db
  integer(i_kind) ,intent(inout) :: isfcalc
  integer(i_kind) ,intent(inout) :: nread
  integer(i_kind),dimension(npe) ,intent(inout) :: nobs
  integer(i_kind) ,intent(  out) :: ndata,nodata
  real(r_kind)    ,intent(in   ) :: rmesh,gstime,twind
  real(r_kind)    ,intent(inout) :: val_tovs
  integer(i_kind) ,intent(in   ) :: mype_root
  integer(i_kind) ,intent(in   ) :: mype_sub
  integer(i_kind) ,intent(in   ) :: npe_sub
  integer(i_kind) ,intent(in   ) :: mpi_comm_sub
  logical         ,intent(in   ) :: dval_use
  type(rad_obs_type),intent(in ) :: radmod

! Declare local parameters

  integer(i_kind),parameter:: n1bhdr=11
  integer(i_kind),parameter:: n2bhdr=4
  integer(i_kind),parameter:: maxobs = 1000000
! Upper bound on the number of replicated channels (1759 in the data)
  integer(i_kind),parameter:: max_bufr_chan = 2048
  integer(i_kind),parameter:: kidsat = 269       ! SAID: Spire Lemur 3U cubesat
  real(r_kind),parameter:: r360=360.0_r_kind
  real(r_kind),parameter:: tbmin=50.0_r_kind
  real(r_kind),parameter:: tbmax=550.0_r_kind
  real(r_kind),parameter:: bmiss_chk=1.0e9_r_kind
  ! The next two are one minute in hours
  real(r_kind),parameter:: one_minute=0.01666667_r_kind
  real(r_kind),parameter:: minus_one_minute=-0.01666667_r_kind

! Declare local variables
  logical outside,iuse,assim

  character(40) :: infile2
  character(8) subset
  character(80) hdr1b,hdr2b

  integer(i_kind) ireadsb,ireadmg
  integer(i_kind) i,j,k,ntest,iob,llll
  integer(i_kind) nscan
  integer(i_kind) iret,idate,nchanl,n,idomsfc(1)
  integer(i_kind) maxinfo
  integer(i_kind) nmind,itx,nreal,nele,itt,num_obs
  integer(i_kind) iskip
  integer(i_kind) lnbufr,ksatid,isflg
  integer(i_kind) ilat,ilon, ifovmod, nadir
  integer(i_kind),dimension(5):: idate5
  integer(i_kind):: ierr
  integer(i_kind):: radedge_min, radedge_max
  integer(i_kind), POINTER :: ifov
  integer(i_kind), TARGET :: ifov_save(maxobs)
  integer(i_kind) :: chan_map(max_bufr_chan)

  real(r_kind) sfcr
  real(r_kind),dimension(0:3):: sfcpct
  real(r_kind),dimension(0:3):: ts
  real(r_kind) :: tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10
  real(r_kind) :: zob,tref,dtw,dtc,tz_tr
  real(r_kind) :: satellite_height, rato

  real(r_kind) dlat,panglr,dlon,tdiff
  real(r_kind) dlon_earth_deg,dlat_earth_deg
  real(r_kind) step,start,dist1
  real(r_kind) lzaest
  real(r_kind),dimension(0:4):: rlndsea
  real(r_kind),allocatable,dimension(:,:):: data_all
  real(r_kind), POINTER :: bt_in(:), crit1,rsat, t4dv, solzen, solazi
  real(r_kind), POINTER :: dlon_earth,dlat_earth,satazi, lza

  integer(i_kind), ALLOCATABLE, TARGET :: it_mesh_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: rsat_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: t4dv_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: dlon_earth_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: dlat_earth_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: crit1_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: lza_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: satazi_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: solzen_save(:) 
  real(r_kind), ALLOCATABLE, TARGET :: solazi_save(:) 
  real(r_kind), ALLOCATABLE, TARGET :: bt_save(:,:)

  integer(i_kind),allocatable,dimension(:):: nrec
  real(r_double),allocatable,dimension(:,:):: data1b8
  real(r_double),dimension(n1bhdr):: bfr1bhdr
  real(r_double),dimension(n2bhdr):: bfr2bhdr

  real(r_kind) cdist,disterr,disterrmax,dlon00,dlat00

  real(r_kind)    :: ptime,timeinflat,crit0
  integer(i_kind) :: ithin_time,n_tbin
  integer(i_kind),pointer :: it_mesh => null()

  integer(i_kind) :: jchan,ibufr
!**************************************************************************
! Initialize variables

  maxinfo=31
  lnbufr = 15
  disterrmax=zero
  ntest=0
  ndata  = 0
  nodata  = 0
  nread  = 0

! Check obstype and satellite before allocating anything (e.g. thinning grids)
  if (obstype /= 'hyms') then
     write(6,*) 'READ_HYMS called for obstype '//obstype//': RETURNING'
     return
  end if

  if (jsatid /= 'lemur3u') then
     write(6,*) 'READ_HYMS: Unrecognized value for jsatid '//jsatid//': RETURNING'
     return
  end if

  ilon=3
  ilat=4
  nscan=0
  if(nst_gsi>0) then
     call gsi_nstcoupler_skindepth(obstype,zob)
  endif

  call radthin_time_info(obstype, jsatid, sis, ptime, ithin_time)
  if( ptime > 0.0_r_kind) then
     n_tbin=nint(2*time_window_max/ptime)
  else
     n_tbin=1
  endif
! Make thinning grids
  call makegrids(rmesh,ithin,n_tbin=n_tbin)

! 40 FOVs per scan (FOVN 0-39 in the BUFR file, used here as ifov=FOVN+1).
! The satellite zenith angle is smallest at FOVN=20, i.e. the scan passes
! nadir between ifov=20 and ifov=21; ifov<=nadir are on the negative side.
  nadir=20

! Select the satinfo channels for this sensor and map BUFR channel numbers
! (CHNM) onto their position in the satinfo list.
  chan_map = 0
  nchanl = 0
  do i=1,jpch_rad
     if (trim(nusis(i)) == trim(sis)) then
        if (nuchan(i) < 1 .or. nuchan(i) > max_bufr_chan) then
           write(6,*) 'READ_HYMS: satinfo channel number out of range: ',nuchan(i)
           cycle
        end if
        nchanl = nchanl + 1
        chan_map(nuchan(i)) = nchanl
     end if
  end do
  if (nchanl == 0) then
     write(6,*) 'READ_HYMS: no satinfo channels found for sis ',trim(sis),': RETURNING'
     call destroygrids
     return
  end if

  radedge_min = 0
  radedge_max = 1000
  do i=1,jpch_rad
     if (trim(nusis(i))==trim(sis)) then
        step  = radstep(i)
        start = radstart(i)
        if (radedge1(i)/=-1 .or. radedge2(i)/=-1) then
           radedge_min=radedge1(i)
           radedge_max=radedge2(i)
        end if
        exit 
     endif
  end do 

! Allocate arrays to hold all data for given satellite
  if(dval_use) maxinfo = maxinfo+2
  nreal = maxinfo + nstinfo
  nele  = nreal   + nchanl
  allocate(data_all(nele,itxmax),nrec(itxmax))
  nrec=999999

! FOV-based surface calculation is not available for HYMS
  if (isfcalc == 1) then
     write(6,*) 'READ_HYMS: isfcalc=1 not supported, using bilinear surface method'
     isfcalc = 0
  endif

!   Set rlndsea for types we would prefer selecting
  rlndsea(0) = zero
  rlndsea(1) = 15._r_kind
  rlndsea(2) = 10._r_kind
  rlndsea(3) = 15._r_kind
  rlndsea(4) = 100._r_kind
     
! If all channels of a given sensor are set to monitor or not
! assimilate mode (iuse_rad<1), reset relative weight to zero.
! We do not want such observations affecting the relative
! weighting between observations within a given thinning group.

  assim=.false.
  search: do i=1,jpch_rad
     if ((nusis(i)==sis) .and. (iuse_rad(i)>0)) then
        assim=.true.
        exit search
     endif
  end do search
  if (.not.assim) val_tovs=zero

! Allocate arrays for BUFR I/O
  ALLOCATE(data1b8(3,max_bufr_chan))
  ALLOCATE(rsat_save(maxobs))
  ALLOCATE(t4dv_save(maxobs))
  ALLOCATE(dlon_earth_save(maxobs))
  ALLOCATE(dlat_earth_save(maxobs))
  ALLOCATE(crit1_save(maxobs))
  ALLOCATE(it_mesh_save(maxobs))
  ALLOCATE(lza_save(maxobs))
  ALLOCATE(satazi_save(maxobs))
  ALLOCATE(solzen_save(maxobs)) 
  ALLOCATE(solazi_save(maxobs)) 
  ALLOCATE(bt_save(nchanl,maxobs))

  iob=1
! Big loop over standard data feed and possible rars/db data
! llll=1 normal feed, llll=2 RARS/EARS data, llll=3 DB/UW data
  ears_db_loop: do llll= 1, 3

     if(llll == 1)then
        infile2 = trim(infile)         ! Set bufr subset names based on type of data to read
     elseif(llll == 2) then
        infile2 = trim(infile)//'ears' ! Set bufr subset names based on type of data to read
     elseif(llll == 3) then
        infile2 = trim(infile)//'_db'  ! Set bufr subset names based on type of data to read
     end if

!    Reopen unit to satellite bufr file
     open(lnbufr,file=trim(infile2),form='unformatted',status = 'old', &
         iostat = ierr)
     if(ierr /= 0) cycle ears_db_loop

     call openbf(lnbufr,'IN',lnbufr)
!    Mnemonics follow BUFR table NC021000 (HYMS brightness temperature data):
!    1 SAID, 2 FOVN, 3-8 YEAR..SECO, 9 CLATH, 10 CLONH, 11 HMSL (altitude, m)
!    SIID is MISSING in the data, so it is not used to select the instrument.
     hdr1b ='SAID FOVN YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH HMSL'
     hdr2b ='SAZA SOZA BEARAZ SOLAZI'
   
!    Loop to read bufr file
     read_subset: do while(ireadmg(lnbufr,subset,idate)>=0 .AND. iob < maxobs)

        read_loop: do while (ireadsb(lnbufr)==0 .and. iob < maxobs)

           rsat       => rsat_save(iob)
           t4dv       => t4dv_save(iob)
           dlon_earth => dlon_earth_save(iob)
           dlat_earth => dlat_earth_save(iob)
           crit1      => crit1_save(iob)
           it_mesh    => it_mesh_save(iob)
           ifov       => ifov_save(iob)
           lza        => lza_save(iob)
           satazi     => satazi_save(iob)
           solzen     => solzen_save(iob)
           solazi     => solazi_save(iob)

!          inflate selection value for ears_db data
           crit0 = 0.01_r_kind
           if ( llll > 1 ) crit0 = crit0 + r100 * real(llll,r_kind)

           call ufbint(lnbufr,bfr1bhdr,n1bhdr,1,iret,hdr1b)

!          Extract satellite id.  If not the one we want, read next message
           if (bfr1bhdr(1) > bmiss_chk) cycle read_subset
           rsat=bfr1bhdr(1) 
           ksatid=nint(bfr1bhdr(1))
           if(ksatid /= kidsat) cycle read_subset

!          Skip subsets with missing FOV number, time or location
           if (bfr1bhdr(2) > bmiss_chk) cycle read_loop
           if (any(bfr1bhdr(3:8) > bmiss_chk)) cycle read_loop

!          Extract observation location and other required information
           if(abs(bfr1bhdr(9)) <= 90._r_kind .and. abs(bfr1bhdr(10)) <= r360)then
              dlat_earth = bfr1bhdr(9)
              dlon_earth = bfr1bhdr(10)
           else
              cycle read_loop
           end if
           if(dlon_earth<zero)  dlon_earth = dlon_earth+r360
           if(dlon_earth>=r360) dlon_earth = dlon_earth-r360


!          Extract date information.  If time outside window, skip this obs
           idate5(1) = bfr1bhdr(3) !year
           idate5(2) = bfr1bhdr(4) !month
           idate5(3) = bfr1bhdr(5) !day
           idate5(4) = bfr1bhdr(6) !hour
           idate5(5) = bfr1bhdr(7) !minute
           call w3fs21(idate5,nmind)
           t4dv= (real((nmind-iwinbgn),r_kind) + bfr1bhdr(8)*r60inv)*r60inv    ! add in seconds
           tdiff=t4dv+(iwinbgn-gstime)*r60inv
           if (l4dvar.or.l4densvar) then
              if (t4dv<minus_one_minute .OR. t4dv>winlen+one_minute) &
                  cycle read_loop
           else
              if(abs(tdiff) > twind+one_minute) cycle read_loop
           endif

           timeinflat=two
           call tdiff2crit(tdiff,ptime,ithin_time,timeinflat,crit0,crit1,it_mesh)
 
           call ufbint(lnbufr,bfr2bhdr,n2bhdr,1,iret,hdr2b)
           if (any(bfr2bhdr > bmiss_chk)) cycle read_loop

           satazi=bfr2bhdr(3)
           if (abs(satazi) > r360) then
              satazi=zero
           endif

!          FOVN is zero-based in the BUFR data
           ifov = nint(bfr1bhdr(2)) + 1
           lza = bfr2bhdr(1)*deg2rad      ! local zenith angle
           if(ifov <= nadir)    lza=-lza
           if(ifov == nadir) nscan = nscan+1 

           panglr=(start+real(ifov-1,r_kind)*step)*deg2rad
           satellite_height=bfr1bhdr(11)
!          Ensure orbit height is reasonable (data are near 530 km)
           if (satellite_height < 400000.0_r_kind .OR. &
              satellite_height > 700000.0_r_kind) satellite_height = 530000.0_r_kind
           rato = one + satellite_height/rearth_equator
           lzaest = asin(rato*sin(panglr))

           if(abs(lza)*rad2deg > max_sensor_zenith_angle) then
              write(6,*)'READ_HYMS WARNING lza error ',lza,panglr
              cycle read_loop
           end if

!          Check for errors in satellite zenith angles 
           if(abs(lzaest-lza)*rad2deg > one) then
              write(6,*)' READ_HYMS WARNING uncertainty in lza ', &
              lza*rad2deg,lzaest*rad2deg,sis,ifov,start,step
              cycle read_loop
           end if

           solzen_save(iob)=bfr2bhdr(2) 
           solazi_save(iob)=bfr2bhdr(4) 

!          Read the HYMSCH replication (CHNM, TMBR, HMSCHQ); iret is the
!          number of replications actually present (1759 in the data).
!          Pre-set to zero so entries not filled by ufbrep are ignored.
           data1b8 = 0.0_r_double
           call ufbrep(lnbufr,data1b8,3,max_bufr_chan,iret,'CHNM TMBR HMSCHQ')
           bt_save(1:nchanl,iob) = tbmin - one
           do ibufr = 1, iret
              if (data1b8(1,ibufr) < 1.0_r_double .or. &
                  data1b8(1,ibufr) > real(max_bufr_chan,r_double)) cycle
              jchan = chan_map(nint(data1b8(1,ibufr)))
              if (jchan < 1) cycle
              if (data1b8(2,ibufr) > bmiss_chk) cycle
              bt_save(jchan,iob) = data1b8(2,ibufr)
           end do

           iob=iob+1

        end do read_loop
     end do read_subset
     call closbf(lnbufr)
     close(lnbufr)
  end do ears_db_loop
  deallocate(data1b8)

  num_obs = iob-1

! Do not return early here: combine_radobs below is collective over the
! sub-communicator, so every task must reach it even with no data.
  if (num_obs <= 0) then
     write(6,*) 'READ_HYMS: No HYMS Data were read in'
     num_obs = 0
  end if
  if( mype_sub==mype_root) write(6,*) 'READ_HYMS: num_obs=', num_obs, ', nscan=', nscan, &
       ', nchanl=', nchanl

! Complete Read_HYMS thinning and QC steps

  ObsLoop: do iob = 1, num_obs  

     rsat       => rsat_save(iob)
     t4dv       => t4dv_save(iob)
     dlon_earth => dlon_earth_save(iob)
     dlat_earth => dlat_earth_save(iob)
     crit1      => crit1_save(iob)
     it_mesh    => it_mesh_save(iob)
     ifov       => ifov_save(iob)
     lza        => lza_save(iob)
     satazi     => satazi_save(iob)
     solzen     => solzen_save(iob)
     solazi     => solazi_save(iob)
     bt_in      => bt_save(1:nchanl,iob)
     
     dlat_earth_deg = dlat_earth
     dlon_earth_deg = dlon_earth
     dlat_earth = dlat_earth*deg2rad
     dlon_earth = dlon_earth*deg2rad   

!    Regional case
     if(regional)then
        call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
        if(diagnostic_reg) then
           call txy2ll(dlon,dlat,dlon00,dlat00)
           ntest=ntest+1
           cdist=sin(dlat_earth)*sin(dlat00)+cos(dlat_earth)*cos(dlat00)* &
                (sin(dlon_earth)*sin(dlon00)+cos(dlon_earth)*cos(dlon00))
           cdist=max(-one,min(cdist,one))
           disterr=acos(cdist)*rad2deg
           disterrmax=max(disterrmax,disterr)
        end if
           
!       Check to see if in domain
        if(outside) cycle ObsLoop
           
!    Global case
     else
        dlat=dlat_earth
        dlon=dlon_earth
        call grdcrd1(dlat,rlats,nlat,1)
        call grdcrd1(dlon,rlons,nlon,1)
     endif

! Check time window
     if (l4dvar.or.l4densvar) then
        if (t4dv<zero .OR. t4dv>winlen) cycle ObsLoop
     else
        tdiff=t4dv+(iwinbgn-gstime)*r60inv
        if(abs(tdiff) > twind) cycle ObsLoop
     endif
 
!    Map obs to thinning grid
     call map2tgrid(dlat_earth,dlon_earth,dist1,crit1,itx,ithin,itt,iuse,sis,it_mesh=it_mesh)
     if(.not. iuse)cycle ObsLoop

!
!    Check FOV and scan-edge usage
     if (.not. use_edges .and. (ifov < radedge_min .OR. ifov > radedge_max )) &
          cycle ObsLoop
     ifovmod=ifov
     if (ifovmod < 1 .or. ifovmod > maxscan) cycle ObsLoop
     nread=nread+nchanl
     
!    If any temperature exceeds limits, the channel is treated as bad.
     iskip=0
     do j=1,nchanl
        if (bt_in(j) < tbmin .or. bt_in(j) > tbmax) iskip = iskip + 1
     end do
     if (iskip >= nchanl) cycle ObsLoop

!    Determine surface properties based on 
!    sst and land/sea/ice mask   
!
!    isflg    - surface flag
!               0 sea
!               1 land
!               2 sea ice
!               3 snow
!               4 mixed                       

     call deter_sfc(dlat,dlon,dlat_earth,dlon_earth,t4dv,isflg, &
          idomsfc(1),sfcpct,ts,tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10,sfcr)

!    Compute "score" for observation.  All scores>=0.0.  Lowest score is "best".
!    No channel-based quality predictor is used for HYMS.
     crit1 = crit1 + rlndsea(isflg) + 10._r_kind*real(iskip,r_kind) + 0.01_r_kind * abs(zz)
     call checkob(dist1,crit1,itx,iuse)
     if(.not. iuse)cycle ObsLoop

     call finalcheck(dist1,crit1,itx,iuse)
     if(.not. iuse)cycle ObsLoop
     
!    interpolate NSST variables to Obs. location and get dtw, dtc, tz_tr
     if(nst_gsi>0) then
        tref  = ts(0)
        dtw   = zero
        dtc   = zero
        tz_tr = one
        if(sfcpct(0)>zero) then
           call gsi_nstcoupler_deter(dlat_earth,dlon_earth,t4dv,zob,tref,dtw,dtc,tz_tr)
        endif
     endif

! Re-calculate look angle
     panglr=(start+real(ifov-1,r_kind)*step)*deg2rad


!     Load selected observation into data array
              
     data_all(1 ,itx)= rsat                      ! satellite ID
     data_all(2 ,itx)= t4dv                      ! time
     data_all(3 ,itx)= dlon                      ! grid relative longitude
     data_all(4 ,itx)= dlat                      ! grid relative latitude
     data_all(5 ,itx)= lza                       ! local zenith angle
     data_all(6 ,itx)= satazi                    ! local azimuth angle
     data_all(7 ,itx)= panglr                    ! look angle
     data_all(8 ,itx)= ifovmod                   ! scan position
     data_all(9 ,itx)= solzen                    ! solar zenith angle
     data_all(10,itx)= solazi                    ! solar azimuth angle
     data_all(11,itx) = sfcpct(0)                ! sea percentage of
     data_all(12,itx) = sfcpct(1)                ! land percentage
     data_all(13,itx) = sfcpct(2)                ! sea ice percentage
     data_all(14,itx) = sfcpct(3)                ! snow percentage
     data_all(15,itx)= ts(0)                     ! ocean skin temperature
     data_all(16,itx)= ts(1)                     ! land skin temperature
     data_all(17,itx)= ts(2)                     ! ice skin temperature
     data_all(18,itx)= ts(3)                     ! snow skin temperature
     data_all(19,itx)= tsavg                     ! average skin temperature
     data_all(20,itx)= vty                       ! vegetation type
     data_all(21,itx)= vfr                       ! vegetation fraction
     data_all(22,itx)= sty                       ! soil type
     data_all(23,itx)= stp                       ! soil temperature
     data_all(24,itx)= sm                        ! soil moisture
     data_all(25,itx)= sn                        ! snow depth
     data_all(26,itx)= zz                        ! surface height
     data_all(27,itx)= idomsfc(1) + 0.001_r_kind ! dominate surface type
     data_all(28,itx)= sfcr                      ! surface roughness
     data_all(29,itx)= ff10                      ! ten meter wind factor
     data_all(30,itx) = dlon_earth_deg           ! earth relative longitude (deg)
     data_all(31,itx) = dlat_earth_deg           ! earth relative latitude (deg)
     
     if(dval_use) then
        data_all(32,itx)= val_tovs
        data_all(33,itx)= itt
     end if
     
     if(nst_gsi>0) then
        data_all(maxinfo+1,itx) = tref            ! foundation temperature
        data_all(maxinfo+2,itx) = dtw             ! dt_warm at zob
        data_all(maxinfo+3,itx) = dtc             ! dt_cool at zob
        data_all(maxinfo+4,itx) = tz_tr           ! d(Tz)/d(Tr)
     endif

     do i=1,nchanl
        data_all(i+nreal,itx)=bt_in(i)
     end do
     nrec(itx)=iob

  end do ObsLoop


! DEAllocate I/O arrays
  DEALLOCATE(rsat_save)
  DEALLOCATE(t4dv_save)
  DEALLOCATE(dlon_earth_save)
  DEALLOCATE(dlat_earth_save)
  DEALLOCATE(crit1_save)
  DEALLOCATE(it_mesh_save)
  DEALLOCATE(lza_save)
  DEALLOCATE(satazi_save)
  DEALLOCATE(solzen_save) 
  DEALLOCATE(solazi_save) 
  DEALLOCATE(bt_save)

  call combine_radobs(mype_sub,mype_root,npe_sub,mpi_comm_sub,&
       nele,itxmax,nread,ndata,data_all,score_crit,nrec)
  if( mype_sub==mype_root) write(6,*) 'READ_HYMS: after combine_obs, nread,ndata is ',nread,ndata
! 
  if(mype_sub==mype_root)then
     do n=1,ndata
        do i=1,nchanl
           if(data_all(i+nreal,n) > tbmin .and. &
                data_all(i+nreal,n) < tbmax)nodata=nodata+1
        end do
     end do

     if(dval_use .and. assim)then
        do n=1,ndata
           itt=nint(data_all(33,n))
           super_val(itt)=super_val(itt)+val_tovs
        end do
     end if
     
!    Write final set of "best" observations to output file
     call count_obs(ndata,nele,ilat,ilon,data_all,nobs)
     write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
     write(lunout) ((data_all(k,n),k=1,nele),n=1,ndata)
  end if
     
! Deallocate local arrays
  deallocate(data_all,nrec)

! Deallocate satthin arrays
  call destroygrids

  if(diagnostic_reg.and.ntest>0) write(6,*)'READ_HYMS:  ',&
     'mype,ntest,disterrmax=',mype,ntest,disterrmax

! End of routine
  return

end subroutine read_hyms
