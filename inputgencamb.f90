program inputgencamb
! This code is to generate the input files of camb which will calculate the inputs
! of the Fisher.

  implicit none

!defining the parameters used for Fisher
integer:: numfid=0, m, n
double precision:: r, ns, nt, nrun, tau, ombh2, omch2, H, omk, As, percentage
double precision:: nu_0, epsilon_0, epsilon_h, epsilon_l
double precision:: paramsave
double precision, dimension(30):: params
character(len=80), dimension(30):: filenamep, filenamem, outfilep, outfilem

INCLUDE 'trunk/filename'

!!!!! section A
!!!!! Add the fiducial parameters in this section
!!!!! NOTE: Also in section B, consistent statements need to be added
!!!!! when add a paramter, type according to following example
!!!!    write(*,*)'r = (e.g 0.1)'
!!!!    read(*,*) r
!!!!    numfid=numfid+1
!!!!	params(numfid)=r


  write(*,*)'Enter the parameters for the fiducial model'
  write(*,*)'r = (e.g 0.01)'
  read(*,*) r
    numfid=numfid+1
	params(numfid)=r
  write(*,*)'ns = (e.g 0.9645)'
  read(*,*) ns
    numfid=numfid+1
	params(numfid)=ns
  write(*,*)'tau = (e.g 0.079)'
  read(*,*) tau
    numfid=numfid+1
	params(numfid)=tau
  write(*,*)'Omegabh^2 = (e.g 0.02225)'
  read(*,*) ombh2
    numfid=numfid+1
	params(numfid)=ombh2
  write(*,*)'Omegach^2 = (e.g 0.1198)'
  read(*,*) omch2
    numfid=numfid+1
	params(numfid)=omch2
  write(*,*)'H (in km/s/Mpc) = (e.g 67.27)'
  read(*,*) H
    numfid=numfid+1
	params(numfid)=H
  write(*,*)'Scalar Amplitude = (e.g 2.2065E-9)'
  read(*,*) As
    numfid=numfid+1
	params(numfid)=As
	write(*,*)'epsilon_0 = (0.5)'
  read(*,*) epsilon_0
    numfid=numfid+1
	params(numfid)=epsilon_0
	
  !!!! The following is for inflation consistency relation
  !!!! If you want to add parameter, add it above
  nt=-r/8.0	
	

!!!!! End section A	
!!!!!!!!!!!!!!!!!!!!	end fiducial parameters


  write(*,*)'Enter fractional change (e.g. 0.1 or 0.05 is good, 0.01 too large)'
  read(*,*) percentage
  write(*,*)'Processing...'

!!!!! section B
!!!!! to add paramter type according to the following
!!!!!     write(9,*)'initial_ratio(1)=',r
!---------------generate input for the fiducial model ------------------
  open(unit=9,file='cambin/cambin_fid.ini',status='unknown')
  rewind 9

  ! for fixed test !!!NOTE: the first statement of the ini.file need to start
  ! with NO space
  write(9,'(a)')'DEFAULT(params.ini)'
  !end for fixed test
  write(9,*)'output_root = ../cambout/camb_fid'
  write(9,*)'initial_ratio(1)=',r
  write(9,*)'scalar_spectral_index(1)=',ns
  write(9,*)'re_optical_depth=',tau
  write(9,*)'ombh2=',ombh2
  write(9,*)'omch2=',omch2
  write(9,*)'hubble=',H
  write(9,*)'scalar_amp(1)=',As
  write(9,*)'epsilon_0=',epsilon_0
  
  !!! The following is for inflation consistency relation
  !!! If want to add parameter, add them above
  write(9,*)'tensor_spectral_index(1)=',nt  
  
  close(9)
!!!!! End section B  
  

!!!! section C
!!!! To add parameter, after changing section A, 
!!!! type according the following example
!!!!    write(10,*)'initial_ratio(1)=',params(n)
!!!!	n=n+1
!!!! and:
!!!!    write(11,*)'initial_ratio(1)=',params(n)
!!!!	n=n+1
!!!! for variant of parameters both plus and minus
  do m=1,numfid
    open(unit=10, file=filenamep(m))
	open(unit=11, file=filenamem(m))
	rewind 10
	rewind 11
	write(10,'(a)')'DEFAULT(params.ini)' !!!for default setup
	write(11,'(a)')'DEFAULT(params.ini)' !!!for default setup
	write(10,*)outfilep(m)
	write(11,*)outfilem(m)
	paramsave=params(m)
	params(m)=paramsave*(1+percentage)
	n=1
    write(10,*)'initial_ratio(1)=',params(n)
	n=n+1
    write(10,*)'scalar_spectral_index(1)=',params(n)
	n=n+1
    write(10,*)'re_optical_depth=',params(n)
	n=n+1
    write(10,*)'ombh2=',params(n)
	n=n+1
    write(10,*)'omch2=',params(n)
	n=n+1
    write(10,*)'hubble=',params(n)
	n=n+1
    write(10,*)'scalar_amp(1)=',params(n)
	n=n+1
	write(10,*)'epsilon_0=',params(n)
	n=n+1
	
    !!! The following is for inflation consistency relation
    !!! If want to add parameter, add them above
	nt=-params(1)/8
    write(10,*) 'tensor_spectral_index(1)=',nt	
	
	close(10)
    !!!!!now for minus
	params(m)=paramsave*(1-percentage)
	n=1
    write(11,*)'initial_ratio(1)=',params(n)
	n=n+1
    write(11,*)'scalar_spectral_index(1)=',params(n)
	n=n+1
    write(11,*)'re_optical_depth=',params(n)
	n=n+1
    write(11,*)'ombh2=',params(n)
	n=n+1
    write(11,*)'omch2=',params(n)
	n=n+1
    write(11,*)'hubble=',params(n)
	n=n+1
    write(11,*)'scalar_amp(1)=',params(n)
	n=n+1
	write(11,*)'epsilon_0',params(n)
	n=n+1
	
    !!! The following is for inflation consistency relation
    !!! If want to add parameter, add them above
	nt=-params(1)/8
    write(11,*) 'tensor_spectral_index(1)=',nt		
	
	params(m)=paramsave	
	close(11)
  end do
!!!! End section C
  

!!!!!!!!!!!saving the fiducial for later use
  open(unit=999,file='trunk/paramfid.dat',status='unknown')
  rewind 999
  write(999,*)numfid
  write(999,*)percentage
  do m=1,numfid
    write(999,*)params(m)
  end do
  
  write(*,*)'Done!'


  

  end program inputgencamb