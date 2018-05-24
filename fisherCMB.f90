! This program is to build CMB fisher matrix using the outputs of CAMB

Program fisherCMB
implicit none

integer:: numfid !number of fiducial parameters
real:: percentage !variant of the fiducial parameters 
double precision, dimension(2:3000,3,3):: R_fid=0, C_fid=0, R_fid_inv=0 !for l=2 to 3000, not neccessarily to 3000 (1,1) for TT (1,2) for TE, (2,2) for EE, (3,3) for BB
double precision:: detR
double precision, dimension(2:3000,3,3,2,30):: C_p, C_m !variant powerspectrum for plus and minus delta
character(len=80):: input_fid='cambout/camb_fid_lensedtotCls.dat'
character(len=80), dimension(2,30):: inputfiles !2 for p and m, 30 just to prepare up 30 params
character(len=80), dimension(80):: store !for storing the variant matrix
integer:: file_unit=9, storefile_unit=900 !unit of the file
integer:: l, m, n, l_BB_max=512 !multiple
double precision:: cltt, clee, clbb, clte !TT,EE,BB,TE
double precision:: clttp, cleep, clbbp, cltep !delta plus TT,EE,BB,TE
double precision:: clttm, cleem, clbbm, cltem !delta minus TT,EE,BB,TE
double precision:: dcltt, dclee, dclbb, dclte
double precision, dimension(30):: params=0
double precision, dimension(30,30):: fisher=0
double precision:: dfisher
double precision:: fsky
double precision:: r1,r2,r3,r12,cm1,cm2,cm3,cm12,cn1,cn2,cn3,cn12
double precision:: total_n_BB, N_TT, N_EE


!!!!!!!!!!!
!prepare for potential input files
INCLUDE 'trunk/storederivate'
!!!!!!!!!!!

write(*,*)'Enter the sky fraction (e.g. 0.65)'
read(*,*)fsky

!!!!!!!!!!read in fiducial parameters
open(unit=800,file='trunk/paramfid.dat',status='unknown')
read(800,*)numfid
read(800,*)percentage
do m=1,numfid
    read(800,*)params(m)
	end do
close(800)

!read in the fiducial power spectrum
open(unit=17,file=input_fid,status='unknown')
open(unit=207,file='total_uncertainty/total_BB_COrE.txt',status='unknown')
open(unit=208,file='total_uncertainty/N_TT_pre.txt',status='unknown')
open(unit=209,file='total_uncertainty/N_EE_pre.txt',status='unknown')

100 read(17,*,end=101)l, cltt, clee, clbb, clte
!!!!!!!!!!!!! Important! C_X in Camb is C_L*l*(l+1)/(2pi)
    C_fid(l,1,1)=cltt/l/(l+1)*2*3.14159265
	C_fid(l,1,2)=clte/l/(l+1)*2*3.14159265
	C_fid(l,2,1)=clte/l/(l+1)*2*3.14159265
	C_fid(l,2,2)=clee/l/(l+1)*2*3.14159265
	C_fid(l,3,3)=clbb/l/(l+1)*2*3.14159265
	read(208,*)N_TT
	read(209,*)N_EE
	R_fid(l,1,1)=C_fid(l,1,1)+N_TT
	R_fid(l,1,2)=C_fid(l,1,2)+N_EE
	R_fid(l,2,1)=C_fid(l,2,1)
	R_fid(l,2,2)=C_fid(l,2,2)
	! Using BB up to l_BB_maxx
	if (l<=l_BB_max) then 
	   read(207,*)total_n_BB
	   R_fid(l,3,3)=C_fid(l,3,3)+total_n_BB
	else
	   R_fid(l,3,3)=1.0
	end if
	!take the inverse of R_fid
	detR=(R_fid(l,1,1)*R_fid(l,2,2)-R_fid(l,1,2)*R_fid(l,2,1))*R_fid(l,3,3)
	R_fid_inv(l,1,1)=R_fid(l,2,2)*R_fid(l,3,3)/detR
	R_fid_inv(l,2,2)=R_fid(l,1,1)*R_fid(l,3,3)/detR
	R_fid_inv(l,1,2)=-R_fid(l,2,1)*R_fid(l,3,3)/detR
	R_fid_inv(l,2,1)=R_fid_inv(l,1,2)
    R_fid_inv(l,3,3)=1/R_fid(l,3,3)  
	!End taking the inverse of R_fid
	open(unit=storefile_unit,file='trunk/R_fid_inv.dat',status='unknown')
	write(storefile_unit,*)l, R_fid_inv(l,1,1),R_fid_inv(l,2,2),R_fid_inv(l,3,3),R_fid_inv(l,1,2)
    goto 100
101 continue
    close(207)
	close(208)
	close(209)
    close(17)
	close(storefile_unit)
    file_unit=file_unit+1
	storefile_unit=storefile_unit+1
	
    


!!!!!!!!!!!!!!!!read in variant inputs and calculate the fist derivative
!!!!!Read in
    do m=1,numfid
	    open(unit=file_unit,file=inputfiles(1,m),status='unknown')
	 	open(unit=file_unit+1,file=inputfiles(2,m),status='unknown')
		open(unit=storefile_unit,file=store(m),status='unknown')
	200 read(file_unit,*,end=201)l, clttp, cleep, clbbp, cltep
	    read(file_unit+1,*)l, clttm, cleem, clbbm, cltem
		!!!!!!!!!!!!! Important! C_X in Camb is C_L*l*(l+1)/(2pi)
		dcltt=(clttp-clttm)/(2*percentage*params(m))/l/(l+1)*2*3.14159265
		dclee=(cleep-cleem)/(2*percentage*params(m))/l/(l+1)*2*3.14159265
		dclbb=(clbbp-clbbm)/(2*percentage*params(m))/l/(l+1)*2*3.14159265
		dclte=(cltep-cltem)/(2*percentage*params(m))/l/(l+1)*2*3.14159265
		write(storefile_unit,*)dcltt, dclee, dclbb, dclte
	    goto 200
	201 continue
	    close(file_unit)
		close(file_unit+1)
		close(storefile_unit)
	    file_unit=file_unit+2
		storefile_unit=storefile_unit+1
	end do
!!!!!!!!!!!!!!!!End read in variant inputs and calculate the fist derivative

!!!!!!!!!!Calculate the fisher matrix
open(unit=600,file='trunk/R_fid_inv.dat')
!!! for diagonal elements
do m=1,numfid
open(unit=601,file=store(m),status='unknown')
rewind 600
301 read(600,*,end=300)l,r1,r2,r3,r12
    read(601,*)cm1,cm2,cm3,cm12
    dfisher=(r1*cm1+r12*cm12)**2
	dfisher=dfisher+2.0*(r1*cm12+r12*cm2)*(r12*cm1+r2*cm12)
	dfisher=dfisher+(r12*cm12+r2*cm2)**2
	if (l<=l_BB_max) then
	  dfisher=dfisher+r3*r3*cm3*cm3
	end if
	dfisher=(l+0.5)*dfisher
	fisher(m,m)=fisher(m,m)+dfisher
    go to 301
300 continue
    close(601)
	fisher(m,m)=fisher(m,m)*fsky
end do

!!!! for off diagonal elements
do m=1,numfid
    do n=m+1,numfid
	rewind 600
    open(unit=601,file=store(m),status='unknown')
    open(unit=602,file=store(n),status='unknown')
	401 read(600,*,end=400)l,r1,r2,r3,r12
	    read(601,*)cm1,cm2,cm3,cm12
	    read(602,*)cn1,cn2,cn3,cn12
		dfisher=(r1*cm1+r12*cm12)*(r1*cn1+r12*cn12)
		dfisher=dfisher+(r1*cm12+r12*cm2)*(r12*cn1+r2*cn12)
		dfisher=dfisher+(r1*cn12+r12*cn2)*(r12*cm1+r2*cm12)
		dfisher=dfisher+(r12*cm12+r2*cm2)*(r12*cn12+r2*cn2)
		if (l<=l_BB_max) then
		dfisher=dfisher+r3*r3*cm3*cn3
		end if
		dfisher=(l+0.5)*dfisher
		fisher(m,n)=fisher(m,n)+dfisher
	    goto 401
	400 continue
	close(601)
	close(602)
	fisher(m,n)=fsky*fisher(m,n)
	fisher(n,m)=fisher(m,n)
	end do
end do
!!!!End calculating fisher matrix


!!!!!!!output the fisher matrix in fisher.txt
open(unit=1000,file='fisher.txt')
do m=1,numfid
write(1000,66)fisher(m,1:numfid)
66 format(es13.6, es18.6, es18.6, es18.6, es18.6, es18.6, es18.6, es18.6)
end do
close(1000)

END program fisherCMB