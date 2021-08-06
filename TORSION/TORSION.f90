!Unless otherwise specified, all intermediate quantities are in J/mol or J/mol/K
program TORSION
use func
use TORANA
implicit real*8 (a-h,o-z)
real*8,parameter :: R=8.3144648D0, Kb=1.3806503D-23 !Ideal gas constant (J/mol/K), Boltzmann constant (J/K)
real*8,parameter :: NA=6.02214179D23 !Avogadro constant
real*8,parameter :: au2kcal_mol=627.51D0,au2KJ_mol=2625.5D0,au2J=4.359744575D-18,au2cm_1=219474.6363D0 !Hartree to various units
real*8,parameter :: cal2J=4.184D0 !This is consistent with G09, though 4.186 is more generally accepted
real*8,parameter :: wave2freq=2.99792458D10 !cm^-1 to s^-1 (Hz)
real*8,parameter :: h=6.62606896D-34 !Planck constant, in J*s
real*8,parameter :: amu2kg=1.66053878D-27
real*8,parameter :: pi=3.141592653589793D0
real*8,parameter :: b2a=0.52917720859D0 !Bohr to Angstrom
real*8,allocatable :: xyz(:,:),freq(:),wavenum(:),wavall(:),wav_proall(:),atmass(:),T(:),q_trans(:),U_trans(:),&
H_trans(:),S_trans(:),q_rot(:),U_rot(:),S_rot(:),U_vib(:),CV_vib(:),S_vib(:),qvib_v0(:),qvib_bot(:)
real*8:: inert1,inert2,inert3, inert11,inert21,inert31
character,allocatable :: elem(:)*2 
character :: c80*80,  Str*80 , filename*200
character*2 :: name2ind(1:118)=(/ "H ","He", &
"Li","Be","B ","C ","N ","O ","F ","Ne", & !3~10
"Na","Mg","Al","Si","P ","S ","Cl","Ar", & !11~18
"K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", & !19~36
"Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe", & !37~54
"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", & !55~71
"Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", & !72~86
"Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr", & !87~103
"Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Ut","Fl","Up","Lv","Us","Uo" /) !104~all  Such as Uuo/Uup is replaced by Uo/Up
! some aguments that read temperature needed    
real :: a , b 
integer :: c , j=0 ,k(10), i, stp
integer*8 :: natm, ntor
logical :: alive


write(*,*) "TORANA: A utility to calculate thermodynamic properties:"
write(*,*) "1.Thermodynamic properties with all of the motion modes"
write(*,*) "2.Thermodynamic properties without internal torsion"
write(*,*) "Programmed by YanLiu(IanPioneee@163.com)"
write(*,*) "First release: 2020-Feb-26  Last update: 2020-May-19"
write(*,*)

isilent=0
call getarg(1,filename)
if (trim(filename)=="") then
	write(*,*) "Input the path of a formated file"
	do while(.true.)
		read(*,"(a)") filename
		inquire(file=filename,exist=alive)
		if (alive) exit
		write(*,*) "Cannot find the file, input again"
	end do
else
	isilent=1
end if

!================load detail informations about molecule===============
open(10,file=filename,status="old")

call loclabel(10,"Temperature(K)")
read(10,*)
Read(10,*) Temperature1, Temperature2, stp
allocate(T(2+stp))
T(1)=298.15
do i=1,(1+stp)
    T(1+i)=Temperature1+(i-1)*Temperature2
end do

call loclabel(10,"Pressure(Atm)")
read(10,*)
read(10,*) P  !Pressure
write(*,"(' Pressure(Atm):',f12.3)") P
P=P*101325 !From Atm to Pa
call loclabel(10,"Rotational symmetry number ")
read(10,*)
read(10,*) nRotsym  !Rotational symmetry number
write(*,"(' Rotational symmetry number:',i2)") abs(nRotsym)

call loclabel(10,"Spin multiplicity")
read(10,*)
read(10,*) nSpinmulti  !Spin multiplicity
write(*,"(' Spin multiplicity:',i2)") nSpinmulti
call loclabel(10,"The number of atoms")
read(10,*)
read(10,*) natm  !The number of atoms
write(*,"(' The number of atoms:',i5)") natm
call loclabel(10,"The number of torsion bonds")
read(10,*)
read(10,*) ntor  !The number of torsion bond
write(*,"(' The number of torsion bonds:',i5)") ntor
!Load elements and XYZ coordinates
call loclabel(10,"Atomic coordinates")
allocate(xyz(natm,3),elem(natm))
read(10,*)
do iatm=1,natm
	read(10,*) elem(iatm),xyz(iatm,:)
end do
!load scale factor
call loclabel(10,"Scale factor")
read(10,*)
read(10,*) factor
write(*,"(' Scale factor:',f12.6)") factor

!Determine masses
call loclabel(10,"Mass for all atoms")
read(10,*)
allocate(atmass(natm))
mass=0
do while(.true.)
	read(10,*,iostat=ierror) c80,tmpmass
	if (ierror/=0) exit
	if (iachar(c80(1:1))>=48.and.iachar(c80(1:1))<=57) then
		read(c80,*) iatm
		atmass(iatm)=tmpmass
	else
		do iatm=1,natm
			if (c80(1:2)==elem(iatm)) atmass(iatm)=tmpmass 
		end do
	end if
end do
if (any(atmass==0)) write(*,*) "Warning: Mass of some atoms are not specified!"
do iatm=1,natm
	write(*,"(' Atom:',i5,'   Mass:',f8.3,' amu')") iatm,atmass(iatm)
end do
totmass=sum(atmass(:))
write(*,"(' Total mass:',f12.6,' amu')") totmass

!================load molecule informations end===============

!==============calculate frequency==========

allocate(wavall(3*natm),wav_proall(3*natm))
call ANALYSIS(natm,atmass,xyz,ntor,wavall,wav_proall,inert1,inert2,inert3)
write(*,*) "**********************************************************************"
write(*,*) "*            Now, thermodynamic properties will be printed           *"
write(*,*) "*Note: firstly, some properties of calculated molecule will be printed*"
write(*,*) "**********************************************************************"
write(*,*)
write(*,*)
	write(*,"(' Moments of inertia: ',3f12.6,' amu*Bohr^2')") inert1,inert2,inert3
	rotcst1=h/(8D0*pi**2*inert1*amu2kg*(b2a*1D-10)**2) !in 1Hz=1/s
	rotcst2=h/(8D0*pi**2*inert2*amu2kg*(b2a*1D-10)**2)
	rotcst3=h/(8D0*pi**2*inert3*amu2kg*(b2a*1D-10)**2)
	write(*,"(' Rotational constant:',3f12.6,' GHz')") rotcst1/1D9,rotcst2/1D9,rotcst3/1D9
	rotT1=rotcst1*h/kb
	rotT2=rotcst2*h/kb
	rotT3=rotcst3*h/kb
	write(*,"(' Rotational temperature:',3f12.6,' K')") rotT1,rotT2,rotT3

close(10)
!==============calculate freq end============

cir=0
outer: do cricle=1,2
    cir=cir+cricle
    if (cir==1) then
        nfreq=3*natm-6
        allocate(freq(nfreq),wavenum(nfreq))
        do ifreq=1, nfreq
            wavenum(ifreq)=wavall(ifreq+6)
            freq(ifreq)=wavenum(ifreq)*wave2freq
        end do
    write(*,*)
    write(*,*) "**************************************************"
    write(*,*) "*    Thermal properties with all motion modes    *"
    write(*,*) "**************************************************"
    write(*,*)
    else
        nfreq=3*natm-6-ntor
        allocate(freq(nfreq),wavenum(nfreq))
        do ifreq=1, nfreq
            wavenum(ifreq)=wav_proall(ifreq+6+ntor)
            freq(ifreq)=wavenum(ifreq)*wave2freq
        end do
    write(*,*)
    write(*,*) "**************************************************"
    write(*,*) "*    Thermal properties without torsion modes    *"
    write(*,*) "**************************************************"
    write(*,*)
    end if

idx=0
write(*,*)"*****  Frequencies and wavenumbers  *****"
do i=1,nfreq
	idx=idx+1
	write(*,"(' Mode',i5,':',E15.5,' Hz',f15.5,' cm-1')") idx,freq(i),wavenum(i)
    freq(i)=freq(i)*factor
    wavenum(i)=wavenum(i)*factor
end do



!Loading finished, now output results
!write(*,"(/,a)") " Note: Only for translation motion, contribution to CV and U are different to CP and H, respectively"
    CV_trans=3D0/2D0*R
    CP_trans=5D0/2D0*R
    Allocate(q_trans(size(T)),U_trans(size(T)),H_trans(size(T)),S_trans(size(T)))
       Do i=1, (size(T))
         q_trans(i)=(2*pi*(totmass*amu2kg)*kb*T(i)/h**2)**(3D0/2D0)*R*T(i)/P
         U_trans(i)=3D0/2D0*R*T(i)
         H_trans(i)=5D0/2D0*R*T(i)
         S_trans(i)=R*(log(q_trans(i)/NA)+5D0/2D0)
       End do

!Convert moment of inertia from a.u.(amu*bohr^2) to kg*m^2
inert11=inert1*amu2kg*(b2a*1D-10)**2
inert21=inert2*amu2kg*(b2a*1D-10)**2
inert31=inert3*amu2kg*(b2a*1D-10)**2
allocate(q_rot(size(T)),U_rot(size(T)),S_rot(size(T)))
      if (nRotsym<0) then !Linear molecule
          CV_rot=R
          Do i=1 , size(T)
          q_rot(i)=8*pi**2*inert1*kb*T(i)/abs(nRotsym)/h**2
	      U_rot(i)=R*T(i)
	      S_rot(i)=R*(log(q_rot(i))+1)
          end do
      else !Non-linear molecule
          CV_rot=3D0*R/2D0
          Do i=1 , size(T)
	      q_rot(i)=8*pi**2/nRotsym/h**3*(2*pi*kb*T(i))**(3D0/2D0)*dsqrt(inert11*inert21*inert31)
	      U_rot(i)=3D0*R*T(i)/2D0
	      S_rot(i)=R*(log(q_rot(i))+3D0/2D0)
          End do
      end if

where (freq<0) freq=0 !Ignore imaginary frequencies
where (wavenum<0) wavenum=0 !Ignore imaginary frequencies
allocate(qvib_v0(size(T)),qvib_bot(size(T)))
   Do j=1 , size(T)
      qvib_v0_heat=1
      qvib_bot_heat=1
      Do i=1,nfreq
         Temp=T(j)
	     if (freq(i)==0) cycle !Ignore imaginary frequency
	     tmpv0=1/( 1-exp( -h*freq(i)/(Kb*Temp) ) )
	     tmpbot=exp(-h*freq(i)/(Kb*2*Temp)) / (1-exp( -h*freq(i)/(Kb*Temp) ))
	     qvib_v0_heat=qvib_v0_heat*tmpv0
	     qvib_bot_heat=qvib_bot_heat*tmpbot
      End do
      qvib_v0(j)=qvib_v0_heat
      qvib_bot(j)=qvib_bot_heat
   End do
   ZPE=h*sum(freq)/2*NA
   allocate(U_vib(size(T)),CV_vib(size(T)),S_vib(size(T)))
   Do j=1 , size(T)
      U_vib_heat=0
      CV_vib_heat=0
      S_vib_heat=0
      Do i=1,nfreq
         Temp=T(j)
	     if (freq(i)==0) cycle !Ignore imaginary frequency
	     prefac=h*freq(i)/(Kb*Temp)
	     term=exp(-h*freq(i)/(Kb*Temp))
         tmp1=R*Temp*prefac*term/(1-term)
	     U_vib_heat=U_vib_heat+tmp1
	     tmp2=R*prefac**2 * term/(1-term)**2
	     CV_vib_heat=CV_vib_heat+tmp2
	     tmp3=R*(prefac*term/(1-term)-log(1-term))
	     S_vib_heat=S_vib_heat+tmp3
      end do
      U_vib(j)=U_vib_heat+ZPE
      CV_vib(j)=CV_vib_heat
      S_vib(j)=S_vib_heat
   End Do

write(*,*)
!write(*,"(a)") " Note: Thermal excitation of electronic states is not taken into account, so electronic contribution to CV and U are zero"

CV_ele=0
U_ele=0
q_ele=nSpinmulti
S_ele=R*log(q_ele)
!write(*,"(' Electronic q: ',f12.6)") q_ele
!write(*,"(' Electronic S: ',f12.6,' cal/mol/K')") S_ele/cal2J

write(*,*)
write(*,*) "****************************************************************"
write(*,*) "*                    Thermodynamic properties                  *"
write(*,*) "****************************************************************"
write(*,*)
write(*,*) " Temperature    S(T)         U(T)        H(T)        G(T)       "
write(*,*) "    K          cal/mol/K    kcal/mol    kcal/mol    kcal/mol    "
Do i=1 ,size(T)
   CV_total=(CV_trans+CV_rot+CV_vib(i))/cal2J
   CP_total=(CP_trans+CV_rot+CV_vib(i))/cal2J
   S_total=(S_trans(i)+S_rot(i)+S_vib(i)+S_ele)/cal2J
   U_total=(U_trans(i)+U_rot(i)+U_vib(i))/cal2J/1000
   H_total=(H_trans(i)+U_rot(i)+U_vib(i))/cal2J/1000
   G_total=H_total-T(i)*S_total/1000
   write(*,"(f12.3,f12.3,f12.3,f12.3,f12.3)") T(i) , S_total , U_total , H_total , G_total
End Do
write(*,*)
write(*,*) " ==============================================================="
write(*,*)
write(*,*) " Temperature     q_v0         q_bot          CV         CP       "
write(*,*) "    K                                     cal/mol/K   cal/mol/K  "
Do i=1 ,size(T)
    q_v0=q_trans(i)*q_rot(i)*qvib_v0(i)*q_ele
    q_bot=q_trans(i)*q_rot(i)*qvib_bot(i)*q_ele
    CV_total=(CV_trans+CV_rot+CV_vib(i))/cal2J
    CP_total=(CP_trans+CV_rot+CV_vib(i))/cal2J
write(*,"(f12.3,E14.6,E14.6,f12.3,f12.3)") T(i) , q_v0 , q_bot , CV_total , CP_total  
End Do

deallocate(freq,wavenum,q_trans,U_trans,H_trans,S_trans)
deallocate(q_rot,U_rot,S_rot)
deallocate(qvib_v0,qvib_bot)
deallocate(U_vib,CV_vib,S_vib)

end do outer

if (isilent==0) then
	write(*,*) "Press ENTER to exit"
	read(*,*)
end if

end program TORSION