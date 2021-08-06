program main
    use func
    implicit real*8 (a-h,o-z)
    real*8,parameter :: NA=6.02214179D23 !Avogadro constant
    real*8,parameter :: h=6.62606896D-34 !Planck constant, in J*s
    real*8,parameter :: wave2freq=2.99792458D10 !cm^-1 to s^-1 (Hz)
    real*8,parameter :: au2kcal_mol=627.51D0,au2KJ_mol=2625.5D0,au2J=4.359744575D-18,au2cm_1=219474.6363D0 !Hartree to various units
    real*8,parameter :: R=8.3144648D0, Kb=1.3806503D-23 !Ideal gas constant (J/mol/K), Boltzmann constant (J/K)
    character :: File_name*80 , c80*80
    real*8,allocatable :: xyz(:,:) , mass(:)
    real*8,allocatable :: angle(:), X(:), E(:) ,EV(:) ,S(:),Parfun(:),U(:),CV(:)
    integer:: bdnum, ith ,nang ,fter,stp,ISN
    integer, allocatable:: bdarr(:,:),bdmar(:,:),atnum(:)
    real*8 ::  Axvec(3), CML(3),CMR(3),Temperature(2)
    real*8 ::  Itot
    character*200 filename
    character,allocatable :: elem(:)*2
    logical :: alive
    
write(*,*)"===================Description================================"
write(*,*)"Program name: TORTHMAL"
write(*,*)"Programmed by YanLiu (IanPioneee@163.com)"
write(*,*) "Release date: 2020-March-18"
write(*,*)"Function: Input a formated file inclouding torsion information,"
write(*,*)"and calculate the thermo properties of each internal rotation"
write(*,*)

supouter: do while (.true.)

    isilent=0
call getarg(1,filename)
if (trim(filename)=="") then
	write(*,*) "Input path of a formated file"
	do while(.true.)
		read(*,"(a)") filename
		inquire(file=filename,exist=alive)
		if (alive) exit
		write(*,*) "Cannot find the file, input again"
	end do
else
	isilent=1
end if

open(10,file=filename,status="old")
!load atom xyz
    call loclabel(10,"Atom Number")
    read(10 ,*)
    read(10,*) natm  !The number of atoms
    write(*,"(' The number of atoms:',i5)") natm
    !load atom xyz
    allocate(xyz(natm,3),elem(natm))
    call loclabel(10,"Atom Coordinate")
    read(10,*)
    do iatm=1,natm
	    read(10,*) elem(iatm),xyz(iatm,:)
    end do
!loading xyz end


!load bond information
    call loclabel(10,"Total Bond Amount")
    read(10,*)
    read(10,*) bdnum
    call loclabel(10,"Bond Array")
    read(10,*)
    read(10,*)
    allocate(bdarr(bdnum,2),bdmar(natm,natm))
    bdmar=0
    do i=1,bdnum
        read(10,*) j, bdarr(i,:)
        bdmar(bdarr(i,1),bdarr(i,2))=1
        bdmar(bdarr(i,2),bdarr(i,1))=1
    end do
call loclabel(10,"Rotation Bound Label")
read(10,*)
read(10,*) ith

!Determine masses
call loclabel(10,"Mass For All Atoms")
read(10,*)
allocate(mass(natm))
mass=0 !allocate initial data,zero,to all atoms
outer:do 
	read(10,*,iostat=ierror) c80,tmpmass
	if (ierror/=0) exit
	if (iachar(c80(1:1))>=48.and.iachar(c80(1:1))<=57) then
		read(c80,*) iatm
		mass(iatm)=tmpmass
	else
		do iatm=1,natm
			if (c80(1:2)==elem(iatm)) mass(iatm)=tmpmass
        end do
	end if
end do outer
if (any(mass==0)) write(*,*) "Warning: Mass of some atoms are not specified!"

totmass=sum(mass(:))
write(*,"(' Total mass:',f12.6,' amu')") totmass
write(*,*) "the input coordinates of all atom and its mass is:"
      do iatm=1,natm
          write(* , "(A5 , f12.6, f12.6, f12.6, f12.6)")elem(iatm) ,xyz(iatm,:) , mass(iatm)
      end do
  
!load potential energy
call loclabel(10,"Discrete Potential Number")
read(10,*)
read(10,*) nang
call loclabel(10,"Number Of Fourier Expansion Terms")
read(10,*)
read(10,*) fter

allocate(angle(nang),E(nang))
call loclabel(10,"Aangle And Energy")
read(10,*)
    do i=1, nang
        read(10,*) angle(i), E(i)
    end do
    
!load FGH information
call loclabel(10,"Numbers Of FGH")
    read(10,*)
    read(10,*) N  !Gride numbers
call loclabel(10,"Solution Area")
    read(10,*)
    read(10,*) L  !Solution arear length
call loclabel(10,"Numbers Of Energy Level")
    read(10,*)
    read(10,*) EL !the number of output energy levels
    
!Load thermal property information
call loclabel(10,"Temperature Span")
read(10,*)
read(10,*) Temperature(1),Temperature(2) , stp
call loclabel(10,"Internal Symmetry Number")
read(10,*)
read(10,*) ISN


close(10)

!Rotation group analysis
call rotana(bdarr,bdmar , natm , ith , atnum)

!Calculate the center of mass of rotation groups
CML=0
GmassL=0
CMR=0
GmassR=0
outer1:do i=1,natm
    inner:do j=1, natm
    if (i==atnum(j)) then
        existence=1
        exit inner
    else if(0==atnum(j)) then
        existence=0
        exit inner
    else
        existence=0
    end if
    end do inner
    if (existence==1) then
    CML(1)=mass(i)*xyz(i,1)+CML(1)
    CML(2)=mass(i)*xyz(i,2)+CML(2)
    CML(3)=mass(i)*xyz(i,3)+CML(3)
    GmassL=mass(i)+GmassL
    else 
    CMR(1)=mass(i)*xyz(i,1)+CMR(1)
    CMR(2)=mass(i)*xyz(i,2)+CMR(2)
    CMR(3)=mass(i)*xyz(i,3)+CMR(3)
    GmassR=mass(i)+GmassR
    end if
end do outer1
write(*,"(f12.6 , f12.6 , f12.6)") GmassL
write(*,"(f12.6 , f12.6 , f12.6)") GmassR
CML=CML/GmassL
CMR=CMR/GmassR
Axvec=CML-CMR
write(*,*) "The coordinates of center of mass of rotation group is"
write(*,"(f12.6 , f12.6 , f12.6)") CML
write(*,"(f12.6 , f12.6 , f12.6)") CMR

!Coordinate transformation for center of mass
write(*,*) "The new coordinate(origin in the center of gravity)"
do iatm=1 , natm
    xyz(iatm , :)=xyz(iatm , :)-CML
    write(*,"(f12.6 , f12.6 , f12.6 )") xyz(iatm , :)
end do

!Coordinate orientation tansformation
call xyzchange(xyz ,Axvec , natm)
write(*,*)"Coordnate tansformation"
do iatm=1,natm
write(* , "(f12.6, f12.6, f12.6)") xyz(iatm,:)
end do

!Calculating reduced moment of inertia
call rmoi(mass, xyz ,atnum , Itot ,natm)
write(*,*) "Reduced moment of inertia is :"
write(*,"(3(f12.6))") Itot

!calculate Fourier expansion term
call FEC(angle , E , nang , fter ,X )
write(*,*)"The form of Fourier expansion terms is"
write(*,*)"V(theta)=Ai*(1-cos(i*theta))+Bi*sin(i*theta)"
write(*,*)"The values of Ai and Bi is"
write(*,*)"         i           A           B"
do i=1, fter
write(*,"(I12,f12.6,f12.6)") i, X(i*2-1) , X(i*2)  !Write Fourier expansion terms coefficient
end do

!Calculate energy levels
!Need input:Fourier expansion terms coefficient, Gride number,
write(*,*) "Energy levels are been calculating"
write(*,*) "Please wait for a moment... ..."
call FGH(N,L,Itot,X,EV)
    write(*,*)"Energy level is:"
    do i=1 , EL
    write(*,*) EV(i)
    end do
    
!!Calculate the thermal property
allocate(Parfun(stp+1),S(stp+1),U(stp+1),CV(stp+1))
loop1:do i=1,(stp+1)
    q=0
    n=0
    T=Temperature(1)+(i-1)*Temperature(2)
    loop2:do j=1, size(EV)
        if(exp(-EV(j)*au2KJ_mol*1000/Kb/T/NA)>=0.000001) then
        q=exp(-EV(j)*au2KJ_mol*1000/Kb/T/NA)+q
        n=n+1
        else
            exit loop2
        end if
        end do loop2
        Parfun(i)=q/ISN
        sup=0
        sdown=0
        CVup1=0
        CVup2=0
        loop3:do k=1,n
         sup=(EV(k)*au2KJ_mol*1000/NA)*exp(-EV(k)*au2KJ_mol*1000/Kb/T/NA)+sup
         sdown=exp(-EV(k)*au2KJ_mol*1000/Kb/T/NA)+sdown
         CVup1=(EV(k)*au2KJ_mol*1000/NA)**2/(Kb*T*T)*exp(-EV(k)*au2KJ_mol*1000/Kb/T/NA)+CVup1
         CVup2=(EV(k)*au2KJ_mol*1000/NA)/(Kb*T*T)*exp(-EV(k)*au2KJ_mol*1000/Kb/T/NA)+CVup2
        end do loop3
        S(i)=(Kb*log(Parfun(i))+sup/sdown/T)*NA/4.184
        U(i)=(sup/sdown)*NA/1000/4.184
        CV(i)=(sdown*CVup1-sup*CVup2)/(sdown**2)*NA/4.184
end do loop1
write(*,*)"Thermal properties:"
write(*,*)"     T(K)          Q      S(cal/K/mol)    U(kcal/mol)   Cv(cal/K/mol)"
do i=1,(stp+1)
write(*,"(f10.3,4(f14.6))")Temperature(1)+(i-1)*Temperature(2), Parfun(i),S(i),U(i),CV(i)
end do     
    
    if (isilent==0) then
    write(*,*) "Input 0 or 1"
    write(*,*)"0: Exit, 1:Continue"
    read(*,*) comm
    end if
    if (comm==0) exit supouter     
    deallocate(xyz, mass)
    deallocate(angle, X, E, EV, S, Parfun, U, CV)
    deallocate(bdarr, bdmar, atnum)
    deallocate(elem)
    
end do supouter

write(*,*) "Press ENTER to exit"
read(*,*)

end program
