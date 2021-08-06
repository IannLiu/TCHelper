Module TORANA
use func
contains

subroutine ANALYSIS(natm,atmass,xyz,ntor,wavall,wav_proall,inert1,inert2,inert3)
implicit real*8 (a-h,o-z)
real*8,parameter :: pi=3.141592653589793D0
integer*8,intent(inout)::natm, ntor
real*8, intent(inout)::atmass(:)
real*8, intent(inout)::inert1,inert2,inert3
real*8,allocatable,intent(inout)::wavall(:), wav_proall(:),xyz(:,:)
real*8,allocatable :: Hess(:,:),eigvecmat(:,:),eigvalarr(:),kmat(:,:),massmat(:,:),freq(:), fre_pro(:), normvec(:,:),tmpvec(:)
real*8,allocatable ::D(:,:) ,E(:,:),P(:,:) ,origin(:) 
integer*8, allocatable :: torgroup(:,:)
real*8 :: VnV(3),V1V2(3)
real*8 :: MIM(3,3) , Meigvecmat(3,3) , Meigvalarr(3)
integer*8 :: nbond,existence, nmode, nmodeall


!Don't support linear molecules
write(*,*) "*****************************************************"
write(*,*) "*          Torsion Projection analysis              * "
write(*,*) "*****************************************************"
write(*,*)


nmode=3*natm-6 !The number of vibrational modes
nmodeall=3*natm
write(*,"(' The number of vibrational modes:',i6)") nmode


write(*,*) "Atomic masses:"
write(*,"(6(f12.6))") atmass

!load xyz and fixe it in the center of mass
xyz=xyz/0.52917720859D0
allocate(origin(3))
write(*,*) "Cartesian coordinate(Bohr)"
do i=1, natm
write(*,"(3(f12.6))") xyz(i,:)
end do
origin=0
do i=1, natm
    origin(1)=atmass(i)*xyz(i,1)+origin(1)
    origin(2)=atmass(i)*xyz(i,2)+origin(2)
    origin(3)=atmass(i)*xyz(i,3)+origin(3)
end do
temp=sum(atmass)
origin=origin/temp
write(*,*) "Origin cartesian coordinate"
write(*,*) origin
do i=1 , natm
    xyz(i,1)=xyz(i,1)-origin(1)
    xyz(i,2)=xyz(i,2)-origin(2)
    xyz(i,3)=xyz(i,3)-origin(3)
end do
write(*,*) "Cartesian coordinate fixed at mass center"
do i=1,natm
write(*,"(3(f12.6))") xyz(i,:)
end do

!Load Hessian matrix 
allocate(Hess(nmodeall,nmodeall))
Hess=0
call loclabel(10,"Cartesian Force Constants")
read(10,*)
read(10,"(5(1PE16.8))") ((Hess(i,j),j=1,i),i=1,nmodeall)
Hess=Hess+transpose(Hess)
do i=1,nmodeall
	Hess(i,i)=Hess(i,i)/2D0
end do
! call showmatgau(Hess,"Hessian matrix",1)

!Load bond array and torsion bond, then analysis torsion group
call loclabel(10,"Total Bond Number")
read(10,*)
read(10,*) nbond
allocate(torgroup(natm,2*ntor))
if (ntor>0) then
call torgana(natm, ntor, nbond, torgroup)
 write(*,*)
 write(*,*) "********************Torsion bond group*****************"
 do i=1,ntor
     write(*,"(' Pivot atom lables where tosion is happening:'I8,I8)") torgroup(1,(i*2-1)),torgroup(1,i*2)
     write(*,*) "One of the torsion groups about mentioned torsion above:"
     write(*,"(10(I5))") torgroup(:,(i*2-1))
     write(*,*)
 end do
else 
    write(*,*) "Warning!!!  Torsion haven't been detected!!!"
end if


!Construct mass matrix
allocate(massmat(nmodeall,nmodeall))
massmat=0
do i=1,natm
	do j=(i-1)*3+1,i*3
		massmat(j,j)=atmass(i)
    end do
end do

!Calculat the moments of inertia 
MIM=0
do iatm=1,natm
    MIM(1,1)=atmass(iatm)*(xyz(iatm,2)**2+xyz(iatm,3)**2)+MIM(1,1)
    MIM(1,2)=-atmass(iatm)*(xyz(iatm,1)*xyz(iatm,2))+MIM(1,2)
    MIM(1,3)=-atmass(iatm)*(xyz(iatm,1)*xyz(iatm,3))+MIM(1,3)
    MIM(2,2)=atmass(iatm)*(xyz(iatm,1)**2+xyz(iatm,3)**2)+MIM(2,2)
    MIM(2,3)=-atmass(iatm)*(xyz(iatm,2)*xyz(iatm,3))+MIM(2,3)
    MIM(3,3)=atmass(iatm)*(xyz(iatm,1)**2+xyz(iatm,2)**2)+MIM(3,3)
end do
MIM(2,1)=MIM(1,2)
MIM(3,1)=MIM(1,3)
MIM(3,2)=MIM(2,3)
write(*,*)
write(*,*) "********************************************************"
write(*,*) "*     Now, output Properties of global rotation        *"
write(*,*) "********************************************************"
write(*,*)
write(*,*) "The moment of inertia(amu*Bohr^2) martix is:"
do i=1,3
write(* , "(f12.6, f12.6, f12.6)") MIM(i,1) , MIM(i,2) , MIM(i,3)
end do
call diagsymat(MIM,Meigvecmat,Meigvalarr,Mistat) !Note the column of eigvecmat has already been automatically normalized
if (istat==0) write(*,*) "Diagonalization passed"
write(*,*) "Principle moment inertia is:"
write(*,"(f12.6,f12.6,f12.6)") Meigvalarr
inert1=Meigvalarr(1)
inert2=Meigvalarr(2)
inert3=Meigvalarr(3)
write(*,*) "And corresponding eigenvect is:"
do i=1,3
    write(*,"(f12.6,f12.6,f12.6)") Meigvecmat(i,:)
end do


!construct translation vector and rotation vector
allocate(D(nmodeall,(6+ntor)),E(nmodeall,nmodeall),P(nmodeall,nmodeall))
D=0
do i=1,natm !global translation vector
    D((i-1)*3+1,1)=1
    D((i-1)*3+2,2)=1
    D((i-1)*3+3,3)=1
end do
do i=1, natm !global rotation vector
    do j=1,3
    D((i-1)*3+j,4)=dot_product(xyz(i,:),Meigvecmat(2,:))*Meigvecmat(j,3)-dot_product(xyz(i,:),Meigvecmat(3,:))*Meigvecmat(j,2)
    D((i-1)*3+j,5)=dot_product(xyz(i,:),Meigvecmat(3,:))*Meigvecmat(j,1)-dot_product(xyz(i,:),Meigvecmat(1,:))*Meigvecmat(j,3)
    D((i-1)*3+j,6)=dot_product(xyz(i,:),Meigvecmat(1,:))*Meigvecmat(j,2)-dot_product(xyz(i,:),Meigvecmat(2,:))*Meigvecmat(j,1)
    end do
end do

do j=7,(6+ntor)
  do i=1,natm !internal rotation vector
     inner:do k=1, natm
         if (i==torgroup(k,((j-6)*2-1))) then
           existence=1
           exit inner
         else if(0==torgroup(k,((j-6)*2-1))) then
           existence=0
           exit inner
         else 
           existence=0
         end if
         end do inner
    if (existence==1)   then
        VnV=xyz(i,:)-xyz(torgroup(1,((j-6)*2-1)),:)
        V1V2=xyz(torgroup(1,((j-6)*2-1)),:)-xyz(torgroup(1,((j-6)*2)),:)
        D((i-1)*3+1,j)=(VnV(2)*V1V2(3)-VnV(3)*V1V2(2))/dsqrt(V1V2(1)**2+V1V2(2)**2+V1V2(3)**2)
        D((i-1)*3+2,j)=(VnV(3)*V1V2(1)-VnV(1)*V1V2(3))/dsqrt(V1V2(1)**2+V1V2(2)**2+V1V2(3)**2)
        D((i-1)*3+3,j)=(VnV(1)*V1V2(2)-VnV(2)*V1V2(1))/dsqrt(V1V2(1)**2+V1V2(2)**2+V1V2(3)**2)
    else
        VnV=xyz(i,:)-xyz(torgroup(1,((j-6)*2)),:)
        V1V2=xyz(torgroup(1,((j-6)*2)),:)-xyz(torgroup(1,((j-6)*2-1)),:)
        D((i-1)*3+1,j)=(VnV(2)*V1V2(3)-VnV(3)*V1V2(2))/dsqrt(V1V2(1)**2+V1V2(2)**2+V1V2(3)**2)
        D((i-1)*3+2,j)=(VnV(3)*V1V2(1)-VnV(1)*V1V2(3))/dsqrt(V1V2(1)**2+V1V2(2)**2+V1V2(3)**2)
        D((i-1)*3+3,j)=(VnV(1)*V1V2(2)-VnV(2)*V1V2(1))/dsqrt(V1V2(1)**2+V1V2(2)**2+V1V2(3)**2)
    end if
  end do
end do

do i=2,(6+ntor) !Schmidt orthogonalization
    do j=1,i-1
    D(:,i)=D(:,i)-(dot_product(D(:,j),D(:,i))/dot_product(D(:,j),D(:,j)))*D(:,j)
    end do
end do

do i=1,(6+ntor)
tmp=dsqrt(sum(D(:,i)**2))
D(:,i)=D(:,i)/tmp
end do
E=0
do i=1,nmodeall
    E(i,i)=1
end do

!Now, calculating frequencies
!===============================================================
!========Calculate original frequenceis at first================
!Construct force constant matrix (mass-weighted Hessian matrix)
write(*,*)
write(*,*)"**************************************************************"
write(*,*)"* Now,frequencies are calculating without torsion projection *"
write(*,*)"*          Note: Normal coordinate will be print             *"
write(*,*)"**************************************************************"
write(*,*)
allocate(kmat(nmodeall,nmodeall))
do i=1,nmodeall
	do j=1,nmodeall
		kmat(i,j)=hess(i,j)/dsqrt(massmat(i,i)*massmat(j,j))
	end do
end do
write(*,*)
call showmatgau(kmat,"Force constant matrix (i.e mass-weighted Hessian)")

!Diagonalization of force constant matrix
allocate(eigvecmat(nmodeall,nmodeall),eigvalarr(nmodeall))
call diagsymat(kmat,eigvecmat,eigvalarr,istat) !Note the column of eigvecmat has already been automatically normalized
if (istat==0) write(*,*) "Diagonalization passed"
write(*,*) eigvalarr

!Convert force constant to harmonic frequencies
amu2kg=1.66053878D-27
b2m=0.529177249D-10
au2J=4.35974434D-18
eigvalarr=eigvalarr*au2J/b2m**2/amu2kg !First convert force constant from a.u. to SI
allocate(freq(nmodeall))
do i=1,nmodeall
	if (eigvalarr(i)<0) then
		freq(i)=-dsqrt(abs(eigvalarr(i)))/(2*pi)
	else
		freq(i)=dsqrt(eigvalarr(i))/(2*pi)
	end if
end do

!Sort all modes according to absolute value of frequencies (from low to high)
allocate(tmpvec(nmodeall))
do i=1,nmodeall
	do j=i+1,nmodeall
		if (abs(freq(i))>abs(freq(j))) then
			temp=freq(i)
			freq(i)=freq(j)
			freq(j)=temp
			tmpvec=eigvecmat(:,i)
			eigvecmat(:,i)=eigvecmat(:,j)
			eigvecmat(:,j)=tmpvec
		end if
	end do
end do
!Now the first six modes are considered as overall motions, we now sort remaining modes from low to high
do i=7,nmodeall
	do j=i+1,nmodeall
		if (freq(i)>freq(j)) then
			temp=freq(i)
			freq(i)=freq(j)
			freq(j)=temp
			tmpvec=eigvecmat(:,i)
			eigvecmat(:,i)=eigvecmat(:,j)
			eigvecmat(:,j)=tmpvec
		end if
	end do
end do
wavall=freq/2.99792458D10

!Convert normal coordinates from mass-weighted basis to Cartesian basis
allocate(normvec(nmodeall,nmodeall))
do i=1,natm
	do j=(i-1)*3+1,i*3
		normvec(j,:)=eigvecmat(j,:)/dsqrt(atmass(i))
	end do
end do
do i=1,nmodeall !Normalization
	tmp=dsqrt(sum(normvec(:,i)**2))
	normvec(:,i)=normvec(:,i)/tmp
end do
write(*,*)
call showmatgau(normvec(:,7:),"Normal coordinates (columns)",form="1x,f9.4,4x")

!Output final frequencies
write(*,*)
write(*,*) "The frequencies (cm-1) corresponding to overall translation and rotation:"
write(*,"(6(f12.5))") wavall(1:6)
write(*,*) "Harmonic vibrational frequencies:"
idx=0
do i=7,nmodeall
	idx=idx+1
	write(*,"(' Mode',i5,':',E15.5,' Hz',f15.5,' cm-1')") idx,freq(i),wavall(i)
end do
deallocate(kmat, eigvecmat,  eigvalarr, tmpvec, normvec)
write(*,*)
write(*,*)"**************************************************************"
write(*,*)"*  Now,frequencies are calculating with torsion projection   *"
write(*,*)"*      Note: Normal coordinate will not be print             *"
write(*,*)"**************************************************************"
write(*,*)
P=E-matmul(D,transpose(D))
Hess=matmul(P,Hess)
Hess=matmul(Hess,P)

!Construct force constant matrix (mass-weighted Hessian matrix)
allocate(kmat(nmodeall,nmodeall))
do i=1,nmodeall
	do j=1,nmodeall
		kmat(i,j)=hess(i,j)/dsqrt(massmat(i,i)*massmat(j,j))
	end do
end do
write(*,*)
call showmatgau(kmat,"Force constant matrix (i.e mass-weighted Hessian)")

!Diagonalization of force constant matrix
allocate(eigvecmat(nmodeall,nmodeall),eigvalarr(nmodeall))
call diagsymat(kmat,eigvecmat,eigvalarr,istat) !Note the column of eigvecmat has already been automatically normalized
if (istat==0) write(*,*) "Diagonalization passed"
write(*,*) eigvalarr

!Convert force constant to harmonic frequencies
amu2kg=1.66053878D-27
b2m=0.529177249D-10
au2J=4.35974434D-18
eigvalarr=eigvalarr*au2J/b2m**2/amu2kg !First convert force constant from a.u. to SI
allocate(fre_pro(nmodeall))
do i=1,nmodeall
	if (eigvalarr(i)<0) then
		fre_pro(i)=-dsqrt(abs(eigvalarr(i)))/(2*pi)
	else
		fre_pro(i)=dsqrt(eigvalarr(i))/(2*pi)
	end if
end do

!Sort all modes according to absolute value of frequencies (from low to high)
allocate(tmpvec(nmodeall))
do i=1,nmodeall
	do j=i+1,nmodeall
		if (abs(fre_pro(i))>abs(fre_pro(j))) then
			temp=fre_pro(i)
			fre_pro(i)=fre_pro(j)
			fre_pro(j)=temp
			tmpvec=eigvecmat(:,i)
			eigvecmat(:,i)=eigvecmat(:,j)
			eigvecmat(:,j)=tmpvec
		end if
	end do
end do
!Now the first six+ntor modes are considered as overall motions, we now sort remaining modes from low to high
do i=7+ntor,nmodeall
	do j=i+1,nmodeall
		if (fre_pro(i)>fre_pro(j)) then
			temp=fre_pro(i)
			fre_pro(i)=fre_pro(j)
			fre_pro(j)=temp
			tmpvec=eigvecmat(:,i)
			eigvecmat(:,i)=eigvecmat(:,j)
			eigvecmat(:,j)=tmpvec
		end if
	end do
end do
wav_proall=fre_pro/2.99792458D10

!Output final frequencies
write(*,*)
    write(*,*) "The frequencies (cm-1) corresponding to overall translation and rotation:"
    write(*,"(6(f12.5))") wav_proall(1:6)
    write(*,*) "Harmonic vibrational frequencies:"
    idx=0
    do i=7,nmodeall
	    idx=idx+1
	    write(*,"(' Mode',i5,':',E15.5,' Hz',f15.5,' cm-1')") idx,freq(i),wav_proall(i)
    end do
    write(*,*)
    end subroutine ANALYSIS

    end module TORANA