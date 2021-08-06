Module func
 use lapack95
 implicit none
 real*8,parameter :: pi=3.141592653589793D0
    contains
    
    
!!!-------- Find the line where the label first appears in fileid
!Return ifound=1 if found the label, else return 0
!Default is rewind, if irewind=0 will not rewind
!If current line already has the label, calling this subroutine will do nothing
subroutine loclabel(fileid,label,ifound,irewind)
integer fileid,ierror
integer,optional :: ifound,irewind
character*200 c200
CHARACTER(LEN=*) label
if ((.not.present(irewind)).or.(present(irewind).and.irewind==1)) rewind(fileid)
do while(.true.)
	read(fileid,"(a)",iostat=ierror) c200
	if (index(c200,label)/=0) then
		backspace(fileid)
		if (present(ifound)) ifound=1 !Found result
		return
	end if
	if (ierror/=0) exit
end do
if (present(ifound)) ifound=0
end subroutine

 
 subroutine rotana(bdarr, bdmar, natm , ith , atnum)
 implicit none
 integer,intent(inout):: bdmar(:,:),bdarr(:,:)
 integer,intent(inout):: natm , ith
 integer,intent(out),allocatable:: atnum(:)
 integer :: m,i,j

 allocate(atnum(natm))
    atnum=0
    bdmar(bdarr(ith,2),:)=0
    bdmar(:,bdarr(ith,2))=0
    
atnum(1)=bdarr(ith,1)
supoutor:do m=1, natm
        if (0/=atnum(m)) then
            outer1:do i=1,natm
                if (bdmar(atnum(m),i)) then
               inner1: do j=1,natm
                             if (0==atnum(j)) then
                                   atnum(j)=i
                                   exit inner1
                                else if(0/=(atnum(j)-i)) then
                                   cycle inner1
                                else
                                   exit inner1
                              end if 
                       end do inner1
                 else 
                       cycle outer1
                 end if
            end do outer1
        else
            exit supoutor
        end if
    bdmar(atnum(m),:)=0
    bdmar(:,atnum(m))=0 
end do supoutor 
write(*,*)"========Rotation atoms number======="
write(*,*) atnum
write(*,*)"===================================="
 end subroutine  rotana
    
subroutine rmoi(mass, xyz ,atnum , Itot ,natm)
implicit none
real*8,intent(out) ::Itot
real*8,intent(in) ::mass(:),xyz(:,:)
integer,intent(in) ::atnum(:) ,natm
real*8 :: tot, IL, IR
integer ::iatm ,i

tot=0
IL=0 
do iatm=1, natm
    tot=(mass(iatm)*(xyz(iatm,1)**2+xyz(iatm,2)**2)+tot)
end do
do i=1, natm 
    if (0/=atnum(i)) then
        IL=mass(atnum(i))*(xyz(atnum(i),1)**2+xyz(atnum(i),2)**2)+IL
    else 
        exit
    end if
end do
IR=tot-IL
Itot=1/((1/IL)+(1/IR))

end subroutine rmoi

subroutine xyzchange(xyz ,Axvec , natm)
implicit none
real*8, intent(inout)::xyz(:,:)
real*8, intent(in)::Axvec(:)
integer,intent(in) ::natm
integer :: iatm
real*8 :: Lvec , costheta , cosfei , sintheta , sinfei

if (Axvec(1)**2+Axvec(2)**2/=0) then
    cosfei=Axvec(1)/sqrt(Axvec(1)**2+Axvec(2)**2)
    sinfei=Axvec(2)/sqrt(Axvec(1)**2+Axvec(2)**2)
   do iatm=1 , natm
       if (xyz(iatm,1)**2+xyz(iatm,2)**2/=0) then
           Lvec=sqrt(xyz(iatm,1)**2+xyz(iatm,2)**2)
           costheta=xyz(iatm,2)/Lvec
           sintheta=xyz(iatm,1)/Lvec
           xyz(iatm,1)=Lvec*(sinfei*costheta+cosfei*sintheta)
           xyz(iatm,2)=Lvec*(cosfei*costheta-sinfei*sintheta)
       else
           xyz(iatm,1)=xyz(iatm,1)
           xyz(iatm,2)=xyz(iatm,2)
       end if
   end do
end if
if (Axvec(1)**2+Axvec(2)**2+Axvec(3)**2/=0) then
    cosfei=Axvec(3)/sqrt(Axvec(1)**2+Axvec(2)**2+Axvec(3)**2)
    sinfei=sqrt(Axvec(1)**2+Axvec(2)**2)/sqrt(Axvec(1)**2+Axvec(2)**2+Axvec(3)**2)
    do iatm=1 , natm
       if (xyz(iatm,1)**2+xyz(iatm,3)**2/=0) then
           Lvec=sqrt(xyz(iatm,1)**2+xyz(iatm,3)**2)
           costheta=xyz(iatm,1)/Lvec
           sintheta=xyz(iatm,3)/Lvec
           xyz(iatm,3)=Lvec*(sintheta*cosfei+costheta*sinfei)
           xyz(iatm,1)=Lvec*(costheta*cosfei-sintheta*sinfei)
       else
           xyz(iatm,1)=xyz(iatm,1)
           xyz(iatm,3)=xyz(iatm,3)
       end if
    end do
else
    write(*,*)"Warnning, error totation axle.Coordnate tansformation failed"
end if
end subroutine xyzchange

Subroutine FEC(ang, Eng , anum , fter ,X )
    implicit None
    real*8,intent(in),allocatable :: ang(:), Eng(:)
    integer,intent(in) :: anum , fter
    real*8,intent(out), allocatable :: X(:)
    real*8,allocatable :: A(: , :), AT_A(:,:),AT_E(:) 
    integer,allocatable:: ipiv(:)
    integer :: i , j
    real*8 ::theta
    
    allocate(ipiv(2*fter),A(anum,2*fter),AT_A(2*fter,2*fter),AT_E(2*fter),X(2*fter))
    do i=1,anum
        theta=ang(i)/180.0*pi
        do j=1, fter
            A(i,2*j-1)=1-cos(j*theta)
            A(i,2*j)=sin(j*theta)
        end do
    end do
    AT_A=matmul(transpose(A),A)
    AT_E=matmul(transpose(A),Eng)
    call getrf(AT_A, ipiv )
    call getri(AT_A, ipiv )
    X=matmul(AT_A,AT_E) 
end Subroutine FEC

Subroutine FGH(N,L,Itot,X,EV)
integer,intent(in) :: N, L
real*8,intent(in), allocatable::X(:)
real*8,intent(in)::Itot
real*8,intent(out), allocatable:: EV(:)
integer :: i ,j ,k ,m,ifo,Length
real*8, allocatable :: H(:,:), wi(:), vr(:,:), vl(:,:)
real*8::delx , delk, T , Eng ,V

allocate(H(N,N),EV(N),vr(N,N))
allocate(wi(N),vl(N,N))
Length=size(X)/2
    delx=(L*pi)/N
    delk=2*pi/(L*pi)
  do i=1, N
        do j=1, N
            T=0
            do k=1,(N-1)/2
                T=cos(k*2*pi*(i-j)/N)*(2/(Itot*1822/(0.52917720859**2))*(pi*k/(N*delx))**2)+T
            end do
            if(i/=j) then
            H(i,j)=2.0/N*T
            else
			V=0
			do m=1, Length
			V=V+X(m*2-1)*(1-cos(m*i*delx))+X(m*2)*sin(m*i*delx)
			end do
            H(i,j)=2.0/N*T+V
            end if
        end do
    end do
    call geev(H , EV, wi ,vl, vr ,ifo)
    if (ifo==0) then
        write(*,*) "Energy levels were successfully obtained"
    else
        write(*,*) "Calculation process of Energy levels failed"
    end if
    !Bubble Sort
    outer: do j=1 ,N-1
        inner:do i=1 , N-1
         if (EV(i)>EV(i+1)) then
             Eng=EV(i)
             EV(i)=EV(i+1)
             EV(i+1)=Eng         
         else
             cycle inner
         end if
        end do inner
    end do outer
end Subroutine FGH
 
end module func