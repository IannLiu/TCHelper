Module func
 use lapack95
 implicit real*8 (a-h,o-z)
 contains
 
 subroutine torgana(natm, tbond, nbond, torgroup)
 implicit none
 integer*8,intent(in):: natm, tbond, nbond
 integer*8,intent(inout),allocatable:: torgroup(:,:)
 integer*8,allocatable:: bdarr(:,:), tbondlab(:), bdmar(:,:) ,workmar(:,:)
 integer :: i,j,k,m

 allocate(bdarr(nbond,2), bdmar(natm,natm), tbondlab(tbond), workmar(natm,natm)) 
 call loclabel(10,"Bound Array")
 read(10,*)
 bdmar=0
    do i=1,nbond
        read(10,*) j, bdarr(i,:)
        bdmar(bdarr(i,1),bdarr(i,2))=1
        bdmar(bdarr(i,2),bdarr(i,1))=1
    end do
 call loclabel(10,"Torsion Bound Label")
 read(10,*)
 read(10,*) tbondlab
 torgroup=0
 loop0:do k=1, tbond
    workmar=bdmar
    workmar(bdarr(tbondlab(k),2),:)=0
    workmar(:,bdarr(tbondlab(k),2))=0
    
torgroup(1,(k*2-1))=bdarr(tbondlab(k),1)
torgroup(1,k*2)=bdarr(tbondlab(k),2)
loop1:do m=1, natm
        if (0/=torgroup(m,(k*2-1))) then
            loop2:do i=1,natm
                if (workmar(torgroup(m,(k*2-1)),i)) then
               loop3: do j=1,natm
                             if (0==torgroup(j,(k*2-1))) then
                                   torgroup(j,(k*2-1))=i
                                   exit loop3
                                else if(0/=(torgroup(j,(k*2-1))-i)) then
                                   cycle loop3
                                else
                                   exit loop3
                              end if 
                       end do loop3
                 else 
                       cycle loop2
                 end if
            end do loop2
        else
            exit loop1
        end if
    workmar(torgroup(m,(k*2-1)),:)=0
    workmar(:,torgroup(m,(k*2-1)))=0 
end do loop1
end do loop0

 end subroutine
 
 !!------------ Diagonalize a symmetry matrix 
!Repack the extremely complex "DSYEV" routine in lapack to terse form
!if istat/=0, means error occurs
subroutine diagsymat(mat,eigvecmat,eigvalarr,istat)
integer istat
real*8 mat(:,:),eigvecmat(:,:),eigvalarr(:)
real*8,allocatable :: lworkvec(:)
isize=size(mat,1)
allocate(lworkvec(3*isize-1))
call DSYEV('V','U',isize,mat,isize,eigvalarr,lworkvec,3*isize-1,istat)
eigvecmat=mat
mat=0D0
forall (i=1:isize) mat(i,i)=eigvalarr(i)
end subroutine

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
end subroutine loclabel

!!------ Display matrix similar to gaussian program, automatically switch to next screen
subroutine showmatgau(mat,label,insemi,form,fileid,useri1,useri2,inncol,titlechar)
!Number of columns is always 5, unadjustable
!"Label" is the title, if content is empty, title will not be printed
!If semi==1, only lower and diagonal element will be shown
!"form" is the format to show data, default is D14.6, can pass into such as "f14.8", total width should be 14 characters
!fildid is destination, 6 corresponds output to screen
!useri1 and useri2 is the dimension of the matrix, default or =-1 is determine automatically
!inncol seems controls spacing between number labels of each frame
!titlechar default is "i8,6x", if you manually set inncol, you also set this to broaden or narrow
implicit real*8(a-h,o-z)
real*8 :: mat(:,:)
character(*),optional :: label,form,titlechar
integer,optional :: insemi,fileid,useri1,useri2,inncol
integer :: semi,ides,ncol
semi=0
ides=6
ncol=5
i1=size(mat,1)
i2=size(mat,2)
if (present(useri1).and.useri1/=-1) i1=useri1
if (present(useri2).and.useri1/=-1) i2=useri2
if (present(insemi)) semi=insemi
if (present(fileid)) ides=fileid
if (present(inncol)) ncol=inncol
if (present(label).and.label/='') write(ides,*) "************ ",label," ************"
nt=ceiling(i2/float(ncol))
do i=1,nt !How many frames
	ns=(i-1)*5+1 !This frame starts from where
	if (i/=nt) ne=(i-1)*ncol+ncol !This frame end to where
	if (i==nt) ne=i2
	!Write basis number in separate line
	write(ides,"(6x)",advance='no')
	do j=ns,ne
		if (present(titlechar)) then
			write(ides,'('//titlechar//')',advance='no') j
		else
			write(ides,"(i8,6x)",advance='no') j
		end if
	end do
	write(ides,*)
	!Write content in each regular line
	do k=1,i1
		if (k<ns.and.semi==1) cycle !The lines have been outputted are skipped
		write(ides,"(i6)",advance='no') k
		do j=ns,ne
			if (semi==1.and.k<j) cycle !Upper trigonal element were passed
			if (present(form)) then
				write(ides,'('//form//')',advance='no') mat(k,j)			
			else
				write(ides,"(D14.6)",advance='no') mat(k,j)
			end if
		end do
		write(ides,*) !Change to next line
	end do
end do
end subroutine showmatgau

 end module func