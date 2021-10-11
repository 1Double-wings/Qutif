module math_funs
    implicit none
    real(kind=8), parameter :: pi = 3.1415926d0

    !张量积的通用接口
        !tensor_v  :   对矢量张量积
        !tensor_o  :   对算符张量积
    interface tensor
        module procedure :: tensor_v, tensor_o
    end interface

    contains

    subroutine linpace(Ni,Nf,N,samples)
        real(kind=8), intent(in) :: Ni,Nf
        integer, intent(in) :: N
        real(kind=8), allocatable, intent(out) :: samples(:)

        real(kind=8) :: dn
        integer :: i

        allocate(samples(N))
        dn = 1d0/real(N)*(Nf - Ni)

        do i = 1,N
            samples(i) = Ni + real(i-1)*dn
        end do

    end subroutine

    subroutine tensor_o(A,B,C)
        complex(8), dimension(:,:), allocatable, intent(in) :: A,B
        complex(8), dimension(:,:), allocatable, intent(out) :: C

        integer :: i,j,m,n,p,q
        m = size(A,dim=1)
        n = m
        p = size(B,dim=1)
        q = p

        if (.not. allocated(C)) allocate(C(size(A,1)*size(B,1),size(A,1)*size(B,1)))
        do j = 1,n
            do i = 1,m
                C(p*(i-1)+1:p*i, q*(j-1)+1:q*j) = B * A(i,j)
            end do
        end do
    end subroutine

    subroutine tensor_v(A,B,C)
        complex(8), dimension(:), allocatable, intent(in) :: A,B
        complex(8), dimension(:), allocatable, intent(out) :: C

        integer :: i

        if (.not. allocated(C)) allocate(C(size(A) * size(B)))
        do i = 1, size(A)
            C((i-1)*size(B)+1 : i*size(B)) = B * A(i)
        end do
    end subroutine

    subroutine eig_o(A,eig_val,eig_vec)
        complex(kind=8), intent(in) :: A(:,:)
        real(kind=8), intent(out) :: eig_val(:)
        complex(kind=8), intent(out), optional :: eig_vec(:,:)

        complex(kind=8), allocatable :: uptri(:,:), tau(:), work(:), z(:,:), work1test(:), &
            work2(:), work2test(:), test1(:,:), test2(:,:)
        doubleprecision, allocatable :: work3(:), diag(:), ndiag(:)
        integer :: info,info2,lwork,dims

        integer :: i,j


        dims = size(A,1)
        lwork = dims

        allocate(uptri(dims,dims))
        allocate(test1(dims,dims))
        allocate(tau(dims-1))
        allocate(diag(dims))
        allocate(ndiag(dims-1))
        allocate(work1test(2))
        allocate(work2test(2))
        allocate(work3(max(1,2*dims-2)))


        uptri = A
        call zhetrd('U',dims,uptri,dims,diag,ndiag,tau,work1test,-1,info)
        lwork = floor(real(work1test(1)))
        allocate(work(lwork))
        call zhetrd('U',dims,uptri,dims,diag,ndiag,tau,work,lwork,info)
        call zungtr('U',dims,uptri,dims,tau,work2test,-1,info)
        lwork = floor(real(work2test(1)))
        allocate(work2(lwork))
        call zungtr('U',dims,uptri,dims,tau,work2,lwork,info)

        if (present(eig_vec)) then
            eig_vec = uptri
            call zsteqr('V',dims,diag,ndiag,eig_vec,dims,work3,info2)
            eig_val = diag
        else
            allocate(z(1,1))
            call zsteqr('N',dims,diag,ndiag,z,1,work3,info2)
            eig_val = diag
            deallocate(z)
        end if

        deallocate(uptri)
        deallocate(tau)
        deallocate(diag)
        deallocate(ndiag)
        deallocate(work)
        deallocate(work2)
        deallocate(work3)
        deallocate(work1test)
        deallocate(work2test)


    end subroutine

    subroutine eig_gene(A,eig_val,eig_vec)
        complex(kind=8), intent(in) :: A(:,:)
        complex(kind=8), intent(out) :: eig_val(:)
        complex(kind=8), intent(out), optional :: eig_vec(:,:)

        complex(kind=8) :: workf(1), vl(1,1), vr(1,1)
        complex(kind=8), allocatable :: workt(:)
        real(kind=8), allocatable :: rwork(:)
        integer :: info, lworkf=-1, lworkt

        allocate(rwork(2*size(A,1)))

        if (present(eig_vec)) then
            call zgeev("N","V",size(A,1),A,size(A,1),eig_val,eig_vec,size(A,1),eig_vec,size(A,1),workf,lworkf,rwork,info)
            lworkt = floor(real(workt(1)))
            allocate(workt(lworkt))
            call zgeev("N","V",size(A,1),A,size(A,1),eig_val,eig_vec,size(A,1),eig_vec,size(A,1),workt,lworkt,rwork,info)
        else
            call zgeev("N","N",size(A,1),A,size(A,1),eig_val,vl,1,vr,1,workf,lworkf,rwork,info)
            lworkt = floor(real(workt(1)))
            allocate(workt(lworkt))
            call zgeev("N","N",size(A,1),A,size(A,1),eig_val,vl,1,vr,1,workt,lworkt,rwork,info)
        end if

        deallocate(workt)
        deallocate(rwork)
    end subroutine


end module
