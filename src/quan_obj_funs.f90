module quan_obj_funs
    !����������Լ�ʸ���ϵĺ���
    use quan_obj
    use math_funs

    implicit none

    !���׹����ͨ�ýӿ�
        !dag_v   :   ��ʸ��dag
        !dag_o   :   �����dag
    interface dag
        module procedure :: dag_v, dag_o
    end interface dag

    interface eig
        module procedure :: eig_hermitian, eig_non_hermitian
    end interface


    contains

    function dag_v(psi)
        !���������
            !psi   :   һ��ʸ��ʵ��psi
        !���������
            !dag_v :   psi^{\dagger}
        type(vecs), intent(in) :: psi
        type(vecs) :: dag_v

        dag_v%lengs = psi%lengs
        allocate(dag_v%subs(size(psi%subs)))
        dag_v%subs = psi%subs

        if (psi%k_b .eq. 'ket') then
            dag_v%k_b = 'bra'
        else
            dag_v%k_b = 'ket'
        end if

        allocate(dag_v%elem(psi%lengs))
        dag_v%elem = conjg(psi%elem)

    end function

    function dag_o(ope)
        !���������
            !ope   :   һ�����ʵ��ope
        !���������
            !dag_o :   o^{\dagger}
        type(opera), intent(in) :: ope
        type(opera) :: dag_o

        dag_o%dims = ope%dims
        allocate(dag_o%subs(size(ope%subs)))
        dag_o%subs = ope%subs

        allocate(dag_o%elem(ope%dims,ope%dims))
        dag_o%elem = transpose(conjg(ope%elem))
    end function


    subroutine isher(A,her)
        !�ж�����Ƿ����
        !���������
            !A   :   һ�����ʵ��
            !her :   ���ص��߼�ֵ
        !���������
            !her :   ����Ϊ�棬�����
        type(opera), intent(in) :: A
        logical, intent(out) :: her
        complex(8), allocatable :: mid(:,:)

        integer :: ii,jj

        allocate(mid(A%dims,A%dims))
        mid = A%elem - transpose(conjg(A%elem))

        her = .true.

        do ii = 1,A%dims
            if (abs(aimag(A%elem(ii,ii))) > 1e-15) then
                her = .false.
                return
            end if
            do jj = 1,A%dims
                if (abs(mid(jj,ii)) > 1e-15) then
                    her = .false.
                    return
                end if
            end do
        end do

    end subroutine


    function density(psi)
        !����psi���ܶȾ�����ʽ
        !���������
            !psi   :   һ��ʸ��ʵ��
        !�������:
            !density  :   �����ܶȾ���pho = |psi><psi|
        type(vecs), intent(in) :: psi
        type(opera) :: density

        integer :: i

        if (psi%lengs == -1) stop 'Wrong density function: State not defined'

        density%dims = psi%lengs
        allocate(density%subs(size(psi%subs)))

        density%subs = psi%subs


        allocate(density%elem(psi%lengs,psi%lengs))
        if (psi%k_b .eq. 'ket') then
            do i = 1,psi%lengs
                density%elem(i,:) = conjg(psi%elem) * psi%elem(i)
            end do
        else
            do i = 1,psi%lengs
                density%elem(i,:) = psi%elem * conjg(psi%elem(i))
            end do
        end if
    end function

    function trace(A)
        type(opera), intent(in) :: A
        complex*8 :: trace

        integer :: i

        if (A%dims == -1) stop 'Wrong call trace: Operator not defined'

        trace = dcmplx(0)
        do i = 1,A%dims
            trace = trace + A%elem(i,i)
        end do
    end function

    subroutine eig_hermitian(A,eigval,eigvec)
        type(opera), intent(in) :: A
        real(kind=8), allocatable, intent(out) :: eigval(:)
        type(vecs), allocatable, intent(out), optional :: eigvec(:)

        complex(kind=8), allocatable :: eigmid(:,:)
        integer :: i

        if (A%dims == -1) stop 'Wrong call eig: Operator not defined'


        if (.not. allocated(eigval)) allocate(eigval(A%dims))
        if (.not. present(eigvec)) then
            call eig_o(A%elem,eigval)
        else
            if (.not. allocated(eigvec)) allocate(eigvec(A%dims))
            allocate(eigmid(A%dims,A%dims))
            eigvec(:)%lengs = A%dims
!            eigvec(:)%k_b = 'ket'
            call eig_o(A%elem,eigval,eigmid)
            do i = 1,A%dims
                allocate(eigvec(i)%subs(size(A%subs)))
                eigvec(i)%subs = A%subs
                allocate(eigvec(i)%elem(A%dims))
                eigvec(i)%elem(:) = eigmid(:,i)
            end do
        end if

        deallocate(eigmid)

    end subroutine

    subroutine eig_non_hermitian(A,eigval,eigvec)
        type(opera), intent(in) :: A
        complex(kind=8), allocatable, intent(out) :: eigval(:)
        type(vecs), allocatable, intent(out), optional :: eigvec(:)

        complex(kind=8), allocatable :: eigmid(:,:)
        integer :: i

        if (A%dims == -1) stop 'Wrong call eig: Operator not defined'

        if (.not. allocated(eigval)) allocate(eigval(A%dims))
        if (.not. present(eigvec)) then
            call eig_gene(A%elem,eigval)
        else
            if (.not. allocated(eigvec)) allocate(eigvec(A%dims))
            allocate(eigmid(A%dims,A%dims))
            eigvec(:)%lengs = A%dims
            call eig_gene(A%elem,eigval,eigmid)
            do i = 1,A%dims
                allocate(eigvec(i)%subs(size(A%subs)))
                allocate(eigvec(i)%elem(A%dims))
                eigvec(i)%elem(:) = eigmid(:,i)
            end do
        end if

        deallocate(eigmid)
    end subroutine
end module
