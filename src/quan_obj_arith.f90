module quan_obj_arith
    use quan_obj
    use math_funs
    use quan_obj_funs

    implicit none

    interface Operator(*)
        module procedure :: nvtime_c, nvtime_r, vntime_c, vntime_r, notime_c, notime_r, ontime_c, &
        ontime_r, vvtime, votime, ovtime, ootime
    end interface

    interface Operator(+)
        module procedure :: ooplus, vvplus
    end interface

    interface Operator(-)
        module procedure :: oominus, vvminus
    end interface

    interface Operator(/)
        module procedure :: vnover_c, vnover_r, onover_c, onover_r
    end interface

    interface assignment(=)
        module procedure :: vveq, mmeq
    end interface

    interface Operator(.ot.)
        module procedure :: tensoroo, tensorvv
    end interface


    contains


    !包含各种结构之间的乘法：
        !nvtime   :   标量矢量相乘
        !vntime/over   :   矢量标量相乘/相除
        !notime/over   :   标量矩阵相乘/相除，即n*inv(o)
        !ontime/over   :   矩阵标量相乘/相除
        !vvtime   :   矢量矢量相乘，结果为矩阵
        !votime   :   矢量矩阵相乘，结果为bra
        !ovtime   :   矩阵矢量相乘，结果为ket
        !ootime   :   矩阵矩阵相乘，结果为矩阵
        !ooplus/minus   :   矩阵矩阵相加/相减
        !vvplus/minus   :   矢量矢量相加/相减
        !vveq   :   矢量矢量相等
        !ooeq   :   矩阵矩阵相等
        !tensoroo   :   矩阵矩阵直积
        !tensorvv   :   矢量矢量直积

    function nvtime_c(A,B)
        complex(8), intent(in) :: A
        type(vecs), intent(in) :: B
        type(vecs) :: nvtime_c

        if (B%lengs == -1) stop 'state is not defined'

        allocate(nvtime_c%subs(size(B%subs)))
        nvtime_c%subs = B%subs

        nvtime_c%k_b = B%k_b
        nvtime_c%lengs = B%lengs

        allocate(nvtime_c%elem(B%lengs))
        nvtime_c%elem = A*B%elem
    end function

    function nvtime_r(A,B)
        real(8), intent(in) :: A
        type(vecs), intent(in) :: B
        type(vecs) :: nvtime_r

        if (B%lengs == -1) stop 'state is not defined'

        allocate(nvtime_r%subs(size(B%subs)))
        nvtime_r%subs = B%subs

        nvtime_r%k_b = B%k_b
        nvtime_r%lengs = B%lengs

        allocate(nvtime_r%elem(B%lengs))
        nvtime_r%elem = A*B%elem
    end function

    function vntime_c(A,B)
        complex(8), intent(in) :: B
        type(vecs), intent(in) :: A
        type(vecs) :: vntime_c

        if (A%lengs == -1) stop 'state is not defined'

        allocate(vntime_c%subs(size(A%subs)))
        vntime_c%subs = A%subs

        vntime_c%k_b = A%k_b
        vntime_c%lengs = A%lengs

        allocate(vntime_c%elem(A%lengs))
        vntime_c%elem = B*A%elem
    end function

    function vntime_r(A,B)
        real(8), intent(in) :: B
        type(vecs), intent(in) :: A
        type(vecs) :: vntime_r

        if (A%lengs == -1) stop 'state is not defined'

        allocate(vntime_r%subs(size(A%subs)))
        vntime_r%subs = A%subs

        vntime_r%k_b = A%k_b
        vntime_r%lengs = A%lengs

        allocate(vntime_r%elem(A%lengs))
        vntime_r%elem = B*A%elem
    end function

    function notime_c(A,B)
        complex(8), intent(in) :: A
        type(opera), intent(in) :: B
        type(opera) :: notime_c

        if (B%dims == -1) stop 'operator is not defined'

        allocate(notime_c%subs(size(B%subs)))
        notime_c%subs = B%subs

        notime_c%dims = B%dims

        allocate(notime_c%elem(B%dims,B%dims))
        notime_c%elem = A*B%elem
    end function

    function notime_r(A,B)
        real(8), intent(in) :: A
        type(opera), intent(in) :: B
        type(opera) :: notime_r

        if (B%dims == -1) stop 'operator is not defined'

        allocate(notime_r%subs(size(B%subs)))
        notime_r%subs = B%subs

        notime_r%dims = B%dims

        allocate(notime_r%elem(B%dims,B%dims))
        notime_r%elem = B%elem(:,:)*A
    end function

    function ontime_c(A,B)
        complex(8), intent(in) :: B
        type(opera), intent(in) :: A
        type(opera) :: ontime_c

        if (A%dims == -1) stop 'operator is not defined'

        allocate(ontime_c%subs(size(A%subs)))
        ontime_c%subs = A%subs

        ontime_c%dims = A%dims

        allocate(ontime_c%elem(A%dims,A%dims))
        ontime_c%elem = B*A%elem
    end function

    function ontime_r(A,B)
        real(8), intent(in) :: B
        type(opera), intent(in) :: A
        type(opera) :: ontime_r

        if (A%dims == -1) stop 'operator is not defined'

        allocate(ontime_r%subs(size(A%subs)))
        ontime_r%subs = A%subs

        ontime_r%dims = A%dims

        allocate(ontime_r%elem(A%dims,A%dims))
        ontime_r%elem = B*A%elem
    end function

    function vvtime(A,B)
        type(vecs), intent(in) :: A
        type(vecs), intent(in) :: B
        type(opera) :: vvtime

        integer :: i


        if (A%lengs == -1 .or. B%lengs == -1) stop 'Wrong vector product: State not defined'
        if (A%k_b .eq. 'bra') stop 'Wrong vector product: ket*bra is required'
        if (B%k_b .eq. 'ket') stop 'Wrong vector product: ket*bra is required'
        if (A%lengs /= B%lengs) stop 'Wrong vector product: Mismatched vector length of states'

        if (size(A%subs) /= size(B%subs)) stop 'Wrong vector product: Mismatched subsystems of states'

        do i = 1,size(A%subs)
            if (A%subs(i) /= B%subs(i)) stop 'Wrong vector product: Mismatched subsystems of states'
        end do

        if (A%k_b .ne. B%k_b) then
            vvtime%dims = A%lengs
            allocate(vvtime%subs(size(A%subs)))
            vvtime%subs = A%subs

            allocate(vvtime%elem(A%lengs,A%lengs))

            do i = 1,A%lengs
                vvtime%elem(:,i) = A%elem * B%elem(i)
            end do
        end if

    end function

    function votime(A,B)
        type(vecs), intent(in) :: A
        type(opera), intent(in) :: B
        type(vecs) :: votime, Adag

        integer :: i,incx=1,incy=1
        complex(8) :: alpha=cmplx(1.0,0.0),beta=cmplx(0.0,0.0)

        if (A%lengs == -1) stop 'Wrong vec*matrix: State not defined'
        if (B%dims == -1) stop 'Wrong vec*matrix: Operator not defined'
        if (A%k_b .eq. 'ket') stop 'Wrong vec*matrix: Type bra is required'
        if (size(A%subs) /= size(B%subs)) stop 'Wrong vec*matrix: Mismatched subsystems between state and operator'
        if (A%lengs /= B%dims) stop 'Wrong state*operator: Mismatched dimensions between state and operator'

        allocate(votime%subs(size(A%subs)))
        do i = 1,size(A%subs)
            votime%subs(i) = A%subs(i)
            if (A%subs(i) /= B%subs(i)) stop 'Wrong vec*matrix: Mismatched subsystems between state and operator'
        end do

        votime%k_b = A%k_b
        votime%lengs = A%lengs

        allocate(votime%elem(A%lengs))

        Adag = dag(A)
        call zgemv('c',B%dims,B%dims,alpha,B%elem,B%dims,Adag%elem,incx,beta,votime%elem,incy)
        votime%elem = conjg(votime%elem)
    end function

    function ovtime(A,B)
        type(vecs), intent(in) :: B
        type(opera), intent(in) :: A
        type(vecs) :: ovtime

        integer :: i,incx=1,incy=1
        complex(8) :: alpha=cmplx(1.0,0.0),beta=cmplx(0.0,0.0)

        if (B%lengs == -1) stop 'Wrong matrix*vec: State not defined'
        if (A%dims == -1) stop 'Wrong matrix*vec: Operator not defined'
        if (B%k_b .eq. 'bra') stop 'Wrong matrix*vec: State type bra is required'
        if (size(A%subs) /= size(B%subs)) stop 'Wrong matrix*vec: Mismatched subsystems between state and operator'
        if (A%dims /= B%lengs) stop 'Wrong operator*state: Mismatched dimensions between state and operator'

        allocate(ovtime%subs(size(B%subs)))
        do i = 1,size(A%subs)
            ovtime%subs(i) = A%subs(i)
            if (A%subs(i) /= B%subs(i)) stop 'Wrong matrix*vec: Mismatched subsystems between state and operator'
        end do

        ovtime%k_b = B%k_b
        ovtime%lengs = B%lengs

        allocate(ovtime%elem(B%lengs))

        call zgemv('n',A%dims,A%dims,alpha,A%elem,A%dims,B%elem,incx,beta,ovtime%elem,incy)
    end function

    function ootime(A,B)
        type(opera), intent(in) :: A,B
        type(opera) :: ootime

        integer :: i
        complex(8) :: alpha=1d0,beta=0d0

        if (A%dims == -1) stop 'Wrong operator*operator: Operator not defined'
        if (B%dims == -1) stop 'Wrong operator*operator: Operator not defined'
        if (size(A%subs) /= size(B%subs)) stop 'Wrong operator*operator: Mismatched subsystems of operators'
        if (A%dims /= B%dims) stop 'Wrong operator*operator: Mismatched dimensions of operators'

        allocate(ootime%subs(size(A%subs)))
        do i = 1,size(A%subs)
            ootime%subs(i) = A%subs(i)
            if (A%subs(i) /= B%subs(i)) stop 'Wrong operator*operator: Mismatched subsystems of operators'
        end do

        ootime%dims = A%dims

        allocate(ootime%elem(A%dims,A%dims))

        call zgemm('n','n',A%dims,A%dims,A%dims,alpha,A%elem,A%dims,B%elem,B%dims,beta,ootime%elem,ootime%dims)
    end function



    function ooplus(A,B)
        type(opera), intent(in) :: A,B
        type(opera) :: ooplus

        integer :: i

        if (A%dims == -1 .or. B%dims == -1) stop 'Wrong operator+operator: Operator not defined'
        if (A%dims /= B%dims) stop 'Wrong operator+operator: Mismatched operators dimension'
        if (size(A%subs) /= size(B%subs)) stop 'Wrong operator+operator: Mismatched subsystems of operator'

        do i = 1,size(A%subs)
            if (A%subs(i) /= B%subs(i)) stop 'Wrong operator+operator: Mismatched subsystems of operator'
        end do

        ooplus%dims = A%dims

        allocate(ooplus%subs(size(A%subs)))
        ooplus%subs = A%subs

        allocate(ooplus%elem(size(A%elem,1),size(A%elem,2)))
        ooplus%elem = A%elem + B%elem
    end function

    function vvplus(A,B)
        type(vecs), intent(in) :: A,B
        type(vecs) :: vvplus

        integer :: i

        if (A%lengs == -1 .or. B%lengs == -1) stop 'Wrong state+state: State not defined'
        if (A%k_b .ne. B%k_b) stop 'Wrong state+state: Mismatched type of state'
        if (A%lengs /= B%lengs) stop 'Wrong state+state: Mismatched states length'
        if (size(A%subs) /= size(B%subs)) stop 'Wrong state+state: Mismatched subsystems of states'

        do i = 1,size(A%subs)
            if (A%subs(i) /= B%subs(i)) stop 'Wrong state+state: Mismatched subsystems of states'
        end do

        vvplus%lengs = A%lengs
        vvplus%k_b = A%k_b

        allocate(vvplus%subs(size(A%subs)))
        vvplus%subs = A%subs

        allocate(vvplus%elem(A%lengs))
        vvplus%elem = A%elem + B%elem
    end function




    function oominus(A,B)
        type(opera), intent(in) :: A,B
        type(opera) :: oominus

        integer :: i

        if (A%dims == -1 .or. B%dims == -1) stop 'Wrong operator+operator: Operator not defined'
        if (A%dims /= B%dims) stop 'Wrong operator+operator: Mismatched operators dimension'
        if (size(A%subs) /= size(B%subs)) stop 'Wrong operator+operator: Mismatched subsystems of operator'

        do i = 1,size(A%subs)
            if (A%subs(i) /= B%subs(i)) stop 'Wrong operator+operator: Mismatched subsystems of operator'
        end do

        oominus%dims = A%dims

        allocate(oominus%subs(size(A%subs)))
        oominus%subs = A%subs

        allocate(oominus%elem(size(A%elem,1),size(A%elem,2)))
        oominus%elem = A%elem - B%elem
    end function

    function vvminus(A,B)
        type(vecs), intent(in) :: A,B
        type(vecs) :: vvminus

        integer :: i

        if (A%lengs == -1 .or. B%lengs == -1) stop 'Wrong state+state: State not defined'
        if (A%k_b .ne. B%k_b) stop 'Wrong state+state: Mismatched type of state'
        if (A%lengs /= B%lengs) stop 'Wrong state+state: Mismatched states length'
        if (size(A%subs) /= size(B%subs)) stop 'Wrong state+state: Mismatched subsystems of states'

        do i = 1,size(A%subs)
            if (A%subs(i) /= B%subs(i)) stop 'Wrong state+state: Mismatched subsystems of states'
        end do

        vvminus%lengs = A%lengs
        vvminus%k_b = A%k_b

        allocate(vvminus%subs(size(A%subs)))
        vvminus%subs = A%subs

        allocate(vvminus%elem(A%lengs))
        vvminus%elem = A%elem - B%elem
    end function




    function vnover_c(A,B)
        type(vecs), intent(in) :: A
        complex(8), intent(in) :: B
        type(vecs) :: vnover_c

        if (A%lengs == -1) stop 'Wrong state/number: State not defined'

        vnover_c%lengs = A%lengs
        vnover_c%k_b = A%k_b

        allocate(vnover_c%subs(size(A%subs)))
        vnover_c%subs = A%subs

        allocate(vnover_c%elem(A%lengs))
        vnover_c%elem = 1/B*A%elem
    end function

    function vnover_r(A,B)
        type(vecs), intent(in) :: A
        real(8), intent(in) :: B
        type(vecs) :: vnover_r

        if (A%lengs == -1) stop 'Wrong state/number: State not defined'

        vnover_r%lengs = A%lengs
        vnover_r%k_b = A%k_b

        allocate(vnover_r%subs(size(A%subs)))
        vnover_r%subs = A%subs

        allocate(vnover_r%elem(A%lengs))
        vnover_r%elem = 1/B*A%elem
    end function

    function onover_c(A,B)
        type(opera), intent(in) :: A
        complex(8), intent(in) :: B
        type(opera) :: onover_c

        if (A%dims == -1) stop 'Wrong operator/number: Operator not defined'

        onover_c%dims = A%dims

        allocate(onover_c%subs(size(A%subs)))
        onover_c%subs = A%subs

        allocate(onover_c%elem(A%dims,A%dims))
        onover_c%elem = 1/B*A%elem
    end function

    function onover_r(A,B)
        type(opera), intent(in) :: A
        real, intent(in) :: B
        type(opera) :: onover_r

        if (A%dims == -1) stop 'Wrong operator/number: Operator not defined'

        onover_r%dims = A%dims

        allocate(onover_r%subs(size(A%subs)))
        onover_r%subs = A%subs

        allocate(onover_r%elem(A%dims,A%dims))
        onover_r%elem = 1/B*A%elem
    end function



    subroutine vveq(res,vin)
        type(vecs), intent(in) :: vin
        type(vecs), intent(out) :: res

        res%lengs = vin%lengs
        res%k_b = vin%k_b

        allocate(res%subs(size(vin%subs)))
        res%subs = vin%subs

        allocate(res%elem(vin%lengs))
        res%elem = vin%elem

    end subroutine

    subroutine mmeq(res,oin)
        type(opera), intent(in) :: oin
        type(opera), intent(out) :: res

        res%dims = oin%dims

        allocate(res%subs(size(oin%subs)))
        res%subs = oin%subs

        allocate(res%elem(oin%dims,oin%dims))
        res%elem = oin%elem

    end subroutine



    function tensoroo(A,B)
        type(opera), intent(in) :: A,B
        type(opera) :: tensoroo


        if (A%dims == -1 .or. B%dims == -1) stop 'Wrong operator tensor product: Operator not defined'


        allocate(tensoroo%subs(size(A%subs) + size(B%subs)))

        tensoroo%dims = A%dims * B%dims
        tensoroo%subs(1 : size(A%subs)) = A%subs
        tensoroo%subs(size(A%subs)+1 : size(B%subs)) = B%subs
        call tensor(A%elem,B%elem,tensoroo%elem)

    end function

    function tensorvv(A,B)
        type(vecs), intent(in) :: A,B
        type(vecs) :: tensorvv


        if (A%lengs == -1 .or. B%lengs == -1) stop 'Wrong state tensor product: State not defined'
        if (A%k_b .ne. B%k_b) stop 'Wrong state tensor product: Mismatched vector type of states'


        allocate(tensorvv%subs(size(A%subs) + size(B%subs)))
        tensorvv%k_b = A%k_b

        tensorvv%lengs = A%lengs * B%lengs
        tensorvv%subs(1 : size(A%subs)) = A%subs
        tensorvv%subs(size(A%subs)+1 : size(B%subs)) = B%subs
        call tensor(A%elem,B%elem,tensorvv%elem)
    end function

end module
