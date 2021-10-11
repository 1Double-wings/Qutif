module time_evl
    use quan_obj
    use quan_obj_funs
    use quan_obj_arith
    implicit none

    contains

    subroutine kpsolve(tlist, H, psi0, colla, expect_o, expect_val, phot)
        real(kind=8), intent(in) :: tlist(:)
        type(opera), intent(in) :: H
        type(vecs), intent(in) :: psi0
        type(opera), intent(in) :: expect_o(:), colla(:)
        type(opera), intent(out), optional :: phot
        real(kind=8), allocatable, intent(out) :: expect_val(:,:)


        type(opera) :: m0, pho0
        logical :: hermi
        integer :: i,j,k
        real(kind=8) :: dt


        if (tlist(1)<0 .or. tlist(size(tlist))<0) stop 'Wrong Kpsolve: Time list must > 0'
        if (H%dims == -1) stop 'Wrong Kpsolve: H not defined'
        call isher(H,hermi)
        if (.not. hermi) stop 'Wrong Kpsolve: H is non-Hermitian'
        if (psi0%lengs == -1) stop 'Wrong Kpsolve: Initial state not defined'
        if (psi0%k_b .ne. 'ket') stop 'Wrong Kpsolve: Initial state is not a ket'
        if (size(H%subs) /= size(psi0%subs)) stop 'Wrong Kpsolve: Mismatched subsystems of operator and state'

        do i = 1,size(H%subs)
            if (H%subs(i) /= psi0%subs(i)) stop 'Wrong Kpsolve: Mismatched subsystems of operator and state'
        end do

        if (.not. allocated(expect_val)) then
            allocate(expect_val(size(expect_o), size(tlist)))
        elseif (size(expect_val,1) /= size(expect_o)) then
            deallocate(expect_val)
            allocate(expect_val(size(expect_o), size(tlist)))
        end if

        pho0 = density(psi0)
        dt = tlist(2) - tlist(1)
        m0 = idt(H%subs) - dcmplx(0d0,dt) * H

        do i = 1,size(colla)
            m0 = m0 - 0.5d0 * dt * dag(colla(i)) * colla(i)
        end do


        do i = 1,size(tlist)
            pho0 = m0*pho0*dag(m0)
            do j = 1,size(colla)
                pho0 = pho0 + dt * colla(j) * pho0 * dag(colla(j))
            end do
            do k = 1,size(expect_o)
                expect_val(k, i) = real(trace(pho0*expect_o(k)))
            end do
        end do

        if (present(phot)) phot = pho0
    end subroutine
end module
