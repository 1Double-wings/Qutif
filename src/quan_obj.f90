module quan_obj
    !定义各种常用算符以及矢量
    implicit none

    !结构类型：算符
    !成员:
        !dims   :   算符(方阵)的维度
        !elem   :   算符的矩阵元
        !subs   :   算符所属子(复合)空间
    type opera
        integer :: dims = -1
        complex(8), allocatable :: elem(:,:)
        integer, allocatable :: subs(:)
    end type
    public :: opera


    !结构类型：矢量
    !成员:
        !lengs   :   矢量(方阵)的长度
        !k_b     :   矢量类型(ket/bra)
        !elem    :   矢量的元素
        !subs    :   矢量所属子(复合)空间
    type vecs
        integer :: lengs = -1
        integer, allocatable :: subs(:)
        character(len=3) :: k_b = 'ket'
        complex(8), allocatable :: elem(:)
    end type
    public :: vecs


    !常用的算符：
        !destory  :  玻色子湮灭算符
        !basis    :  生成ket
        !idt      :  生成单位矩阵
        !sigmax,y,z  :  Pauli矩阵
        !sigmam,p :  升降(二分量)算符
    contains

    function destroy(N)
        !输入参数：
            !N   :   算符的维度
        !输出参数:
            !destroy   :   返回destroy为湮灭算符
        type(opera) :: destroy
        integer, intent(in) :: N

        integer :: ii
        complex(8) :: jj

        destroy%dims = N
        allocate(destroy%subs(1))
        destroy%subs = N

        allocate(destroy%elem(N,N))
        do ii = 2,N
            jj = cmplx(ii-1,0)
            destroy%elem(ii-1,ii) = sqrt(jj)
        end do
    end function destroy

    function basis(N,num)
        !输入参数：
            !N     :   矢量的长度
            !num   :   第num个基矢
        !输出参数：
            !basis   :   返回basis为一个ket
        type(vecs) :: basis
        integer, intent(in) :: N, num

        if (N<=0) stop 'N is an integer >0'
        if (num>N) stop 'num is an integer<=N'

        basis%lengs = N
        basis%k_b = 'ket'

        allocate(basis%subs(1))
        basis%subs = N

        allocate(basis%elem(N))
        basis%elem(:) = cmplx(0,0)
        if (num .ne. 0) basis%elem(num) = cmplx(1,0)
    end function basis

    function idt(N)
        !输入参数：
            !N       :   算符维度
        !输出参数：
            !idt   :   返回一个恒等算符
        type(opera) :: idt
        integer, intent(in) :: N(:)

        integer :: i

        idt%dims = 1
        allocate(idt%subs(size(N)))
        idt%subs = N

        do i = 1,size(N)
            idt%dims = idt%dims * N(i)
        end do

        allocate(idt%elem(idt%dims,idt%dims))
        idt%elem = 0
        do i = 1,idt%dims
            idt%elem(i,i) = 1
        end do
    end function

    function sigmax()
        !输出参数：
            !sigmax   :   返回Pauli算符Sx
        type(opera) :: sigmax

        sigmax%dims = 2
        allocate(sigmax%subs(1))
        sigmax%subs = 2

        allocate(sigmax%elem(2,2))

        sigmax%elem = 0
        sigmax%elem(1,2) = 1
        sigmax%elem(2,1) = 1
        end function sigmax

    function sigmay()
        type(opera) :: sigmay

        sigmay%dims = 2

        allocate(sigmay%subs(1))
        sigmay%subs = 2


        allocate(sigmay%elem(2,2))

        sigmay%elem = 0
        sigmay%elem(1,2) = cmplx(0,-1)
        sigmay%elem(2,1) = cmplx(0,1)
        end function sigmay

    function sigmaz()
        type(opera) :: sigmaz

        sigmaz%dims = 2

        allocate(sigmaz%subs(1))
        sigmaz%subs = 2

        allocate(sigmaz%elem(2,2))

        sigmaz%elem = 0
        sigmaz%elem(1,1) = 1
        sigmaz%elem(2,2) = -1
        end function sigmaz

    function sigmap()
        type(opera) :: sigmap

        sigmap%dims = 2

        allocate(sigmap%subs(1))
        sigmap%subs = 2

        allocate(sigmap%elem(2,2))

        sigmap%elem = 0
        sigmap%elem(1,2) = 1
        end function sigmap

    function sigmam()
        type(opera) :: sigmam

        sigmam%dims = 2

        allocate(sigmam%subs(1))
        sigmam%subs = 2

        allocate(sigmam%elem(2,2))

        sigmam%elem = 0
        sigmam%elem(2,1) = 1
        end function sigmam
end module
