module quan_obj
    !������ֳ�������Լ�ʸ��
    implicit none

    !�ṹ���ͣ����
    !��Ա:
        !dims   :   ���(����)��ά��
        !elem   :   ����ľ���Ԫ
        !subs   :   ���������(����)�ռ�
    type opera
        integer :: dims = -1
        complex(8), allocatable :: elem(:,:)
        integer, allocatable :: subs(:)
    end type
    public :: opera


    !�ṹ���ͣ�ʸ��
    !��Ա:
        !lengs   :   ʸ��(����)�ĳ���
        !k_b     :   ʸ������(ket/bra)
        !elem    :   ʸ����Ԫ��
        !subs    :   ʸ��������(����)�ռ�
    type vecs
        integer :: lengs = -1
        integer, allocatable :: subs(:)
        character(len=3) :: k_b = 'ket'
        complex(8), allocatable :: elem(:)
    end type
    public :: vecs


    !���õ������
        !destory  :  ��ɫ���������
        !basis    :  ����ket
        !idt      :  ���ɵ�λ����
        !sigmax,y,z  :  Pauli����
        !sigmam,p :  ����(������)���
    contains

    function destroy(N)
        !���������
            !N   :   �����ά��
        !�������:
            !destroy   :   ����destroyΪ�������
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
        !���������
            !N     :   ʸ���ĳ���
            !num   :   ��num����ʸ
        !���������
            !basis   :   ����basisΪһ��ket
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
        !���������
            !N       :   ���ά��
        !���������
            !idt   :   ����һ��������
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
        !���������
            !sigmax   :   ����Pauli���Sx
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
