! ===========================================================================================================
! Variable structurée Mesh et fonction relatives
!
! Maillage avec interfaces pour schémas de type volumes finis. Sous programmes :
!
! newMesh
! rmMesh
! ===========================================================================================================

MODULE maillage

    use math

    implicit none

    type :: Mesh
        ! vecteur de longueur l + 2 : x_{1} ... x_{l+2}
        integer :: l

        ! points du maillage x_{i}
        real(rp), dimension(:), allocatable :: x
        ! points milieu des mailles x_{i+1/2}
        real(rp), dimension(:), allocatable :: x2

        ! distance h_{i} autour de x_{i}, ie entre x_{i-1/2} et x_{i+1/2}
        real(rp), dimension(:), allocatable :: h
        ! distance h_{i+1/2} autour de x_{i+1/2}, ie entre x_{i} et x_{i+1}
        real(rp), dimension(:), allocatable :: h2
    end type Mesh

contains

    ! -------------------------------------------------------------------------------------------------------
    ! constructeur Mesh
    ! -------------------------------------------------------------------------------------------------------
    subroutine newMesh(x, m)
        ! paramètres
        real(rp), dimension(:) :: x
        type(Mesh) :: m

        ! variables locales
        integer :: l, i

        l = size(x) - 2
        allocate(m%x(l + 2), m%x2(l + 1), m%h(l), m%h2(l + 1))

        ! affectation de l
        m%l = l

        ! affectation de x
        m%x = x

        ! affectation de h2 et x2
        do i = 1, l + 1
            m%h2(i) = x(i + 1) - x(i)
            m%x2(i) = (x(i) + x(i + 1)) / 2.0_rp
        end do
        m%x2(1) = x(1)
        m%x2(l + 1) = x(l + 2)

        ! affectation de h
        do i = 1, l
            m%h(i) = m%x2(i + 1) - m%x2(i)
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! destructeur Mesh
    ! -------------------------------------------------------------------------------------------------------
    subroutine rmMesh(m)
        ! paramètres
        type(Mesh) :: m

        if (allocated(m%x)) then
            deallocate(m%x, m%h, m%x2, m%h2)
        end if
    end subroutine

END MODULE maillage
