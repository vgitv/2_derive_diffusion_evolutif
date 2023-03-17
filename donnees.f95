! ===========================================================================================================
! ===========================================================================================================

MODULE donnees

    use maillage

    implicit none

    real(rp), dimension(:, :), allocatable :: A
    !real(rp), parameter :: delta_n = 1.0_rp
    !real(rp), parameter :: delta_p = 1.0_rp
    type(Mesh), save :: maill

    ! cste telle que nr * pr = nl * pl = C0 pour CB compatibles avec éq. thermique
    real(rp), parameter :: C0 = 1.0_rp

contains

    ! -------------------------------------------------------------------------------------------------------
    ! dopage
    ! -------------------------------------------------------------------------------------------------------
    function c_dopage(x)
        ! paramètres
        real(rp), intent(in) :: x

        ! return
        real(rp) :: c_dopage

        c_dopage = 10.0_rp
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Conditions initiales
    ! -------------------------------------------------------------------------------------------------------
    function n0(x)
        real(rp), intent(in) :: x
        real(rp) :: n0
        n0 = 0.5_rp
    end function

    function p0(x)
        real(rp), intent(in) :: x
        real(rp) :: p0
        p0 = 1.0_rp
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Conditions au bord
    ! -------------------------------------------------------------------------------------------------------
    function n_l(t)
        real(rp), intent(in) :: t
        real(rp) :: n_l
        n_l = 1.0_rp
    end function

    function n_r(t)
        real(rp), intent(in) :: t
        real(rp) :: n_r
        n_r = 1.0_rp
    end function

    function p_l(t)
        real(rp), intent(in) :: t
        real(rp) :: p_l
        p_l = C0 / n_l(t)
    end function

    function p_r(t)
        real(rp), intent(in) :: t
        real(rp) :: p_r
        p_r = C0 / n_r(t)
    end function

    function psi_l(t)
        real(rp), intent(in) :: t
        real(rp) :: psi_l
        psi_l = 0.5_rp * (log(n_l(t)) - log(p_l(t)))
    end function

    function psi_r(t)
        real(rp), intent(in) :: t
        real(rp) :: psi_r
        psi_r = 0.5_rp * (log(n_r(t)) - log(p_r(t)))
    end function

END MODULE donnees
