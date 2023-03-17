MODULE donnees

    use math
    use maillage
    use variables

    implicit none

contains

    ! =======================================================================================================
    ! PROBLÈME DE POISSON
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! second membre du problème de Poisson
    ! -------------------------------------------------------------------------------------------------------
    function f(x)
        ! paramètres
        real(rp), intent(in) :: x

        ! return
        real(rp) :: f

        f = 1.0_rp
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Solution de l'équation de Poisson avec f = 1 et ul = ur = 0
    ! -------------------------------------------------------------------------------------------------------
    function u(x)
        ! paramètres
        real(rp), intent(in) :: x

        ! return
        real(rp) :: u

        u = x * (1.0_rp - x) / 2.0_rp
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Construction matrice A pour résolution A x = b dans le cas de l'équation de Poisson (idem éq thermique)
    ! -------------------------------------------------------------------------------------------------------
    subroutine build_A(m, A)
        ! paramètres
        type(Mesh), intent(in) :: m
        real(rp), dimension(:, :), intent(out) :: A

        ! variables locales
        integer :: i

        A = 0.0_rp
        do i = 1, m%l - 1
            A(i, i) = 1.0_rp / m%h2(i) + 1.0_rp / m%h2(i + 1)
            A(i, i + 1) = -1.0_rp / m%h2(i + 1)
            A(i + 1, i) = A(i, i + 1)
        end do
        A(m%l, m%l) = 1.0_rp / m%h2(m%l) + 1.0_rp / m%h2(m%l + 1)
    end subroutine build_A



    ! -------------------------------------------------------------------------------------------------------
    ! Construction vecteur b pour résolution A x = b dans le cas de l'équation de Poisson
    ! -------------------------------------------------------------------------------------------------------
    subroutine build_b(m, fc, ul, ur, b)
        ! paramètres
        type(Mesh), intent(in) :: m
        real(rp), dimension(:), intent(out) :: b
        real(rp), external :: fc
        real(rp), intent(in) :: ul, ur

        ! variables locales
        integer :: i

        do i = 1, m%l
            b(i) = m%h(i) * f(m%x(i + 1))
        end do

        b(1) = b(1) + ul / m%h2(1)
        b(m%l) = b(m%l) + ur / m%h2(m%l + 1)
    end subroutine build_b



    ! =======================================================================================================
    ! NEWTON 1D
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! fonction test pour algo newton 1D : une fonction g et sa dérivée gp ( = gprime)
    ! -------------------------------------------------------------------------------------------------------
    function g(x)
        ! paramètres
        real(rp), intent(in) :: x

        ! return
        real(rp) :: g

        g = (x - 1.0_rp) * (x + 1.0_rp) * (x - 3.0_rp)
    end function

    function gp(x)
        ! paramètres
        real(rp), intent(in) :: x

        ! return
        real(rp) :: gp

        gp = (x + 1.0_rp) * (x - 3.0_rp) + (x - 1.0_rp) * (x - 3.0_rp) + (x - 1.0_rp) * (x + 1.0_rp)
    end function



    ! =======================================================================================================
    ! Schéma explicite
    ! =======================================================================================================
    ! -------------------------------------------------------------------------------------------------------
    ! Calcul second membre A_psi psi = b_psi. La sol donne psi pour n et p donnés (temps fixé)
    ! -------------------------------------------------------------------------------------------------------
    ! m : maillage d'espace
    ! n : approx de n à un temps donné sur les elts du maillage d'espace
    ! p : approx de p à un temps donné sur les elts du maillage d'espace
    ! c : dopage sur les pts du maillage
    ! psi_gauche : valeur à gauche (condition bord) pour un temps donné
    ! psi_droit  : valeur à droite (condition bord) pour un temps donné
    ! b_psi : vecteur retourné
    subroutine build_b_psi(m, n, p, c, psi_gauche, psi_droit, b_psi)
        ! paramètres
        type(Mesh), intent(in) :: m
        real(rp), dimension(:), intent(in) :: n, p, c
        real(rp), intent(in) :: psi_gauche, psi_droit
        real(rp), dimension(:), intent(out) :: b_psi

        ! variables locales
        integer :: i

        do i = 1, m%l
            b_psi(i) = -m%h(i) * (n(i) - p(i) - c(i))
        end do

        b_psi(1) = b_psi(1) + psi_gauche / m%h2(1)
        b_psi(m%l) = b_psi(m%l) + psi_droit / m%h2(m%l + 1)
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Pour n et p données (temps fixé), calcul psi au même temps
    ! -------------------------------------------------------------------------------------------------------
    ! m : maillage d'espace
    ! n : approx de n à un temps donné sur les elts du maillage d'espace
    ! p : approx de p à un temps donné sur les elts du maillage d'espace
    ! c : dopage sur les pts du maillage
    ! psi_gauche : valeur à gauche (condition bord) pour un temps donné
    ! psi_droit  : valeur à droite (condition bord) pour un temps donné
    ! psi : vecteur retourné
    subroutine potentiel(m, n, p, c, psi_gauche, psi_droit, psi)
        ! paramètres
        type(Mesh), intent(in) :: m
        real(rp), dimension(:), intent(in) :: n, p, c
        real(rp), intent(in) :: psi_gauche, psi_droit
        real(rp), dimension(:), intent(out) :: psi

        ! variables locales
        real(rp), dimension(m%l) :: b_psi
        real(rp), dimension(m%l, m%l) :: A_psi

        call build_A(m, A_psi)
        call build_b_psi(m, n, p, c, psi_gauche, psi_droit, b_psi)
        call linSolve("plu", m%l, A_psi, b_psi, psi)
    end subroutine

END MODULE donnees
