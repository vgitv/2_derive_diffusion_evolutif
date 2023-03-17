MODULE donnees

    use math
    use maillage

    implicit none

contains

    ! =======================================================================================================
    ! Schéma explicite
    ! =======================================================================================================
    ! -------------------------------------------------------------------------------------------------------
    ! Construction matrice A_psi : A_psi psi = b_psi. La sol donne psi pour n et p donnés (temps fixé)
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



    ! -------------------------------------------------------------------------------------------------------
    ! Calcul le vecteur n au pas de temps suivant.
    ! -------------------------------------------------------------------------------------------------------
    subroutine iter_n(m, n, n_suiv)
        ! paramètres
        type(Mesh), intent(in) :: m
        real(rp), dimension(:), intent(in) :: n
        real(rp), dimension(:), intent(out) :: n_suiv

        ! variables locales
        !++!
    end subroutine

END MODULE donnees
