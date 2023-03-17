! ===========================================================================================================
! Module pour schéma résolution dérive diffusion.
! ===========================================================================================================

MODULE dd

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
    ! m : maillage d'espace
    ! A : matrice retournée
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
    ! flux numériques
    ! -------------------------------------------------------------------------------------------------------
    ! schéma centré
    function B1(x) result(B)
        real(rp), intent(in) :: x
        real(rp) :: B
        B = 1.0_rp - 0.5_rp * x
    end function

    ! schéma décentré amont
    function B2(x) result(B)
        real(rp), intent(in) :: x
        real(rp) :: B
        B = 1.0_rp - min(x, 0.0_rp)
    end function

    ! schéma Scharfetter-Gummel
    function B3(x) result(B)
        real(rp), intent(in) :: x
        real(rp) :: B
        B = x / (exp(x) - 1.0_rp)
    end function

    ! flux pour n
    function flux_G(h, psi_gauche, psi_droit, n_gauche, n_droit, B) result(flux)
        real(rp), intent(in) :: h, psi_gauche, psi_droit, n_gauche, n_droit
        real(rp), external :: B
        real(rp) :: flux
        flux = (B(psi_gauche - psi_droit) * n_gauche - B(psi_droit - psi_gauche) * n_droit) / h
    end function

    ! flux pour p
    function flux_H(h, psi_gauche, psi_droit, p_gauche, p_droit, B) result(flux)
        real(rp), intent(in) :: h, psi_gauche, psi_droit, p_gauche, p_droit
        real(rp), external :: B
        real(rp) :: flux
        flux = (B(psi_droit - psi_gauche) * p_gauche - B(psi_gauche - psi_droit) * p_droit) / h
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Calcul le vecteur n au pas de temps suivant.
    ! -------------------------------------------------------------------------------------------------------
    subroutine iter_n(m, n, B, n_suiv)
        ! paramètres
        type(Mesh), intent(in) :: m
        real(rp), dimension(:), intent(in) :: n
        real(rp), external :: B
        real(rp), dimension(:), intent(out) :: n_suiv

        ! variables locales
        !++!
    end subroutine

END MODULE dd
