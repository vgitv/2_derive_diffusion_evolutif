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
    ! Pour n et p données (temps fixé), calcule psi au même temps
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
    ! t : temps courant
    ! dt : pas de temps
    ! m : maillage
    ! n : approx de n à un temps donné
    ! n_l : condition au bord gauche
    ! n_r : condition au bord droit
    ! psi : approx de psi au même temps
    ! psi_l : condition au bord gauche
    ! psi_r : condition au bord droit
    ! B : schéma
    ! n_suiv : approx de n au temps suivant (return)
    subroutine iter_n(t, dt, m, n, n_r, n_l, psi, psi_l, psi_r, B, n_suiv)
        ! paramètres
        real(rp), intent(in) :: t
        real(rp), intent(in) :: dt
        type(Mesh), intent(in) :: m
        real(rp), dimension(:), intent(in) :: n
        real(rp), external :: n_r, n_l
        real(rp), dimension(:), intent(in) :: psi
        real(rp), external :: psi_l, psi_r
        real(rp), external :: B
        real(rp), dimension(:), intent(out) :: n_suiv

        ! variables locales
        integer :: i
        real(rp) :: flux1, flux2

        ! premier élément
        flux1 = flux_G(m%h2(1), psi_l(t), psi(1), n_l(t), n(1), B)
        flux2 = flux_G(m%h2(2), psi(1), psi(2), n(1), n(2), B)
        n_suiv(1) = (flux1 - flux2) * dt / m%h(1) + n(1)

        do i = 2, m%l - 1
            !flux1 = flux_G(m%h2(i), psi(i - 1), psi(i), n(i - 1), n(i), B)
            flux1 = flux2
            flux2 = flux_G(m%h2(i + 1), psi(i), psi(i + 1), n(i), n(i + 1), B)
            n_suiv(i) = (flux1 - flux2) * dt / m%h(i) + n(i)
        end do

        ! dernier élément
        flux1 = flux_G(m%h2(m%l), psi(m%l - 1), psi(m%l), n(m%l - 1), n(m%l), B)
        flux2 = flux_G(m%h2(m%l + 1), psi(m%l), psi_r(t), n(m%l), n_r(t), B)
        n_suiv(m%l) = (flux1 - flux2) * dt / m%h(i) + n(i)
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Calcul le vecteur p au pas de temps suivant.
    ! -------------------------------------------------------------------------------------------------------
    ! t : temps courant
    ! dt : pas de temps
    ! m : maillage
    ! p : approx de n à un temps donné
    ! p_l : condition au bord gauche
    ! p_r : condition au bord droit
    ! psi : approx de psi au même temps
    ! psi_l : condition au bord gauche
    ! psi_r : condition au bord droit
    ! B : schéma
    ! p_suiv : approx de n au temps suivant (return)
    subroutine iter_p(t, dt, m, p, p_r, p_l, psi, psi_l, psi_r, B, p_suiv)
        ! paramètres
        real(rp), intent(in) :: t
        real(rp), intent(in) :: dt
        type(Mesh), intent(in) :: m
        real(rp), dimension(:), intent(in) :: p
        real(rp), external :: p_r, p_l
        real(rp), dimension(:), intent(in) :: psi
        real(rp), external :: psi_l, psi_r
        real(rp), external :: B
        real(rp), dimension(:), intent(out) :: p_suiv

        ! variables locales
        integer :: i
        real(rp) :: flux1, flux2

        ! premier élément
        flux1 = flux_H(m%h2(1), psi_l(t), psi(1), p_l(t), p(1), B)
        flux2 = flux_H(m%h2(2), psi(1), psi(2), p(1), p(2), B)
        p_suiv(1) = (flux1 - flux2) * dt / m%h(1) + p(1)

        do i = 2, m%l - 1
            !flux1 = flux_H(m%h2(i), psi(i - 1), psi(i), p(i - 1), p(i), B)
            flux1 = flux2
            flux2 = flux_H(m%h2(i + 1), psi(i), psi(i + 1), p(i), p(i + 1), B)
            p_suiv(i) = (flux1 - flux2) * dt / m%h(i) + p(i)
        end do

        ! dernier élément
        flux1 = flux_H(m%h2(m%l), psi(m%l - 1), psi(m%l), p(m%l - 1), p(m%l), B)
        flux2 = flux_H(m%h2(m%l + 1), psi(m%l), psi_r(t), p(m%l), p_r(t), B)
        p_suiv(m%l) = (flux1 - flux2) * dt / m%h(i) + p(i)
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! schéma volumes finis explicite dérive diffusion évolutif
    ! -------------------------------------------------------------------------------------------------------
    ! IN :
    ! t : maillage temps
    ! m : maillage espace
    ! B : schéma parmi B1 (centré), B2 (décentré amont) et B3 (schéma de Scharfetter-Gummel)
    ! c : dopage sur les pts du maillage
    ! n0 : condition initiale (CI) pour n
    ! p0 : CI pour p
    ! n_l : condition limite (CL) left pour n
    ! n_r : CL right pour n
    ! p_l, p_r, psi_l, psi_r : idem
    ! OUT :
    ! n : électrons
    ! p : trous
    ! psi : potentiel
    subroutine vf_explicite_dd(t, m, B, c, n0, p0, n_l, n_r, p_l, p_r, psi_l, psi_r, n, p, psi)
        ! paramètres
        real(rp), dimension(:), intent(in) :: t
        type(Mesh), intent(in) :: m
        real(rp), external :: B, c, n0, p0, n_l, n_r, p_l, p_r, psi_l, psi_r
        real(rp), dimension(:), intent(out) :: n, p, psi

        ! variables locales
        real(rp) :: dt
        integer :: i
        real(rp), dimension(m%l) :: n_suiv, p_suiv, psi_suiv, ci

        dt = t(2) - t(1)

        ! initialisation n0 p0
        do i = 1, m%l
            n(i) = n0(m%x(i + 1))
            p(i) = p0(m%x(i + 1))
            ci(i) = c(m%x(i + 1))
        end do

        ! calcul psi0 avec subroutine potentiel
        call potentiel(m, n, p, ci, psi_l(0), psi_r(0), psi)

        do i = 1, size(t)
            ! calcul itéré n suivant
            call iter_n(t(i), dt, m, n, n_r, n_l, psi, psi_l, psi_r, B, n_suiv)

            ! calcul itéré p suivant
            call iter_p(t(i), dt, m, p, p_r, p_l, psi, psi_l, psi_r, B, p_suiv)

            ! déduction itéré psi suivant
            call potentiel(m, n_suiv, p_suiv, ci, psi_l(t(i)), psi_r(t(i)), psi_suiv)

            ! remplacement itérés
            n = n_suiv
            p = p_suiv
            psi = psi_suiv
        end do
    end subroutine



    ! =======================================================================================================
    ! SCHÉMA SEMI IMPLICITE
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! Matrice A_N (A semi implicite)
    ! -------------------------------------------------------------------------------------------------------
    ! t_courant : temps donné
    ! dt : pas de temps
    ! m : maillage d'espace
    ! psi : vecteur itéré kième
    ! psi_l : fonction CL bord gauche
    ! psi_r : CL droit
    ! B : fonction schéma parmi B1, B2, B3
    ! A : matrice retournée
    subroutine build_A_si(t_courant, dt, m, psi, psi_l, psi_r, B, A)
        ! paramètres
        real(rp) :: t_courant
        real(rp) :: dt
        type(Mesh), intent(in) :: m
        real(rp), dimension(m%l), intent(in) :: psi
        real(rp), external :: psi_l, psi_r
        real(rp), external :: B
        real(rp), dimension(m%l, m%l), intent(out) :: A

        ! variables locales
        integer :: i

        A(1, 1) = m%h(1) / dt + B( psi(1)-psi(2) ) / m%h2(2) + B( psi(1)-psi_l(t_courant) ) / m%h2(1)
        A(1, 2) = -B( psi(2)-psi(1) ) / m%h2(2)
        A(2, 1) = -B( psi(1)-psi(2) ) / m%h2(2)
        do i = 2, m%l - 1
            A(i, i) = m%h(i) / dt + B( psi(i)-psi(i+1) ) / m%h2(i+1) + B( psi(i)-psi(i-1) ) / m%h2(i)
            A(i, i+1) = -B( psi(i+1)-psi(i) ) / m%h2(i+1)
            A(i+1, i) = -B( psi(i)-psi(i+1) ) / m%h2(i+1)
        end do
        A(m%l, m%l) = m%h(m%l) / dt + B( psi(m%l)-psi_r(t_courant) ) / m%h2(m%l+1) + &
            B( psi(m%l)-psi(m%l-1) ) / m%h2(m%l)
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Matrice b_N (b semi implicite)
    ! -------------------------------------------------------------------------------------------------------
    ! t_courant : temps donné
    ! dt : pas de temps
    ! m : maillage d'espace
    ! n : vecteur itéré kième
    ! n_l : fonction CL bord gauche
    ! n_r : CL droit
    ! B : fonction schéma parmi B1, B2, B3
    ! b_N : matrice retournée
    subroutine build_b_si(t_courant, dt, m, n, n_l, n_r, B, b_N)
        ! paramètres
        real(rp), intent(in) :: t_courant
        real(rp), intent(in) :: dt
        type(Mesh), intent(in) :: m
        real(rp), dimension(m%l), intent(in) :: n
        real(rp), external :: n_l, n_r
        real(rp), external :: B
        real(rp), dimension(m%l), intent(out) :: b_N

        ! variables locales
        integer :: i

        do i = 1, m%l
            !++!
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! schéma volumes finis semi implicite dérive diffusion évolutif
    ! -------------------------------------------------------------------------------------------------------
    ! IN :
    ! t : maillage temps
    ! m : maillage espace
    ! B : schéma parmi B1 (centré), B2 (décentré amont) et B3 (schéma de Scharfetter-Gummel)
    ! c : dopage sur les pts du maillage
    ! n0 : condition initiale (CI) pour n
    ! p0 : CI pour p
    ! n_l : condition limite (CL) left pour n
    ! n_r : CL right pour n
    ! p_l, p_r, psi_l, psi_r : idem
    ! OUT :
    ! n : électrons
    ! p : trous
    ! psi : potentiel
    subroutine vf_semi_implicite_dd(t, m, B, c, n0, p0, n_l, n_r, p_l, p_r, psi_l, psi_r, n, p, psi)
        ! paramètres
        real(rp), dimension(:), intent(in) :: t
        type(Mesh), intent(in) :: m
        real(rp), external :: B, c, n0, p0, n_l, n_r, p_l, p_r, psi_l, psi_r
        real(rp), dimension(:), intent(out) :: n, p, psi

        ! variables locales
        real(rp) :: dt
        integer :: i
        real(rp), dimension(m%l) :: n_suiv, p_suiv, psi_suiv, ci
        real(rp), dimension(m%l, m%l) :: A_N
        real(rp), dimension(m%l) :: b_N

        dt = t(2) - t(1)

        ! initialisation n0 p0
        do i = 1, m%l
            n(i) = n0(m%x(i + 1))
            p(i) = p0(m%x(i + 1))
            ci(i) = c(m%x(i + 1))
        end do

        ! calcul psi0 avec subroutine potentiel
        call potentiel(m, n, p, ci, psi_l(0), psi_r(0), psi)

        do i = 1, size(t)
            ! calcul itéré n suivant
            call iter_n(t(i), dt, m, n, n_r, n_l, psi, psi_l, psi_r, B, n_suiv)

            ! calcul itéré p suivant
            call iter_p(t(i), dt, m, p, p_r, p_l, psi, psi_l, psi_r, B, p_suiv)

            ! déduction itéré psi suivant
            call potentiel(m, n_suiv, p_suiv, ci, psi_l(t(i)), psi_r(t(i)), psi_suiv)

            ! remplacement itérés
            n = n_suiv
            p = p_suiv
            psi = psi_suiv
        end do
    end subroutine

END MODULE dd
