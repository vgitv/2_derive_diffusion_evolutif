PROGRAM main

    use math
    use maillage
    use donnees
    use dd

    implicit none

    integer :: n, m
    real(rp) :: xmin, xmax, Tf

    real(rp), dimension(:), allocatable :: x2, t
    real(rp) :: cfl
    real(rp), dimension(:), allocatable :: approx_p, approx_n, approx_psi

    open(unit = 1, file = "entrees/constantes")
    read (1, *) n
    read (1, *) xmin
    read (1, *) xmax
    read (1, *) m
    read (1, *) Tf
    close(1)

    allocate(x2(n + 1), t(m))
    x2 = linspace(xmin, xmax, n + 1)
    call newMesh(x2, maill)

    t = linspace(0.0_rp, Tf, m)

    allocate(approx_n(maill%l), approx_p(maill%l), approx_psi(maill%l))



    ! -------------------------------------------------------------------------------------------------------
    ! Schéma explicite
    ! -------------------------------------------------------------------------------------------------------
    call vf_dd(t, maill, B2, c_dopage, n0, p0, n_l, n_r, p_l, p_r, psi_l, psi_r, &
        approx_n, approx_p, approx_psi)

    call saveSol(maill%x, (/ n_l(maill%x(1)), approx_n, n_r(maill%x(maill%l + 2)) /), "sorties/approx_n.dat")
    call saveSol(maill%x, (/ p_l(maill%x(1)), approx_p, p_r(maill%x(maill%l + 2)) /), "sorties/approx_p.dat")
    call saveSol(maill%x, &
        (/ psi_l(maill%x(1)), approx_psi, psi_r(maill%x(maill%l + 2)) /), &
        "sorties/approx_psi.dat")



    ! *******************************************************************************************************
    deallocate(approx_n, approx_p, approx_psi)
    call rmMesh(maill)
    deallocate(x2)

END PROGRAM main
