PROGRAM main

    use math
    use maillage
    use donnees
    use dd

    implicit none

    integer :: n, m
    real(rp) :: xmin, xmax, Tf

    real(rp), dimension(:), allocatable :: x, t

    open(unit = 1, file = "entrees/constantes")
    read (1, *) n
    read (1, *) xmin
    read (1, *) xmax
    read (1, *) m
    read (1, *) Tf
    close(1)

    allocate(x(n), t(m))
    x = linspace(xmin, xmax, n + 2)
    call newMesh(x, maill)

    t = linspace(0.0_rp, Tf, m)

    ! -------------------------------------------------------------------------------------------------------
    ! Schéma explicite
    ! -------------------------------------------------------------------------------------------------------

END PROGRAM main
