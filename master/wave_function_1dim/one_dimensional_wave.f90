program evolucion_temporal
    implicit none
    integer, parameter :: nx = 1001
    real :: time, dx, dt, inv2dx, A, sigma, x0, alpha, beta
    integer :: lastiter, iter, i
    real, dimension(nx) :: xh, phi, ppi, psi, ppi1, psi1, ppi2, psi2, ppi3, psi3

    ! Uso de la función para borrar el contenido del archivo por si tenía algo
    call borrar_contenido('phi.txt')

    ! DECLARAMOS NUESTRAS VARIABLES
    nx = 1001  ! número de puntos total
    time = 0.0
    dx = 0.002
    dt = 0.0005
    inv2dx = 1.0 / (2.0 * dx)
    A = 1.0
    sigma = 0.1
    x0 = 0.0
    alpha = 1.0
    beta = 1.0
    lastiter = 3000

    ! Inicializamos las funciones: de esta forma las funciones van de 0 a nx-1, con nx elementos
    xh = 0.0
    phi = 0.0
    ppi = 0.0
    psi = 0.0
    ppi1 = 0.0
    psi1 = 0.0
    ppi2 = 0.0
    psi2 = 0.0
    ppi3 = 0.0
    psi3 = 0.0

    print *, '    |-------------------------|'
    print *, '    |  STARTING 1D EVOLUTION  |'
    print *, '    |-------------------------|'

    ! DEFINIMOS LA MALLA ESPACIAL, que va de -1 a 1 con un salto de 0.002. Lo puedes mostrar para estar seguro
    do i = 1, nx
        xh(i) = -1.0 + (i-1) * dx
    end do

    ! DATOS INICIALES del CAMPO ESCALAR
    do i = 1, nx
        phi(i) = A * exp(-(xh(i) - x0)**2 / sigma**2)
        psi(i) = -2.0 * (xh(i) - x0) * phi(i) / sigma**2
    end do

    print *, ' *****************************************'
    print *, '  Starting Time Evolution'
    print *, ' *****************************************'

    open(1, file='phi.txt', status='unknown')
    do i = 1, nx
        write(1, *) time, xh(i), phi(i)
    end do
    write(1, *)

    ! EMPEZAMOS EL ALGORITMO PARA LA EVOLUCION TEMPORAL
    do iter = 1, lastiter
        time = time + dt

        ! RUNGE-KUTTA
        ! PRIMER PASO ppi1 y psi1, son las estrellas
        do i = 2, nx-1
            ppi1(i) = ppi(i) + dt * inv2dx * (alpha * (psi(i + 1) - psi(i - 1)) + beta * (ppi(i + 1) - ppi(i - 1)))
            psi1(i) = psi(i) + dt * inv2dx * (alpha * (ppi(i + 1) - ppi(i - 1)) + beta * (psi(i + 1) - psi(i - 1)))
        end do
        ! CONDICIONES DE FRONTERA ppi1 y psi1
        ppi1(1) = ppi1(3) + (xh(1) - xh(3)) / (xh(2) - xh(3)) * (ppi1(2) - ppi1(3))
        psi1(1) = ppi1(1)
        ppi1(nx) = ppi1(nx - 2) + (xh(nx) - xh(nx - 2)) / (xh(nx - 1) - xh(nx - 2)) * (ppi1(nx - 1) - ppi1(nx - 2))
        psi1(nx) = -ppi1(nx)

        ! SEGUNDO PASO para las dos estrellas
        do i = 2, nx-1
            ppi2(i) = 0.75 * ppi(i) + 0.25 * ppi1(i) + 0.25 * dt * inv2dx * (alpha * (psi1(i + 1) - psi1(i - 1)) + beta * (ppi1(i + 1) - ppi1(i - 1)))
            psi2(i) = 0.75 * psi(i) + 0.25 * psi1(i) + 0.25 * dt * inv2dx * (alpha * (ppi1(i + 1) - ppi1(i - 1)) + beta * (psi1(i + 1) - psi1(i - 1)))
        end do
        ! CONDICIONES DE FRONTERA para las dos estrellas
        ppi2(1) = ppi2(3) + (xh(1) - xh(3)) / (xh(2) - xh(3)) * (ppi2(2) - ppi2(3))
        psi2(1) = ppi2(1)
        ppi2(nx) = ppi2(nx - 2) + (xh(nx) - xh(nx - 2)) / (xh(nx - 1) - xh(nx - 2)) * (ppi1(nx - 1) - ppi1(nx - 2))
        psi2(nx) = -ppi2(nx)

        ! TERCER PASO
        do i = 2, nx-1
            ppi3(i) = (1.0/3.0) * ppi(i) + (2.0/3.0) * ppi2(i) + (2.0/3.0) * dt * inv2dx * (alpha * (psi2(i + 1) - psi2(i - 1)) + beta * (ppi2(i + 1) - ppi2(i - 1)))
            psi3(i) = (1.0/3.0) * psi(i) + (2.0/3.0) * psi2(i) + (2.0/3.0) * dt * inv2dx * (alpha * (ppi2(i + 1) - ppi2(i - 1)) + beta * (psi2(i + 1) - psi2(i - 1)))
        end do
        ppi3(1) = ppi3(3) + (xh(1) - xh(3)) / (xh(2) - xh(3)) * (ppi3(2) - ppi3(3))
        psi3(1) = ppi3(1)
        ppi3(nx) = ppi3(nx - 2) + (xh(nx) - xh(nx - 2)) / (xh(nx - 1) - xh(nx - 2)) * (ppi3(nx - 1) - ppi3(nx - 2))
        psi3(nx) = -ppi3(nx)

        ! La solución final no será más que:
        phi = phi + dt * (alpha * ppi3 + beta * psi3)
        ppi = ppi3
        psi = psi3

        if (mod(iter, 50) == 0) then
            open(1, file='phi.txt', status='unknown', position='append')
            do i = 1, nx
                write(1, *) time, xh(i), phi(i)
            end do
            write(1, *)
            print *, 'Output time,it: ', time, iter
        end if
    end do

    ! Cargar los datos del archivo
    open(1, file='phi.txt', status='old')
    do i = 1, nx
        read(1, *) time, xh(i), phi(i)
    end do

    ! Separar los datos en diferentes variables
    tiempo = time
    posicion = xh
    funcion = phi

    ! Crear una figura 3D
    call scatter_plot(tiempo, posicion, funcion)

contains

    subroutine borrar_contenido(nombre_archivo)
        character(len=*), intent(in) :: nombre_archivo
        integer :: unidad

        open(newunit=unidad, file=nombre_archivo, status='unknown')
        close(unidad)
    end subroutine borrar_contenido

    subroutine scatter_plot(tiempo, posicion, funcion)
        real, dimension(:), intent(in) :: tiempo, posicion, funcion
        integer :: nx
        integer :: i
        real :: time, xh, phi
        character(len=10) :: color

        ! Crear una figura 3D
        call plot3d_init()

        ! H
