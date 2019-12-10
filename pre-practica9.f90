! ---------------------------------- Pre-pràctica 9 ------------------------------------- !
! Autor: Javier Rozalén Sarmiento
! Grup: B1B
! Data: 17/12/2019
!
! Funcionalitat: es resol l'equació de Poisson en 2-D per trobar el camp de temperatures
! en una geometria rectangular amb condicions de contorn de Dirichlet.
!
! Comentaris: els apartats 1 i 2 es resumeixen en la subrutina "solver", la funció "rho" i
! en els primers blocs del programa principal.
!
! Nota: el programa principal i les subrutines tenen indentacions tals que, en alguns
! editors de text (SublimeText, p.ex), es poden plegar i desplegar per facilitar la lectura.

program pre_practica9
    implicit none
    integer Nx,Ny,k,l,N
    double precision Lx,Ly,h,rho
    double precision, allocatable :: T_old(:,:),T_new(:,:)
    external rho
    common/cts/Lx,Ly,Nx,Ny

    Lx=32.5d0
    Ly=16.5d0
    h=0.5d0
    Nx=int(Lx/h)
    Ny=int(Ly/h)
    allocate(T_old(Nx,Ny))
    allocate(T_new(Nx,Ny))
    
    ! Condicions inicials
    do k=1,Nx
        do l=1,Ny
            T_old(k,l)=0.d0 
        enddo
    enddo
    do k=1,Nx
        T_old(1,k)=4.d0
        T_old(k,1)=3.36d0
        T_old(k,Ny)=23.1d0
        T_old(Nx,k)=25.d0
    enddo

    ! -------------------------------- Apartat 3 --------------------------------------- !
    open(11,file="aux1.dat")
    write(11,*) "# Niteracions, Temperatura"
    N=2000
    do k=1,N
        call solver(T_old,T_new,h,rho,0)
        write(11,*) k,T_new(int(18.d0/h),int(12.5d0/h)) ! index 0 (gnuplot)
        T_old=T_new
    enddo

    call write(11)

    ! Condicions inicials
    do k=1,Nx
        do l=1,Ny
            T_old(k,l)=0.d0 
        enddo
    enddo
    do k=1,Nx
        T_old(1,k)=4.d0
        T_old(k,1)=3.36d0
        T_old(k,Ny)=23.1d0
        T_old(Nx,k)=25.d0
    enddo
    
    do k=1,N
        call solver(T_old,T_new,h,rho,1)
        write(11,*) k,T_new(int(18.d0/h),int(12.5d0/h)) ! index 1 (gnuplot)
        T_old=T_new
    enddo
    close(11)
end program pre_practica9

! Subrutina Jacobi --> Calcula una iteració per resoldre l'equació de Poisson
subroutine solver(T_old,T_new,h,funci,icontrol)
    ! T_old --> Matriu amb els valors a recalcular
    ! T_new --> Matriu amb els valors nous calculats
    ! h --> Interval de particionat de la malla
    ! funci --> Funció igualada al laplacià de la funció incògnita
    ! icontrol --> Variable de control del mètode de resolució
    implicit none
    integer i,j,Nx,Ny,icontrol
    double precision T_old(Nx,Ny),T_new(Nx,Ny),h,funci,Lx,Ly
    common/cts/Lx,Ly,Nx,Ny

    if (icontrol.eq.0) then ! Jacobi
        do i=2,Nx-1
            do j=2,Ny-1
                T_new(i,j)=0.25d0*(T_old(i,j+1)+T_old(i,j-1)+T_old(i+1,j)+T_old(i-1,j)+(h**2.d0)*funci(i*h,j*h))
            enddo
        enddo
    else if (icontrol.eq.1) then ! Gauss-Seidel
        do i=2,Nx-1
            do j=2,Ny-1
                T_new(i,j)=0.25d0*(T_old(i,j+1)+T_old(i,j-1)+T_old(i+1,j)+T_old(i-1,j)+(h**2.d0)*funci(i*h,j*h))
                T_old(i,j)=T_new(i,j)
            enddo
        enddo
    else
        print*,"A viam, xato..."
    endif

    return 
end subroutine solver

! Subrutina write --> Escriu dues línies en blanc en un arxiu
subroutine write(arxiu)
    ! arxiu --> número de l'arxiu
    implicit none
    integer arxiu

    write(arxiu,*) ""
    write(arxiu,*) ""

    return
end subroutine

double precision function rho(x,y)
    implicit none
    double precision x,y,rho1,rho2,r

    r=dsqrt((x-7.d0)**2.d0+(y-8.d0)**2.d0)
    rho1=1.33d0*dexp(-((r-4.d0)/0.4d0)**2.d0)
    if ((x.ge.20.d0).and.(x.le.24.d0).and.(y.ge.13.d0).and.(y.le.15.d0)) then
        rho2=1.3d0
    else
        rho2=0.d0
    endif
    rho=rho1+rho2

    return
end function rho