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
    integer Nx,Ny,k,l,m,counter,criteri
    double precision Lx,Ly,h,rho,t_ini(3),epsilon
    double precision, allocatable :: T_old(:,:),T_new(:,:),T_check(:,:)
    external rho
    common/cts/Lx,Ly,Nx,Ny

    Lx=32.5d0
    Ly=16.5d0
    h=0.25d0
    Nx=int(Lx/h)
    Ny=int(Ly/h)
    epsilon=0.0001d0
    allocate(T_old(Nx,Ny))
    allocate(T_new(Nx,Ny))
    allocate(T_check(Nx,Ny))
    

    ! -------------------------------- Apartat 3 --------------------------------------- !
    open(11,file="aux1.dat") ! Convergència de Jacobi
    open(12,file="aux2.dat") ! Convergència de Gauss-Seidel
    open(14,file="aux4.dat") ! Convergència Sobrerelaxació
    write(11,*) "# Niteracions, Temperatura"
    write(12,*) "# Niteracions, Temperatura"
    t_ini=(/6.d0,19.d0,320.d0/) ! Vector amb les configuracions inicials de temperatura

    ! Convergència de Jacobi
    do m=1,3
        call ci(T_old,t_ini(m))
        T_new=T_old
        counter=1
        criteri=0
        do while (criteri.ne.1)
            call solver(T_old,T_new,h,rho,0)
            write(11,*) counter,T_new(int(18.d0/h),int(12.5d0/h)) ! index m-1 (gnuplot) 
            counter=counter+1
            call check(T_new,T_old,criteri,epsilon)
            T_old=T_new
        enddo
        call write(11)
    enddo
 
    ! Convergència de Gauss-Seidel
    do m=1,3
        call ci(T_old,t_ini(m))
        T_new=T_old
        counter=1
        criteri=0
        do while (criteri.ne.1)
            T_check=T_old
            call solver(T_old,T_new,h,rho,1)
            write(12,*) counter,T_new(int(18.d0/h),int(12.5d0/h)) ! index m-1 (gnuplot) 
            counter=counter+1
            call check(T_new,T_check,criteri,epsilon)
            T_old=T_new
        enddo
        call write(12)
    enddo
 
    ! Convergència de Sobrerelaxació
    do m=1,3
        call ci(T_old,t_ini(m))
        T_new=T_old
        counter=1
        criteri=0
        do while (criteri.ne.1)
            T_check=T_old
            call solver(T_old,T_new,h,rho,2)
            write(14,*) counter,T_new(int(18.d0/h),int(12.5d0/h)) ! index m-1 (gnuplot) 
            counter=counter+1
            call check(T_new,T_check,criteri,epsilon)
            T_old=T_new
        enddo
        call write(14)
    enddo

    ! Escriptura en arxiu pel mapejat gràfic 2D i 3D
    open(13,file="aux3.dat")
    write(13,*) "#X, Y, Temperatura"

    ! Cridem el nostre solver preferit per fer els mapes
    call ci(T_old,6.d0)
    T_new=T_old
    counter=1
    criteri=0
    do while (criteri.ne.1)
        T_check=T_old
        call solver(T_old,T_new,h,rho,2) 
        counter=counter+1
        call check(T_new,T_check,criteri,epsilon)
        T_old=T_new
    enddo

    ! Escrivim el resultat en un arxiu
    do k=1,Nx
        do l=1,Ny
            write(13,*) k*h,l*h,T_new(k,l)
        enddo
        write(13,*) ""
    enddo

    close(11)
    close(12)
    close(13)
    close(14)
    
    ! -------------------------------- Apartat 4 --------------------------------------- !
    ! Els programes "script_map2D.gnu" i "script_map3D.gnu" s'encarreguen de crear els mapes
    ! amb escala de colors a partir de l'arxiu en què acabem d'escriure.
end program pre_practica9

! Subrutina solver --> Calcula una iteració per resoldre l'equació de Poisson
subroutine solver(T_old,T_new,h,funci,icontrol)
    ! T_old --> Matriu amb els valors a recalcular
    ! T_new --> Matriu amb els valors nous calculats
    ! h --> Interval de particionat de la malla
    ! funci --> Funció igualada al laplacià de la funció incògnita
    ! icontrol --> Variable de control del mètode de resolució
        ! 0 --> Jacobi
        ! 1 --> Gauss-Seidel
        ! 2 --> Sobrerelaxació (pròximament)
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
    else if (icontrol.eq.2) then ! Sobrerelaxació
        do i=2,Nx-1
            do j=2,Ny-1
                T_new(i,j)=T_old(i,j)+0.25d0*1.73d0*(T_old(i,j+1)+T_old(i,j-1)+T_old(i+1,j)+T_old(i-1,j)-4.d0*T_old(i,j))
                T_new(i,j)=T_new(i,j)+0.25d0*1.73d0*(h**2.d0)*funci(i*h,j*h) ! Aquesta línea és la continuació de la de dalt
                T_old(i,j)=T_new(i,j)
            enddo
        enddo
    endif

    return 
end subroutine solver

! Subrutina ci --> Retorna una matriu amb les condicions inicials per la temperatura
subroutine ci(T_old,t_int)
    ! T_old --> matriu amb les condicions inicials
    ! t_int --> temperatura (homogènia) a l'interior
    implicit none
    double precision T_old(Nx,Ny),t_int
    double precision Lx,Ly
    integer Nx,Ny
    integer k,l
    common/cts/Lx,Ly,Nx,Ny

    ! Especifiquem la temperatura de l'interior
    do k=2,Nx-1
        do l=2,Ny-1
            T_old(k,l)=t_int 
        enddo
    enddo
    ! Especifiquem la temperatura de les parets horitzontals
    do k=1,Nx
        T_old(k,1)=3.36d0
        T_old(k,Ny)=23.1d0
    enddo
    ! Especifiquem la temperatura de les parets verticals
    do l=1,Ny
        T_old(1,l)=4.d0
        T_old(Nx,l)=25.d0
    enddo

    return
end subroutine

! Subrutina check --> Retorna 1 o 0 en funció de si es compleix cert criteri o no
subroutine check(T_1,T_2,criteri,epsilon)
    ! T_1,T_2 --> Matrius a comparar
    ! Criteri --> És 0 si no es compleix el criteri, i 1 si ho fa 
    ! epsilon --> Variable de control del criteri
    implicit none
    double precision T_1(Nx,Ny),T_2(Nx,Ny),Lx,Ly,epsilon
    integer i,j,Nx,Ny,criteri
    common/cts/Lx,Ly,Nx,Ny

    criteri=1
    do i=2,Nx-1
        do j=2,Ny-1
            if (abs(T_1(i,j)-T_2(i,j)).gt.epsilon) then
                criteri=0
            endif
        enddo
    enddo

    return
end subroutine check

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
