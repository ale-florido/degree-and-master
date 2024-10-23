 program autogravitante

      implicit none

      !definimos las variables

      INTEGER, PARAMETER :: N=10000, M=120000,final=20  !10.000, 360.000
      INTEGER :: i, j,iter
      REAL, PARAMETER :: PI = 4.0 * ATAN(1.0)
      REAL (KIND=8) :: B,C,eta,eta2,sigma,sigma2,dr,dt,r=1000.d0,inc,mass=1.d0,t_min=0.d0,t_max=3000.d0,omega=1.d0,lambda=0.d0
      REAL(KIND=8)::particulas1,particulas2,l2norm_ligadura
      REAL(KIND=8),DIMENSION(-1:N) :: rh, phi1, pi1, psi1, phi2, pi2, psi2, phi3, pi3, psi3, phi4, pi4, psi4, a, alpha,a_old,xh
      REAL(KIND=8),DIMENSION(0:M) :: th
      REAL(KIND=8),DIMENSION(-1:N) :: a_1, a_2, alpha_1, alpha_2, a_aux, ligadura, a_aux2,a_aux3
      REAL(KIND=8),DIMENSION(-1:N) :: phi1_1, phi1_2, phi2_3, psi1_1, psi1_2, psi1_3, pi1_1, pi1_2, pi1_3
      REAL(KIND=8),DIMENSION(-1:N) :: phi2_1, phi2_2, phi1_3, psi2_1, psi2_2, psi2_3, pi2_1, pi2_2, pi2_3
      REAL(KIND=8),DIMENSION(-1:N) :: phi3_1, phi3_2, phi3_3, psi3_1, psi3_2, psi3_3, pi3_1, pi3_2, pi3_3
      REAL(KIND=8),DIMENSION(-1:N) :: phi4_1, phi4_2, phi4_3, psi4_1, psi4_2, psi4_3, pi4_1, pi4_2, pi4_3

      !definimos nuestras malla de puntos radiales
      dr=r/N
      WRITE(*,*) 'dr', dr
      inc=0.d0
      i=-1
      do while (i<=N)
          rh(i)=inc-dr/2.d0
          inc=inc+dr
          i=i+1
      end do

      !creamos un vector temporal
      dt=(t_max-t_min)/M
      WRITE(*,*) 'dt', dt
      inc=0.d0
      i=0
      do while (i<=M)
          th(i)=inc
          inc=inc+dt
          i=i+1
      end do


      !condiciones iniciales
      eta=0.2d0 !relación entre los dos campos, campo 2 entre 1 ('excitado entre gs')
      !sqrt(3)
      B=0.005d0  !0.197d0/sqrt(4.d0*PI)
      C=B*eta

      sigma=80.d0
      !Cambiar valor eta2 para que no sea una superposición sencilla entre los campos
      eta2=90.d0/80.d0
      sigma2=eta2*sigma
      

      !campos en t=0 suponiendo minkowski para las partes imaginarias

      open(2, file='data/campo.txt',STATUS='REPLACE',ACTION= 'WRITE',POSITION='append')
      i=0
      do while (i<=N)
          phi1(i)=B*exp(-(rh(i)/sigma)**(2.d0))
          psi1(i)=-(2.d0*rh(i)/sigma**(2.d0))*phi1(i)
          pi1(i)=0.d0
          phi2(i)=0.d0
          psi2(i)=0.d0
          !PPI2 PORQUE USAMOS MINKOWSKI?
          pi2(i)=-omega*phi1(i)  

          phi3(i)=C*exp(-(rh(i)/sigma2)**(2.d0))
          psi3(i)=-(2.d0*rh(i)/sigma2**(2.d0))*phi3(i)
          pi3(i)=0.d0
          phi4(i)=0.d0
          psi4(i)=0.d0
          !PPI2 PORQUE USAMOS MINKOWSKI?
          pi4(i)=-omega*phi3(i) 

          !Vamos añadiendo todos los valores del campo hasta que lleguemos a r=200, que sólo añadiremos algunos
          IF (rh(i)<200) THEN
            write(2,*) rh(i), th(0), phi1(i), phi2(i), phi3(i),phi4(i)
          END IF

          IF (rh(i)>200) THEN
            IF (MOD(i,100).eq.0) THEN
            write(2,*) rh(i), th(0), phi1(i), phi2(i), phi3(i),phi4(i)
            END IF
          END IF
          i=i+1
      end do
      !En Fortran, se define el -1 como un ghost point para ayudarnos con las CC
      phi1(-1)=phi1(0)
      psi1(-1)=-psi1(0)
      pi1(-1)=pi1(0)
      phi2(-1)=phi2(0)
      psi2(-1)=-psi2(0)
      pi2(-1)=pi2(0)

      phi3(-1)=phi3(0)
      psi3(-1)=-psi3(0)
      pi3(-1)=pi3(0)
      phi4(-1)=phi4(0)
      psi4(-1)=-psi4(0)
      pi4(-1)=pi4(0)
      close(2)

      !Runge-Kutta de tres iteraciones para obtener los coefientes m�tricos


      !intentemos comenzar método iterativo para guardar los valores de a, alpha y pi2, entre otros
      !pasando así de minkowski al caso real
      iter=0
      !lo hace muy rápido, así que puse final=10 pq converge en 4 o 5 pasos. Ver que iría pasando para diversos valores
    
      do while (iter<=final)
      
      a_old=a

      !en el origen del dominio num�rico suponemos metrica plana
      a(-1)=1 !En este punto fantasma le ponemos también que verifique esta condición
      a(0)=1
      i=0 !integramos del origen hacia la frontera
     
      do while(i<=N-1)
        !multiplico por dos por como hemos definido la derivada espacial
        a_1(i)=a(i)+ dr*a(i)*((1.d0-a(i)**(2.d0))/(rh(i))+&
        (rh(i)/2.d0)*(psi1(i)**(2.d0)+psi2(i)**(2.d0)+pi1(i)**(2.d0)+pi2(i)**(2.d0)&
        +2.d0*0.5d0*(mass*a(i))**(2.d0)*(phi1(i)**(2.d0)+phi2(i)**(2.d0))+&
        psi3(i)**(2.d0)+psi4(i)**(2.d0)+pi3(i)**(2.d0)+pi4(i)**(2.d0)&
        +2.d0*0.5d0*(mass*a(i))**(2.d0)*(phi3(i)**(2.d0)+phi4(i)**(2.d0))+&
        (a(i))**(2.d0)*lambda*((phi1(i)**(2.d0)+phi2(i)**(2.d0))**(2.d0)&
        & +(phi3(i)**(2.d0)+phi4(i)**(2.d0))**(2.d0)) ))

        xh(i) = rh(i) + dr

        a_2(i)=a(i)+0.5d0*dr*(a(i)*((1.d0-a(i)**(2.d0))/(rh(i))+&
        (rh(i)/2.d0)*(psi1(i)**(2.d0)+psi2(i)**(2.d0)+pi1(i)**(2.d0)+pi2(i)**(2.d0)&
        +(mass*a(i))**(2.d0)*(phi1(i)**(2.d0)+phi2(i)**(2.d0))+&
        psi3(i)**(2.d0)+psi4(i)**(2.d0)+pi3(i)**(2.d0)+pi4(i)**(2.d0)&
        +(mass*a(i))**(2.d0)*(phi3(i)**(2.d0)+phi4(i)**(2.d0))+&
        (a(i))**(2.d0)*lambda*((phi1(i)**(2.d0)+phi2(i)**(2.d0))**(2.d0)&
        & +(phi3(i)**(2.d0)+phi4(i)**(2.d0))**(2.d0)) )) +&
        &(a_1(i)*((1.d0-a_1(i)**(2.d0))/(xh(i))&
        +(xh(i)/2.d0)*(psi1(i+1)**(2.d0)+psi2(i+1)**(2.d0)+pi1(i+1)**(2.d0)+pi2(i+1)**(2.d0)&
        +(mass*a_1(i))**(2.d0)*(phi1(i+1)**(2.d0)+phi2(i+1)**(2.d0))+psi3(i+1)**(2.d0)&
        +psi4(i+1)**(2.d0)+pi3(i+1)**(2.d0)+pi4(i+1)**(2.d0)&
        +(mass*a_1(i))**(2.d0)*(phi3(i+1)**(2.d0)+phi4(i+1)**(2.d0))+&
        (a_1(i))**(2.d0)*lambda*((phi1(i+1)**(2.d0)+phi2(i+1)**(2.d0))**(2.d0)&
        +(phi3(i+1)**(2.d0)+phi4(i+1)**(2.d0))**(2.d0)) ) )) )

        a(i+1)=a_2(i)

    i=i+1
  end do


      
      !calculamos el otro coeficiente m trico desde la frontera hacia el origen del dominio computacional
      alpha(N)=1.d0/a(N)  
      !alpha(N)=1.d0  
      i=N

      do while(i>0)
        alpha_1(i)=alpha(i)-dr*alpha(i)*((1.d0-a(i)**(2.d0))/(-rh(i))+&
        (rh(i)/2.d0)*(psi1(i)**(2.d0)+psi2(i)**(2.d0)+pi1(i)**(2.d0)+pi2(i)**(2.d0)&
        -(mass*a(i))**(2.d0)*(phi1(i)**(2.d0)&
        +phi2(i)**(2.d0))+psi3(i)**(2.d0)+psi4(i)**(2.d0)+pi3(i)**(2.d0)+pi4(i)**(2.d0)&
        -(mass*a(i))**(2.d0)*(phi3(i)**(2.d0)+phi4(i)**(2.d0))-&
        (a(i))**(2.d0)*lambda*((phi1(i)**(2.d0)+phi2(i)**(2.d0))**(2.d0)&
        & +(phi3(i)**(2.d0)+phi4(i)**(2.d0))**(2.d0)) ))

        xh(i) = rh(i) - dr/1.d0
        !2\cdot0.5 (2 de la derivada, 0.5 del RK2)
        alpha_2(i)=alpha(i)-dr*(alpha(i)*((1.d0-a(i)**(2.d0))/(-2.d0*rh(i))+&
        (rh(i)/4.d0)*(psi1(i)**(2.d0)+psi2(i)**(2.d0)+pi1(i)**(2.d0)+pi2(i)**(2.d0)&
        -(mass*a(i))**(2.d0)*(phi1(i)**(2.d0)&
        +phi2(i)**(2.d0))+psi3(i)**(2.d0)+psi4(i)**(2.d0)+pi3(i)**(2.d0)+pi4(i)**(2.d0)&
        -(mass*a(i))**(2.d0)*(phi3(i)**(2.d0)+phi4(i)**(2.d0))-&
        (a(i))**(2.d0)*lambda*((phi1(i)**(2.d0)+phi2(i)**(2.d0))**(2.d0)&
        & +(phi3(i)**(2.d0)+phi4(i)**(2.d0))**(2.d0)) )) + &
        alpha_1(i)*((1.d0-a(i-1)**(2.d0))/(-2.d0*xh(i))+&
        (xh(i)/4.d0)*(psi1(i-1)**(2.d0)+psi2(i-1)**(2.d0)+pi1(i-1)**(2.d0)+pi2(i-1)**(2.d0)&
        -(mass*a(i-1))**(2.d0)*(phi1(i-1)**(2.d0)&
        +phi2(i-1)**(2.d0))+psi3(i-1)**(2.d0)+psi4(i-1)**(2.d0)+pi3(i-1)**(2.d0)+pi4(i-1)**(2.d0)&
        -(mass*a(i-1))**(2.d0)*(phi3(i-1)**(2.d0)+phi4(i-1)**(2.d0))-&
        (a(i-1))**(2.d0)*lambda*((phi1(i-1)**(2.d0)+phi2(i-1)**(2.d0))**(2.d0)&
        & +(phi3(i-1)**(2.d0)+phi4(i-1)**(2.d0))**(2.d0)) )))

        alpha(i-1)=alpha_2(i)
    i=i-1
  end do

    
      alpha(-1)=alpha(0)

      do i=0,N
            pi2(i) = -a(i)*phi1(i)/alpha(i)*omega  
            pi4(i)=-a(i)*phi3(i)/alpha(i)*omega  
      end do
      pi2(-1)=pi2(0)
      pi4(-1)=pi4(0)


      IF (MOD(iter,1).eq.0) THEN
        WRITE(*,*) 'variación a en un punto cualquiera =',a_old(N-1000)-a(N-1000)      
      END IF

     iter=iter+1
    end do

      !en este punto ya tenemos tenemos todas las condiciones iniciales, es decir,
      !los campos phi(r,0), psi(r,0) y pi(r,0). Adem�s, de los coefientes m�tricos
      !a(r,0) y alpha(r,0)

      !guardamos en fichero los valores iniciales de coeficentes metricos
      open(1, file='data/metrica.txt',STATUS='REPLACE',ACTION= 'WRITE',POSITION='append')

      i=0
      do while(i<=N)
        IF (rh(i)<300) THEN
          write(1,*) a(i), alpha(i), rh(i), th(0), 0, sqrt(rh(i)*alpha(i)*(alpha(i+1)-alpha(i))/(2.d0*dr)) 
        END IF

        IF (rh(i)>300) THEN
          IF (MOD(i,50).eq.0) THEN
          write(1,*) a(i), alpha(i), rh(i), th(0), 0, sqrt(rh(i)*alpha(i)*(alpha(i+1)-alpha(i))/(2.d0*dr)) 
          END IF
        END IF
        i=i+1
      end do
      close(1)


      open(3, file='data/max_min.txt',STATUS='REPLACE',ACTION= 'WRITE',POSITION='append')
      write(3,*) sqrt(phi1(1)*phi1(1) + phi2(1)*phi2(1)),sqrt(phi3(1)*phi3(1) + phi4(1)*phi4(1)), &
                minval(alpha), th(0),phi1(1),phi2(1), phi3(1),phi4(1)
      close(3)


      !EVOLUCIÓN TEMPORAL
      !ahora vamos a desarollar la evoluci�n, para ello usaremos un integrador Runge-Kutta de tres iteraciones




      j=1
      do while(j<=M)

      !primera iteraci�n
    IF (MOD(j,100).eq.0) THEN
        WRITE(*,*) 'Tiempo ', th(j)
    END IF

      i=0
      
      do while(i<=N-1)
          !CAMPO 1
            phi1_1(i)=phi1(i)+dt*((alpha(i)*pi1(i))/a(i))

            pi1_1(i)=pi1(i)+dt*(3.d0*(((rh(i+1)**(2.d0)*alpha(i+1)*psi1(i+1))/a(i+1)&
            -(rh(i-1)**(2.d0)*alpha(i-1)*psi1(i-1))/a(i-1))/(rh(i+1)**(3.d0)-rh(i-1)**(3.d0)))&
            -a(i)*(mass**(2.d0)+2.d0*lambda*(phi1(i)**(2.d0)+phi2(i)**(2.d0)))*alpha(i)*phi1(i))

            psi1_1(i)=psi1(i)+dt/(2.d0*dr)*(alpha(i+1)*pi1(i+1)/a(i+1)-alpha(i-1)*pi1(i-1)/a(i-1))
            
            phi2_1(i)=phi2(i)+dt*((alpha(i)*pi2(i))/a(i))

            pi2_1(i)=pi2(i)+dt*(3.d0*(((rh(i+1)**(2.d0)*alpha(i+1)*psi2(i+1))/a(i+1)&
            -(rh(i-1)**(2.d0)*alpha(i-1)*psi2(i-1))/a(i-1))/(rh(i+1)**(3.d0)-rh(i-1)**(3.d0)))&
            -a(i)*(mass**(2.d0)+2.d0*lambda*(phi1(i)**(2.d0)+phi2(i)**(2.d0)))*alpha(i)*phi2(i))

            psi2_1(i)=psi2(i)+dt/(2.d0*dr)*(alpha(i+1)*pi2(i+1)/a(i+1)-alpha(i-1)*pi2(i-1)/a(i-1))
            

            !CAMPO 2

            phi3_1(i)=phi3(i)+dt*((alpha(i)*pi3(i))/a(i))

            pi3_1(i)=pi3(i)+dt*(3.d0*(((rh(i+1)**(2.d0)*alpha(i+1)*psi3(i+1))/a(i+1)&
            -(rh(i-1)**(2.d0)*alpha(i-1)*psi3(i-1))/a(i-1))/(rh(i+1)**(3.d0)-rh(i-1)**(3.d0)))&
            -a(i)*(mass**(2.d0)+(2.d0)*lambda*(phi3(i)**(2.d0)+phi4(i)**(2.d0)))*alpha(i)*phi3(i))

            psi3_1(i)=psi3(i)+dt/(2.d0*dr)*(alpha(i+1)*pi3(i+1)/a(i+1)-alpha(i-1)*pi3(i-1)/a(i-1))
            

            phi4_1(i)=phi4(i)+dt*((alpha(i)*pi4(i))/a(i))

            pi4_1(i)=pi4(i)+dt*(3.d0*(((rh(i+1)**(2.d0)*alpha(i+1)*psi4(i+1))/a(i+1)&
            -(rh(i-1)**(2.d0)*alpha(i-1)*psi4(i-1))/a(i-1))/(rh(i+1)**(3.d0)-rh(i-1)**(3.d0)))&
            -a(i)*(mass**(2.d0)+(2.d0)*lambda*(phi3(i)**(2.d0)+phi4(i)**(2.d0)))*alpha(i)*phi4(i))

            psi4_1(i)=psi4(i)+dt/(2.d0*dr)*(alpha(i+1)*pi4(i+1)/a(i+1)-alpha(i-1)*pi4(i-1)/a(i-1))
            

            i=i+1
      end do

      !condiciones de frontera para ghost point - AHORA LO ENTIENDO XD
      phi1_1(-1)=phi1_1(0)
      psi1_1(-1)=-psi1_1(0)
      pi1_1(-1)=pi1_1(0)

      phi2_1(-1)=phi2_1(0)
      psi2_1(-1)=-psi2_1(0)
      pi2_1(-1)=pi2_1(0)

      phi3_1(-1)=phi3_1(0)
      psi3_1(-1)=-psi3_1(0)
      pi3_1(-1)=pi3_1(0)

      phi4_1(-1)=phi4_1(0)
      psi4_1(-1)=-psi4_1(0)
      pi4_1(-1)=pi4_1(0)

      !condiciones de frontera para la malla radial
      pi1_1(N)=pi1(N)-dt*((3.d0*pi1(N)-4.d0*pi1(N-1)+pi1(N-2))/(2.d0*dr)+pi1(N)/rh(N))
      phi1_1(N)=phi1(N)+dt*((alpha(N)*pi1(N))/a(N))
      psi1_1(N)=-pi1_1(N)-phi1_1(N)/rh(N)

      pi2_1(N)=pi2(N)-dt*((3.d0*pi2(N)-4.d0*pi2(N-1)+pi2(N-2))/(2.d0*dr)+pi2(N)/rh(N))
      phi2_1(N)=phi2(N)+dt*((alpha(N)*pi2(N))/a(N))
      psi2_1(N)=-pi2_1(N)-phi2_1(N)/rh(N)

      pi3_1(N)=pi3(N)-dt*((3.d0*pi3(N)-4.d0*pi3(N-1)+pi3(N-2))/(2.d0*dr)+pi3(N)/rh(N))
      phi3_1(N)=phi3(N)+dt*((alpha(N)*pi3(N))/a(N))
      psi3_1(N)=-pi3_1(N)-phi3_1(N)/rh(N)

      pi4_1(N)=pi4(N)-dt*((3.d0*pi4(N)-4.d0*pi4(N-1)+pi4(N-2))/(2.d0*dr)+pi4(N)/rh(N))
      phi4_1(N)=phi4(N)+dt*((alpha(N)*pi4(N))/a(N))
      psi4_1(N)=-pi4_1(N)-phi4_1(N)/rh(N)

      !segunda iteracion 

      i=0
      do while(i<=N-1)
            !CAMPO 1
            phi1_2(i)=(3.d0/4.d0)*phi1(i)+(1.d0/4.d0)*phi1_1(i)+(dt/4.d0)*((alpha(i)*pi1_1(i))/a(i))

            pi1_2(i)=(3.d0/4.d0)*pi1(i)+(1.d0/4.d0)*pi1_1(i)+(dt/4.d0)*(3.d0*(((rh(i+1)**(2.d0)*alpha(i+1)*psi1_1(i+1))/a(i+1)&
            -(rh(i-1)**(2.d0)*alpha(i-1)*psi1_1(i-1))/a(i-1))/(rh(i+1)**(3.d0)-rh(i-1)**(3.d0)))&
            -a(i)*(mass**(2.d0)+(2.d0)*lambda*(phi1_1(i)**(2.d0)+phi2_1(i)**(2.d0)))*alpha(i)*phi1_1(i))

            psi1_2(i)=(3.d0/4.d0)*psi1(i)+(1.d0/4.d0)*psi1_1(i)+&
            (dt/4.d0)/(2.d0*dr)*(alpha(i+1)*pi1_1(i+1)/a(i+1)-alpha(i-1)*pi1_1(i-1)/a(i-1))

            phi2_2(i)=(3.d0/4.d0)*phi2(i)+(1.d0/4.d0)*phi2_1(i)+(dt/4.d0)*((alpha(i)*pi2_1(i))/a(i))

            pi2_2(i)=(3.d0/4.d0)*pi2(i)+(1.d0/4.d0)*pi2_1(i)+(dt/4.d0)*(3.d0*(((rh(i+1)**(2.d0)*alpha(i+1)*psi2_1(i+1))/a(i+1)&
            -(rh(i-1)**(2.d0)*alpha(i-1)*psi2_1(i-1))/a(i-1))/(rh(i+1)**(3.d0)-rh(i-1)**(3.d0)))&
            -a(i)*(mass**(2.d0)+(2.d0)*lambda*(phi1_1(i)**(2.d0)+phi2_1(i)**(2.d0)))*alpha(i)*phi2_1(i))

            psi2_2(i)=(3.d0/4.d0)*psi2(i)+(1.d0/4.d0)*psi2_1(i)+&
            (dt/4.d0)/(2.d0*dr)*(alpha(i+1)*pi2_1(i+1)/a(i+1)-alpha(i-1)*pi2_1(i-1)/a(i-1))
            

            !CAMPO 2
            phi3_2(i)=(3.d0/4.d0)*phi3(i)+(1.d0/4.d0)*phi3_1(i)+(dt/4.d0)*((alpha(i)*pi3_1(i))/a(i))

            pi3_2(i)=(3.d0/4.d0)*pi3(i)+(1.d0/4.d0)*pi3_1(i)+(dt/4.d0)*(3.d0*(((rh(i+1)**(2.d0)*alpha(i+1)*psi3_1(i+1))/a(i+1)&
            -(rh(i-1)**(2.d0)*alpha(i-1)*psi3_1(i-1))/a(i-1))/(rh(i+1)**(3.d0)-rh(i-1)**(3.d0)))&
            -a(i)*(mass**(2.d0)+(2.d0)*lambda*(phi3_1(i)**(2.d0)+phi4_1(i)**(2.d0)))*alpha(i)*phi3_1(i))

            psi3_2(i)=(3.d0/4.d0)*psi3(i)+(1.d0/4.d0)*psi3_1(i)+&
            (dt/4.d0)/(2.d0*dr)*(alpha(i+1)*pi3_1(i+1)/a(i+1)-alpha(i-1)*pi3_1(i-1)/a(i-1))
            

            phi4_2(i)=(3.d0/4.d0)*phi4(i)+(1.d0/4.d0)*phi4_1(i)+(dt/4.d0)*((alpha(i)*pi4_1(i))/a(i))

            pi4_2(i)=(3.d0/4.d0)*pi4(i)+(1.d0/4.d0)*pi4_1(i)+(dt/4.d0)*(3.d0*(((rh(i+1)**(2.d0)*alpha(i+1)*psi4_1(i+1))/a(i+1)&
            -(rh(i-1)**(2.d0)*alpha(i-1)*psi4_1(i-1))/a(i-1))/(rh(i+1)**(3.d0)-rh(i-1)**(3.d0)))&
            -a(i)*(mass**(2.d0)+(2.d0)*lambda*(phi3_1(i)**(2.d0)+phi4_1(i)**(2.d0)))*alpha(i)*phi4_1(i))

            psi4_2(i)=(3.d0/4.d0)*psi4(i)+(1.d0/4.d0)*psi4_1(i)+&
            (dt/4.d0)/(2.d0*dr)*(alpha(i+1)*pi4_1(i+1)/a(i+1)-alpha(i-1)*pi4_1(i-1)/a(i-1))
            


            i=i+1
      end do

      !condiciones de frontera para ghost point
      phi1_2(-1)=phi1_2(0)
      psi1_2(-1)=-psi1_2(0)
      pi1_2(-1)=pi1_2(0)

      phi2_2(-1)=phi2_2(0)
      psi2_2(-1)=-psi2_2(0)
      pi2_2(-1)=pi2_2(0)

      phi3_2(-1)=phi3_2(0)
      psi3_2(-1)=-psi3_2(0)
      pi3_2(-1)=pi3_2(0)

      phi4_2(-1)=phi4_2(0)
      psi4_2(-1)=-psi4_2(0)
      pi4_2(-1)=pi4_2(0)

      !condiciones de frontera para la malla radial
      pi1_2(N)=pi1_1(N)-dt*((3.d0*pi1_1(N)-4.d0*pi1_1(N-1)+pi1_1(N-2))/(2.d0*dr)+pi1_1(N)/rh(N))
      phi1_2(N)=phi1_1(N)+dt*((alpha(N)*pi1_1(N))/a(N))
      psi1_2(N)=-pi1_2(N)-phi1_2(N)/rh(N)

      pi2_2(N)=pi2_1(N)-dt*((3.d0*pi2_1(N)-4.d0*pi2_1(N-1)+pi2_1(N-2))/(2.d0*dr)+pi2_1(N)/rh(N))
      phi2_2(N)=phi2_1(N)+dt*((alpha(N)*pi2_1(N))/a(N))
      psi2_2(N)=-pi2_2(N)-phi2_2(N)/rh(N)

      pi3_2(N)=pi3_1(N)-dt*((3.d0*pi3_1(N)-4.d0*pi3_1(N-1)+pi3_1(N-2))/(2.d0*dr)+pi3_1(N)/rh(N))
      phi3_2(N)=phi3_1(N)+dt*((alpha(N)*pi3_1(N))/a(N))
      psi3_2(N)=-pi3_2(N)-phi3_2(N)/rh(N)

      pi4_2(N)=pi4_1(N)-dt*((3.d0*pi4_1(N)-4.d0*pi4_1(N-1)+pi4_1(N-2))/(2.d0*dr)+pi4_1(N)/rh(N))
      phi4_2(N)=phi4_1(N)+dt*((alpha(N)*pi4_1(N))/a(N))
      psi4_2(N)=-pi4_2(N)-phi4_2(N)/rh(N)

      !tercera iteracion
      i=0
      do while(i<=N-1)
            phi1_3(i)=(1.d0/3.d0)*phi1(i)+(2.d0/3.d0)*phi1_2(i)+(2.d0/3.d0)*dt*((alpha(i)*pi1_2(i))/a(i))

            pi1_3(i)=(1.d0/3.d0)*pi1(i)+(2.d0/3.d0)*pi1_2(i)+(2.d0/3.d0)*dt*(3.d0*(((rh(i+1)**(2.d0)*alpha(i+1)*psi1_2(i+1))/a(i+1)&
            -(rh(i-1)**(2.d0)*alpha(i-1)*psi1_2(i-1))/a(i-1))/(rh(i+1)**(3.d0)-rh(i-1)**(3.d0)))&
            -a(i)*(mass**(2.d0)+(2.d0)*lambda*(phi1_2(i)**(2.d0)+phi2_2(i)**(2.d0)))*alpha(i)*phi1_2(i))

            psi1_3(i)=(1.d0/3.d0)*psi1(i)+(2.d0/3.d0)*psi1_2(i)+(2.d0/3.d0)*dt/&
            (2.d0*dr)*(alpha(i+1)*pi1_2(i+1)/a(i+1)-alpha(i-1)*pi1_2(i-1)/a(i-1))

            phi2_3(i)=(1.d0/3.d0)*phi2(i)+(2.d0/3.d0)*phi2_2(i)+(2.d0/3.d0)*dt*((alpha(i)*pi2_2(i))/a(i))

            pi2_3(i)=(1.d0/3.d0)*pi2(i)+(2.d0/3.d0)*pi2_2(i)+(2.d0/3.d0)*dt*(3.d0*(((rh(i+1)**(2.d0)*alpha(i+1)*psi2_2(i+1))/a(i+1)&
            -(rh(i-1)**(2.d0)*alpha(i-1)*psi2_2(i-1))/a(i-1))/(rh(i+1)**(3.d0)-rh(i-1)**(3.d0)))&
            -a(i)*(mass**(2.d0)+(2.d0)*lambda*(phi1_2(i)**(2.d0)+phi2_2(i)**(2.d0)))*alpha(i)*phi2_2(i))

            psi2_3(i)=(1.d0/3.d0)*psi2(i)+(2.d0/3.d0)*psi2_2(i)+(2.d0/3.d0)*dt/&
            (2.d0*dr)*(alpha(i+1)*pi2_2(i+1)/a(i+1)-alpha(i-1)*pi2_2(i-1)/a(i-1))
            
            phi3_3(i)=(1.d0/3.d0)*phi3(i)+(2.d0/3.d0)*phi3_2(i)+(2.d0/3.d0)*dt*((alpha(i)*pi3_2(i))/a(i))

            pi3_3(i)=(1.d0/3.d0)*pi3(i)+(2.d0/3.d0)*pi3_2(i)+(2.d0/3.d0)*dt*(3.d0*(((rh(i+1)**(2.d0)*alpha(i+1)*psi3_2(i+1))/a(i+1)&
            -(rh(i-1)**(2.d0)*alpha(i-1)*psi3_2(i-1))/a(i-1))/(rh(i+1)**(3.d0)-rh(i-1)**(3.d0)))&
            -a(i)*(mass**(2.d0)+(2.d0)*lambda*(phi3_2(i)**(2.d0)+phi4_2(i)**(2.d0)))*alpha(i)*phi3_2(i))

            psi3_3(i)=(1.d0/3.d0)*psi3(i)+(2.d0/3.d0)*psi3_2(i)+(2.d0/3.d0)*dt/&
            (2.d0*dr)*(alpha(i+1)*pi3_2(i+1)/a(i+1)-alpha(i-1)*pi3_2(i-1)/a(i-1))

            phi4_3(i)=(1.d0/3.d0)*phi4(i)+(2.d0/3.d0)*phi4_2(i)+(2.d0/3.d0)*dt*((alpha(i)*pi4_2(i))/a(i))

            pi4_3(i)=(1.d0/3.d0)*pi4(i)+(2.d0/3.d0)*pi4_2(i)+(2.d0/3.d0)*dt*(3.d0*(((rh(i+1)**(2.d0)*alpha(i+1)*psi4_2(i+1))/a(i+1)&
            -(rh(i-1)**(2.d0)*alpha(i-1)*psi4_2(i-1))/a(i-1))/(rh(i+1)**(3.d0)-rh(i-1)**(3.d0)))&
            -a(i)*(mass**(2.d0)+(2.d0)*lambda*(phi3_2(i)**(2.d0)+phi4_2(i)**(2.d0)))*alpha(i)*phi4_2(i))

            psi4_3(i)=(1.d0/3.d0)*psi4(i)+(2.d0/3.d0)*psi4_2(i)+(2.d0/3.d0)*dt/&
            (2.d0*dr)*(alpha(i+1)*pi4_2(i+1)/a(i+1)-alpha(i-1)*pi4_2(i-1)/a(i-1))

            !Y renombramos a nuestras variables originales

            phi1(i)=phi1_3(i)
            pi1(i)=pi1_3(i)
            psi1(i)=psi1_3(i)

            phi2(i)=phi2_3(i)
            pi2(i)=pi2_3(i)
            psi2(i)=psi2_3(i)

            phi3(i)=phi3_3(i)
            pi3(i)=pi3_3(i)
            psi3(i)=psi3_3(i)

            phi4(i)=phi4_3(i)
            pi4(i)=pi4_3(i)
            psi4(i)=psi4_3(i)

            i=i+1
      end do

      !condiciones frontera para la malla radial
      pi1(N)=pi1_2(N)-dt*((3.d0*pi1_2(N)-4.d0*pi1_2(N-1)+pi1_2(N-2))/(2.d0*dr)+pi1_2(N)/rh(N))
      phi1(N)=phi1_2(N)+dt*((alpha(N)*pi1_2(N))/a(N))
      psi1(N)=-pi1(N)-phi1(N)/rh(N)

      pi2(N)=pi2_2(N)-dt*((3.d0*pi2_2(N)-4.d0*pi2_2(N-1)+pi2_2(N-2))/(2.d0*dr)+pi2_2(N)/rh(N))
      phi2(N)=phi2_2(N)+dt*((alpha(N)*pi2_2(N))/a(N))
      psi2(N)=-pi2(N)-phi2(N)/rh(N)

      pi3(N)=pi3_2(N)-dt*((3.d0*pi3_2(N)-4.d0*pi3_2(N-1)+pi3_2(N-2))/(2.d0*dr)+pi3_2(N)/rh(N))
      phi3(N)=phi3_2(N)+dt*((alpha(N)*pi3_2(N))/a(N))
      psi3(N)=-pi3(N)-phi3(N)/rh(N)

      pi4(N)=pi4_2(N)-dt*((3.d0*pi4_2(N)-4.d0*pi4_2(N-1)+pi4_2(N-2))/(2.d0*dr)+pi4_2(N)/rh(N))
      phi4(N)=phi4_2(N)+dt*((alpha(N)*pi4_2(N))/a(N))
      psi4(N)=-pi4(N)-phi4(N)/rh(N)

      !Para el ghost point
      phi1(-1)=phi1(0)
      pi1(-1)=pi1(0)
      psi1(-1)=-psi1(0)

      phi2(-1)=phi2(0)
      pi2(-1)=pi2(0)
      psi2(-1)=-psi2(0)

      phi3(-1)=phi3(0)
      pi3(-1)=pi3(0)
      psi3(-1)=-psi3(0)

      phi4(-1)=phi4(0)
      pi4(-1)=pi4(0)
      psi4(-1)=-psi4(0)


      !definimos una variables auxiliares para verificar la ligadura de momento
      IF (j>=3) THEN
        i=0
        do while(i<=N)
        a_aux3(i)=a_aux2(i)
        i=i+1
        end do
        END IF


        IF (j>=2) THEN
        i=0
        do while(i<=N)
        a_aux2(i)=a_aux(i)
        i=i+1
        end do
        END IF

      i=0
      do while(i<=N)
         a_aux(i)=a(i)
         i=i+1
      end do

      !en este punto ya hemos actualizado todas las variables referentes al campo escalar
      !ahora calculamos como los nuevos valores del campo modifican los coeficientes metricos

      a(-1)=1
      a(0)=1
      i=0 !integramos del origen hacia la frontera
      do while(i<=N-1)
        a_1(i)=a(i)+ dr*a(i)*((1.d0-a(i)**(2.d0))/(rh(i))+&
        (rh(i)/2.d0)*(psi1(i)**(2.d0)+psi2(i)**(2.d0)+pi1(i)**(2.d0)+pi2(i)**(2.d0)&
        +2.d0*0.5d0*(mass*a(i))**(2.d0)*(phi1(i)**(2.d0)+phi2(i)**(2.d0))+&
        psi3(i)**(2.d0)+psi4(i)**(2.d0)+pi3(i)**(2.d0)+pi4(i)**(2.d0)&
        +2.d0*0.5d0*(mass*a(i))**(2.d0)*(phi3(i)**(2.d0)+phi4(i)**(2.d0))+&
        (a(i))**(2.d0)*lambda*((phi1(i)**(2.d0)+phi2(i)**(2.d0))**(2.d0)&
        & +(phi3(i)**(2.d0)+phi4(i)**(2.d0))**(2.d0)) ))

        xh(i) = rh(i) + dr

        a_2(i)=a(i)+dr*(a(i)*((1.d0-a(i)**(2.d0))/(2.d0*rh(i))+&
        (rh(i)/4.d0)*(psi1(i)**(2.d0)+psi2(i)**(2.d0)+pi1(i)**(2.d0)+pi2(i)**(2.d0)&
        +2.d0*0.5d0*(mass*a(i))**(2.d0)*(phi1(i)**(2.d0)+phi2(i)**(2.d0))+&
        psi3(i)**(2.d0)+psi4(i)**(2.d0)+pi3(i)**(2.d0)+pi4(i)**(2.d0)&
        +2.d0*0.5d0*(mass*a(i))**(2.d0)*(phi3(i)**(2.d0)+phi4(i)**(2.d0))+&
        (a(i))**(2.d0)*lambda*((phi1(i)**(2.d0)+phi2(i)**(2.d0))**(2.d0)&
        & +(phi3(i)**(2.d0)+phi4(i)**(2.d0))**(2.d0)) )) +&
        &(a_1(i)*((1.d0-a_1(i)**(2.d0))/(2.d0*xh(i))&
        +(xh(i)/4.d0)*(psi1(i+1)**(2.d0)+psi2(i+1)**(2.d0)+pi1(i+1)**(2.d0)+pi2(i+1)**(2.d0)&
        +2.d0*0.5d0*(mass*a_1(i))**(2.d0)*(phi1(i+1)**(2.d0)+phi2(i+1)**(2.d0))+psi3(i+1)**(2.d0)&
        +psi4(i+1)**(2.d0)+pi3(i+1)**(2.d0)+pi4(i+1)**(2.d0)&
        +2.d0*0.5d0*(mass*a_1(i))**(2.d0)*(phi3(i+1)**(2.d0)+phi4(i+1)**(2.d0))+&
        (a_1(i))**(2.d0)*lambda*((phi1(i+1)**(2.d0)+phi2(i+1)**(2.d0))**(2.d0)&
        +(phi3(i+1)**(2.d0)+phi4(i+1)**(2.d0))**(2.d0)) ) )) )

        a(i+1)=a_2(i)

    i=i+1
  end do


      
      !calculamos el otro coeficiente m trico desde la frontera hacia el origen del dominio computacional
      alpha(N)=1.d0/a(N)  
      !alpha(N)=1.d0  
      i=N

      do while(i>0)
        alpha_1(i)=alpha(i)-dr*alpha(i)*((1.d0-a(i)**(2.d0))/(-rh(i))+&
        (rh(i)/2.d0)*(psi1(i)**(2.d0)+psi2(i)**(2.d0)+pi1(i)**(2.d0)+pi2(i)**(2.d0)&
        -2.d0*0.5d0*(mass*a(i))**(2.d0)*(phi1(i)**(2.d0)&
        +phi2(i)**(2.d0))+psi3(i)**(2.d0)+psi4(i)**(2.d0)+pi3(i)**(2.d0)+pi4(i)**(2.d0)&
        -2.d0*0.5d0*(mass*a(i))**(2.d0)*(phi3(i)**(2.d0)+phi4(i)**(2.d0))-&
        (a(i))**(2.d0)*lambda*((phi1(i)**(2.d0)+phi2(i)**(2.d0))**(2.d0)&
        & +(phi3(i)**(2.d0)+phi4(i)**(2.d0))**(2.d0)) ))

        xh(i) = rh(i) - dr/1.d0

        alpha_2(i)=alpha(i)-dr*(alpha(i)*((1.d0-a(i)**(2.d0))/(-2.d0*rh(i))+&
        (rh(i)/4.d0)*(psi1(i)**(2.d0)+psi2(i)**(2.d0)+pi1(i)**(2.d0)+pi2(i)**(2.d0)&
        -2.d0*0.5d0*(mass*a(i))**(2.d0)*(phi1(i)**(2.d0)&
        +phi2(i)**(2.d0))+psi3(i)**(2.d0)+psi4(i)**(2.d0)+pi3(i)**(2.d0)+pi4(i)**(2.d0)&
        -2.d0*0.5d0*(mass*a(i))**(2.d0)*(phi3(i)**(2.d0)+phi4(i)**(2.d0))-&
        (a(i))**(2.d0)*lambda*((phi1(i)**(2.d0)+phi2(i)**(2.d0))**(2.d0)&
        & +(phi3(i)**(2.d0)+phi4(i)**(2.d0))**(2.d0)) )) + &
        alpha_1(i)*((1.d0-a(i-1)**(2.d0))/(-2.d0*xh(i))+&
        (xh(i)/4.d0)*(psi1(i-1)**(2.d0)+psi2(i-1)**(2.d0)+pi1(i-1)**(2.d0)+pi2(i-1)**(2.d0)&
        -2.d0*0.5d0*(mass*a(i-1))**(2.d0)*(phi1(i-1)**(2.d0)&
        +phi2(i-1)**(2.d0))+psi3(i-1)**(2.d0)+psi4(i-1)**(2.d0)+pi3(i-1)**(2.d0)+pi4(i-1)**(2.d0)&
        -2.d0*0.5d0*(mass*a(i-1))**(2.d0)*(phi3(i-1)**(2.d0)+phi4(i-1)**(2.d0))-&
        (a(i-1))**(2.d0)*lambda*((phi1(i-1)**(2.d0)+phi2(i-1)**(2.d0))**(2.d0)&
        & +(phi3(i-1)**(2.d0)+phi4(i-1)**(2.d0))**(2.d0)) )))

        alpha(i-1)=alpha_2(i)

    i=i-1
  end do
      alpha(-1)=alpha(0)

      l2norm_ligadura = 0.d0
      IF (j>=2) THEN
         i=0 !verificamos la ligadura de momento
         do while(i<N)
         !ligadura(i)=(a(i)-a_aux2(i))/((2.0)*dt)-0.5*rh(i)*alpha(i)*(psi1(i)*pi1(i)&
         !     +psi2(i)*pi2(i)+psi3(i)*pi3(i)+psi4(i)*pi4(i))
         !ligadura(i)=((3.d0*a(i)-4.d0*a_aux(i)+a_aux2(i))/(2.d0*dt))&
         !     -0.5*rh(i)*alpha(i)*(psi1(i)*pi1(i)+psi2(i)*pi2(i)+psi3(i)*pi3(i)+psi4(i)*pi4(i))
          ligadura(i)=((11.d0/3.d0*a(i)-6.d0*a_aux(i) &
          & + 3.d0*a_aux2(i)-2.d0/3.d0*a_aux3(i))/(2.d0*dt))-&
          0.5d0*rh(i)*alpha(i)*(psi1(i)*pi1(i)+psi2(i)*pi2(i)+psi3(i)*pi3(i)+psi4(i)*pi4(i))
          
          l2norm_ligadura = l2norm_ligadura + ligadura(i)**2
          i=i+1
          
         end do
      END IF


      !cálculo del número de partículas en cada instante de tiempo para cada campo. usando la regla del trapecio:
      particulas1=0.d0
      particulas2=0.d0
      i=0
      do while(i<=N)   !en vez de llegar hasta r_max=1000, llegaremos hasta el punto 1000/5=200
        particulas1=particulas1+PI**(2.d0)*(dr/2.d0)*((phi1(i-1)*pi2(i-1)-phi2(i-1)*pi1(i-1))*rh(i-1)**(4.d0)&
            +(phi1(i)*pi2(i)-phi2(i)*pi1(i))*rh(i)**(4.d0))
        particulas2=particulas2+PI**(2.d0)*(dr/2.d0)*((phi3(i-1)*pi4(i-1)-phi4(i-1)*pi3(i-1))*rh(i-1)**(4.d0)&
        +(phi3(i)*pi4(i)-phi4(i)*pi3(i))*rh(i)**(4.d0))

        i=i+1
      end do

      IF (MOD(j,1500).eq.0) THEN


      OPEN(UNIT=1, FILE='data/metrica.txt',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append')
        DO i = 0,N
          IF (rh(i)<200) THEN
            IF (MOD(i,5).eq.0) THEN
              write(1,*) a(i), alpha(i), rh(i), th(j), ligadura(i), sqrt(rh(i)*alpha(i)*(alpha(i+1)-alpha(i))/dr)  
              !HE AÑADIDO NUEVO AQUÍ LA VELOCIDAD DE ROTACIÓN
            END IF
          END IF
  
          IF (rh(i)>200) THEN
            IF (MOD(i,50).eq.0) THEN
            write(1,*) a(i), alpha(i), rh(i), th(j), ligadura(i), sqrt(rh(i)*alpha(i)*(alpha(i+1)-alpha(i))/dr) 
            END IF
          END IF
        ENDDO
      CLOSE(unit=1)

      OPEN(UNIT=2, FILE='data/campo.txt',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append')
            DO i = 0,N
              IF (rh(i)<200) THEN
                IF (MOD(i,5).eq.0) THEN
                  write(2,*) rh(i), th(j), phi1(i), phi2(i), phi3(i),phi4(i)
                END IF
              END IF
    
              IF (rh(i)>200) THEN
                IF (MOD(i,100).eq.0) THEN
                write(2,*) rh(i), th(j), phi1(i), phi2(i), phi3(i),phi4(i)
                END IF
              END IF
            ENDDO
        CLOSE(unit=2)

      END IF


      open(3, file='data/max_min.txt',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append')
      IF (MOD(j,20).eq.0) THEN
        write(3,*)  sqrt(phi1(1)*phi1(1) + phi2(1)*phi2(1)),sqrt(phi3(1)*phi3(1) + phi4(1)*phi4(1)), &
              minval(alpha), th(j),phi1(1),phi2(1), phi3(1),phi4(1)
      END IF
      close(3)
      
      IF (j>=4) THEN
        open(4, file='data/ligadura.txt',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append')
          IF (MOD(j,1).eq.0) THEN
          write(4,*) th(j), sqrt(l2norm_ligadura/N)
          END IF
        close(4)
          END IF
      
          open(5, file='data/particulas.txt',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append')
          IF (MOD(j,5).eq.0) THEN
            write(5,*)  th(j),particulas1,particulas2,particulas2/particulas1
          END IF
          close(5)
      
      j=j+1
      end do





end program autogravitante
