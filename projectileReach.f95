MODULE VARIABLES
	REAL*8 :: A, B, GAMMA
	REAL*8 :: U0, Ux, Uy, ANG_VENTO, DENSITY
	REAL*8 :: V0, V0y, V0x, DIAM_MEDIO, ANG_PROJ, MASS, Y0, Cf
	REAL*8 :: ERRO
	REAL*8 :: PI = 3.1415926535897932384626433832795, g=9.8
END MODULE


PROGRAM CALCULADORA_DE_ALCANCE
USE VARIABLES
IMPLICIT NONE
LOGICAL :: EXISTS
CHARACTER (LEN=1) :: ASW


!O MÉTODO CONSISTE EM DETERMINAR O TEMPO DE QUEDA DE UM PROJÉTIL ATRAVÉ DAS EQUAÇÕES De MOVIMENTO 
!QUE LEVAM EM CONTA O ARRASTO DO VENTO. OU SEJA, DETERMINAR QUANDO A ALTURA É ZERO QUANDO T NÃO É ZERO.
!COM ESSE TEMPO DETERMINAMOS O ALCANCE DO PROJÉTIL

WRITE (*,*)""
WRITE (*,*)"---------------------------------------------------------------------------"
WRITE (*,*) "-CALCULADORA DE ALCANCE DE LANCAMENTO DE PROJETEIS COM RESISTENCIA DO AR E VENTO-"
WRITE (*,*)"---------------------------------------------------------------------------"
WRITE (*,*)""
20 CONTINUE
inquire(file="log.txt", exist=EXISTS)
  if (EXISTS) then
    open(12, file="log.txt", status="old", position="append", action="write")
  else
    open(12, file="log.txt", status="new", action="write")
  end if
WRITE (*,*) "ADOTE TODAS SUAS UNIDADES NO S.I."
CALL COLECT_DADOS_PROJETIL()
CALL COLECT_DADOS_AMBIENTE()
GAMMA = ((PI*DIAM_MEDIO*DIAM_MEDIO*DENSITY)/MASS)*Cf
A=(V0y - Uy + G/GAMMA)/GAMMA
B= Uy - G/GAMMA


WRITE (*,*)""
WRITE(*,*) "SEU COEFICIENTE DE ARRASTO E:", GAMMA
WRITE(12,*) "SEU COEFICIENTE DE ARRASTO E:", GAMMA
WRITE (*,*)""
WRITE (*,*) "PRIMEIRO SERA CALCULADO O TEMPO DE QUEDA E APARTIR DESTE VALOR, O ALCANCE"
WRITE (*,*)""
WRITE (*,*) "NOTA: VALORES NEGATIVOS DE ALCANCE -> QUEDA ATRAS DO PONTO DE LANCAMENTO"
WRITE (*,*)""
10 CONTINUE



WRITE (*,*) "QUAL METODO DESEJA USAR, NEWTOW-RAPHSON (n) OU FALSA POSICAO (p)?"
READ (*,*) ASW
WRITE (*,*) "QUAL O ERRO ABSOLUTO?"
READ (*,*) ERRO
IF (ASW=="n" .or. ASW == 'N') THEN
	CALL NEWTON()
ELSE IF (ASW=="p" .or. ASW == "P") THEN
  	CALL FALSA()
END IF


PRINT *, 'DESEJA TENTAR DE NOVO? (Y/N)'
READ *, ASW
IF(ASW == 'Y' .or. ASW == 'y')THEN
  GOTO 10
ELSE IF (ASW == 'N' .or. ASW == 'n') THEN
  CONTINUE
END IF


PRINT *, 'DESEJA CALCULAR OUTRO LANCAMENTO? (Y/N)'
READ *, ASW
IF(ASW == 'Y' .or. ASW == 'y')THEN
  GOTO 20
ELSE IF (ASW == 'N' .or. ASW == 'n') THEN
  STOP
END IF

END PROGRAM
!----------------
SUBROUTINE NEWTON()
	USE VARIABLES
	INTEGER I, m
    REAL*8 :: Y, DY, t, x, PHI
    write(12,*) "Tentativa com Newton"
	WRITE (*,*) "DIGITE O CHUTE INICIAL"
	READ (*,*) t
    Y= PHI(T)
	WRITE (*,*) Y
	WRITE (*,*) "DIGITE O NUMERO MAXIMO DE ITERACOES"
	READ (*,*) m
	I=0
	DO WHILE((I.ne.m) .and. (abs(Y).GE.ERRO))
	 	Y= PHI(t)
	 	DY= ((((GAMMA*t + 1)*EXP(-GAMMA*t))-1)/t*t)*A - Y0/t*t
	 	t=t-(Y/DY)
	 	I=I+1
	END DO
    X=((1-EXP(-GAMMA*t))/GAMMA)*(V0x-Ux) + Ux*t
    WRITE(*,*) "ITERACAO", i
    WRITE(*,*)"|PHI(T)|", ABS(Y) 
    WRITE(*,*)"TEMPO DE QUEDA", t 
    WRITE(*,*)"ALCANCE", X
    WRITE(12,*) "ITERACAO", i
    WRITE(12,*)"|PHI(T)|", ABS(Y) 
    WRITE(12,*)"TEMPO DE QUEDA", t 
    WRITE(12,*)"ALCANCE", X
    
END
SUBROUTINE FALSA()
	USE VARIABLES
    IMPLICIT NONE
	INTEGER I, m
    REAL*8 ::  Y, t, x, PHI
    REAL*8 :: F1,F2, e,f
    write(12,*) "Tentativa com falsa posicao"
	WRITE (*,*) "DIGITE UM INTERVALO [a,b] (SEPARANDO COM ENTER CADA VALOR) PARA SE PROCURAR DOIS PONTOS BONS INICIAIS"
    READ(*,*) e
    READ(*,*) f
    WRITE(*,*) PHI(e), PHI(f)
    CALL SEARCH_BEST_INTERVAL(E,F)
    WRITE (*,*) "NOVO INTERVALO"
    WRITE(*,*) '[',E,',',F,']'
    WRITE(*,*) "DIGITE O NUMERO MÁXIMO DE ITERAÇÕES"
    READ (*,*) m
    I=0
    F1= PHI(e)
    F2= PHI(f)
    T=(e*F2 - f*F1)/(F2-F1)
    Y=PHI(T)
    DO WHILE ((I.NE.m) .AND. (ABS(Y).GE.ERRO))
	  	F1= PHI(e)
    	F2= PHI(f)
	    T=(e*F2 - f*F1)/(F2-F1)
        IF(PHI(T)*PHI(e) .LT. 0) THEN
          f=T
        ELSE
          e=T
        END IF
       	Y=PHI(T)
        I=I+1
	ENDDO
    X=((1-EXP(-GAMMA*t))/GAMMA)*(V0x-Ux) + Ux*t
    WRITE(*,*) "ITERACAO", i
    WRITE(*,*)"|PHI(T)|", ABS(Y) 
    WRITE(*,*)"TEMPO DE QUEDA", t 
    WRITE(*,*)"ALCANCE", X
    WRITE(12,*) "ITERACAO", i
    WRITE(12,*)"|PHI(T)|", ABS(Y) 
    WRITE(12,*)"TEMPO DE QUEDA", t 
    WRITE(12,*)"ALCANCE", X
END
SUBROUTINE COLECT_DADOS_PROJETIL ()
	USE VARIABLES
    IMPLICIT NONE
	WRITE (*,*) "--------------------------------DADOS DO PROJETIL-----------------------------"
	write(12,*)"--------------------------------DADOS DO PROJETIL-----------------------------"
	WRITE (*,*) "QUAL O DIAMETRO MEDIO DAS SECCOES TRANSVERSAIS DE SEU PROJETIL?"
	READ (*,*) DIAM_MEDIO
    write(12,*) "diametro médio", DIAM_MEDIO
	WRITE (*,*) "QUAL A MASSA DESTE PROJETIL"
	READ (*,*) MASS
    write(12,*) "massa", MASS
	WRITE (*,*) "QUAL O MODULO DA VELOCIDADE INICIAL DO PROJETIL?"
	READ (*,*) V0
    write(12,*) "velocidade inicial", v0
	WRITE (*,*) "QUAL O ANGULO DE LANCAMENTO COM A HORIZONTAL? (GRAUS)(0<X<90)"
	READ (*,*) ANG_PROJ
    write(12,*) "angulo de lancamente", ANG_PROJ
    WRITE (*,*) "QUAL A ALTURA INICIAL?"
    READ (*,*) Y0
    write(12,*) "altura inicial", Y0
    WRITE (*,*) "QUAL O COEFICIENTE DE FORMA? (BALA : 0.288), (ESFERA: 0.7), (CUBO: 1.4), (GOLFINHO: 0.02)"
    READ (*,*) Cf
    write(12,*) "coeficiente de forma", Cf
	ANG_PROJ = ANG_PROJ*PI/180
	V0x = V0*COS(ANG_PROJ)
	V0Y = V0*SIN(ANG_PROJ)
END

SUBROUTINE COLECT_DADOS_AMBIENTE ()
	USE VARIABLES
    IMPLICIT NONE
	WRITE (*,*) "-------------------------------DADOS DO AMBIENTE------------------------------"
    write(12,*)"--------------------------------DADOS DO AMBIENTE-----------------------------"
	WRITE (*,*) "QUAL A DENSIDADE DO FLUIDO EM QUE OCORRE O LANCAMENTO?"
    WRITE(*,*) "(DENSIDADE DO AR A 1atm COM 25C = 1,28 Kg/m^3)"
	READ (*,*) DENSITY
    write(12,*) "densidade do ar", DENSITY
	WRITE (*,*) "QUAL O MODULO DA VELOCIDADE DE MOVIMENTO DO FLUIDO"
	READ (*,*) U0
    write(12,*) "velocidade do ar", U0
	WRITE (*,*) "QUAL O ANGULO COM A HORIZONTAL DA VELOCIDADE DO AR (GRAUS)"
	READ (*,*) ANG_VENTO
    write(12,*) "angulo do vento", ANG_VENTO
	ANG_VENTO=ANG_VENTO*PI/180
	Ux = U0*COS(ANG_VENTO)
	Uy = U0*SIN(ANG_VENTO)
END
SUBROUTINE SEARCH_BEST_INTERVAL (E,F)
    REAL*8 :: DELTA, E, F, a,b, PHI
   	INTEGER :: m
    DELTA = (F-E)/10
    a = E
    M=0
    DO WHILE(m.LT.10)
    	b = a + DELTA
        IF(PHI(a)*PHI(b) .LT. 0) THEN
			E =  a
            F =  b
            M=10
        END IF
        a = b
        M = M+1
	ENDDO
END 
REAL*8 FUNCTION PHI(Z)
	USE VARIABLES
	REAL*8 :: Z
	PHI =  ((1-EXP(-GAMMA*z))/z)*A + B + Y0/z
END