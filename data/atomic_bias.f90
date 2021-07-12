! gfortran atomic_bias.f90 -o atomic_bias.x
! ./atomic_bias.x
program atomic_bias
      implicit none
      double precision, allocatable   :: R(:,:)
      character(len=2), allocatable   :: ele(:)
      integer                         :: iat, Nat, imol, Nmol, b, a, x
      integer, allocatable            :: method(:)
      integer                         :: nmethod, i, method_list(1:23)
      character(len=500)              :: f1, f2, t1, method_names(1:23), acc
      character(len=100)              :: comment, ac
      character(len=500), allocatable :: title(:)
      logical                         :: present1
      double precision, allocatable   :: HOF(:)
      double precision    :: Ga_c(1:23), Cl_c(1:23), As_c(1:23), H_c(1:23)
      double precision    :: Se_c(1:23), F_c(1:23), Be_c(1:23), O_c(1:23)
      double precision    :: K_c(1:23), Br_c(1:23), Mg_c(1:23), Li_c(1:23)
      double precision    :: B_c(1:23), Al_c(1:23), Ge_c(1:23), Ca_c(1:23)
      double precision    :: Na_c(1:23), P_c(1:23), Si_c(1:23), S_c(1:23)
      double precision    :: N_c1(1:23),  C_c1(1:23)
      integer    :: n_Ga, n_Cl, n_As, n_H
      integer    :: n_Se, n_F, n_Be, n_O
      integer    :: n_K, n_Br, n_Mg, n_Li
      integer    :: n_B, n_Al, n_Ge, n_Ca
      integer    :: n_Na, n_P, n_Si, n_S
      integer    :: n_N,  n_C

      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) "This program calculates the atomic correction needed "
      write(*,*) "for DFT when calculating a molecules heat of formation. "
      write(*,*) "The program has been parameterized for the following atoms"
      write(*,*) "     ____                            ____ " 
      write(*,*) "     |H |                            |He| "
      write(*,*) "     --------    ------------------------ "
      write(*,*) "     |Li| Be|    |B  |C  |N  |O  |F  |Ne| "
      write(*,*) "     |Na| Mg|    |Al |Si |P  |S  |Cl |Ar| "
      write(*,*) "     |K | Ca|    |Ga |Ge |As |Se |Br |Kr| "
      write(*,*) "     --------    ------------------------ "
      write(*,*) "                                          "
      write(*,*) "The program currently provides corrections for 23 popular"
      write(*,*) "functionals -- "
      write(*,*) " 1.  BLYP"
      write(*,*) " 2.  PW91"
      write(*,*) " 3.   PBE"
      write(*,*) " 4.  TPSS"
      write(*,*) " 5. B3LYP"
      write(*,*) " 6. B3LYP-D3"
      write(*,*) " 7. O3LYP"
      write(*,*) " 8. X3LYP"
      write(*,*) " 9. PBE0"
      write(*,*) "10. TPSS0"
      write(*,*) "11. TPSS0-D3"
      write(*,*) "12. M06-2X"
      write(*,*) "13. M06-2X-D3"
      write(*,*) "14. wB97X"
      write(*,*) "15. wB97X-D3"
      write(*,*) "16. wB97X-V"
      write(*,*) "17. wB97X-D3BJ"
      write(*,*) "18. CAM-B3LYP"
      write(*,*) "19. B2PLYP"
      write(*,*) "20. B2PLYP-D"
      write(*,*) "21. mPW2PLYP-D"
      write(*,*) "22. wB97M-D3BJ (coefficients set to 0 for Li, Na, K due to convergence issues.)"
      write(*,*) "23. wB97M-V    (coefficients set to 0 for Li, Na, K due to convergence issues.)"
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

      write(*,*) "Please enter the name of the xyz file (w/o .xyz)"
       read(*,*) f1
      write(f2,'(2a)') trim(f1), '.xyz'
      inquire(file=trim(f2), exist=present1)
      if ( present1 ) then
         a = 2
         write(*,*) "File found."
      else
         a = 4
         write(*,*) "File not found."
         write(*,*) "Please try again."
         do while (a == 4) 
          write(*,*) "Please enter the name of the xyz file (w/o .xyz)."
          read(*,*) f1
          write(f2,'(2a)') trim(f1), '.xyz'
          inquire(file=trim(f2), exist=present1)
          if ( present1 ) a = 2
          enddo
      endif

      write(*,*) "Please enter the number of molecules in xyz"
       read(*,*) Nmol
       allocate(title(1:Nmol))
      write(*,*) "Do the geometries have titles? [Y/N]"
       read(*,*) comment
       ac = comment(1:1)
       call To_upper(ac)
       if (trim(ac) .eq. 'Y') then
          write(*,*) "File contains titles."
       endif
       if (trim(ac) .eq. 'N') then
          write(*,*) "File does not contain titles."
       endif

       if (trim(ac) .ne. 'Y' .and. trim(ac) .ne. 'N') then
         write(*,*) "Please try again"
         b = 4
          do while (b == 4 )
           write(*,*) "Do the geometries have titles? [Y/N]"
           read(*,*) comment
           ac = comment(1:1)
           call To_upper(ac)
           if (trim(ac) .eq. 'Y') then
              b = 2
            write(*,*) "File contains titles."
           endif
           if (trim(ac) .eq. 'N') then
               b = 2
             write(*,*) "File does not contain titles."
           endif
          enddo
      endif

      method_list = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, &
                    14, 15, 16, 17, 18, 19, 20, 21, 22, 23 /)
      method_names = (/'      BLYP', '      PW91', '       PBE', '      TPSS', '     B3LYP', &
                       '  B3LYP-D3', '     O3LYP', '     X3LYP', '      PBE0', '     TPSS0', &
                       '  TPSS0-D3', '    M06-2X', ' M06-2X-D3', '     wB97X', '  wB97X-D3', &
                       '   wB97X-V', 'wB97X-D3BJ', ' CAM-B3LYP', '    B2PLYP', '   B2PLYPD', &
                       ' mPW2PLYPD', 'wB97M-D3BJ', '   wB97M-V'/)

      write(*,*) "Please enter number of methods."
       read(*,*) nmethod
       allocate(method(1:nmethod), HOF(1:nmethod))
       write(*,*) "Please choose the type of method(s) from the above list."
       write(*,*) "             (Press enter after each entry.)"
       x = 1
       do while (x .le. nmethod)
       read(*,*) method(x)
       if ( any(method_list==method(x)) ) then
       write(*,*) "Method chosen.", method(x), trim(method_names(method(x)))
       x = x + 1
       else 
       write(*,*) "Retry"
       endif       
       enddo
!       write(*,'(a,23i7)') "Methods chosen: ", (method(i),i=1,nmethod)
       write(*,'(a,23a)') "Method(s) chosen: ", (trim(method_names(method(i))),i=1,nmethod)
       write(*,*) '----------XXXXXXXXXX----------'
       write(*,*) '    Results start here        '
       write(acc,'(2a)') trim(f1), '_aocorrections_output.txt'
       open(unit=102,file=trim(acc))
       write(102,'(a,23a)') "Method(s) chosen: ", (trim(method_names(method(i))),i=1,nmethod)

!BLYP        PW91        PBE       TPSS      B3LYP   B3LYP-D3      O3LYP      X3LYP      PBE0   TPSS0      TPSS0-D3     M062X    M062XD3
 Al_c(1:13)=(/ -8.2968354d0, -7.6065042d0, -9.2784338d0, -3.1025452d0, -3.5061830d0, -1.8847274d0, -5.3762124d0,&
 -3.6208504d0, -4.6165503d0, -0.0149350d0,  1.2305921d0, 1.5383649d0,  1.5674670d0 /)
 As_c(1:13)=(/ 11.2751361d0,  8.4388393d0,  7.8439821d0,  8.2000612d0,  6.5839776d0,  7.1612106d0,  5.8758398d0,&
 6.0823045d0,  1.3042929d0,  1.7110734d0,  2.1320495d0, -4.8928758d0, -4.9033298d0 /)
 B_c(1:13)=(/ -0.0893373d0,  2.3393754d0,  1.9346766d0, -1.3042352d0,  0.8568218d0,  0.7252314d0,  3.6661200d0,&
 0.7198660d0,  1.3310846d0, -1.2152098d0, -1.3215324d0, 1.9644597d0,  1.9573132d0 /)
 Be_c(1:13) =(/  5.9378693d0,  2.2326284d0,  0.9985192d0,  1.2005685d0,  5.6244552d0,  5.0950281d0,  1.1625514d0,&
 5.2101446d0,  0.4780005d0,  0.7379510d0,  0.3406290d0, 2.3238554d0,  2.3177475d0 /)
 Br_c(1:13)=(/ -3.4568230d0,  2.0380216d0,  1.8232396d0, -0.1005123d0, -4.3389188d0, -2.2706558d0,  0.9254921d0,&
 -4.1752798d0, -0.9052391d0, -2.4087643d0, -0.8748097d0, -4.0414350d0, -4.0546447d0 /)
 C_c1(1:13) =(/ -0.4140257d0,  9.0467476d0,  9.1603347d0, -0.1861754d0, -1.3671356d0,  0.4670170d0,  9.7697832d0,&
 -1.3194874d0,  4.1129342d0, -3.3672929d0, -1.9787509d0, 0.5924039d0,  0.6125757d0 /)
 Ca_c(1:13)=(/  6.3206291d0,  7.5174590d0,  5.4879277d0, 11.6061688d0,  5.0751712d0,  5.2819474d0,  5.7141702d0, &
 3.8851431d0,  3.6413551d0,  8.3477748d0,  8.5169011d0, -2.8827750d0, -2.8777232d0 /)
 Cl_c(1:13) =(/ -4.7781441d0,  2.0353839d0,  1.9222771d0, -1.1367080d0, -4.7371379d0, -2.4501837d0,  2.4748680d0,&
 -4.4638473d0, -0.2062167d0, -2.6002269d0, -0.8915539d0, -0.1436146d0, -0.1464540d0 /)
 F_c(1:13) =(/  0.0855756d0,  5.2413345d0,  4.6649364d0,  0.2018914d0, -2.0493479d0, -1.6422433d0,  6.9449316d0,&
 -1.8757495d0, -0.8521549d0, -4.3171805d0, -4.0063944d0, 0.5132423d0,  0.5215380d0 /)
 Ga_c(1:13)=(/ -8.5863914d0, -7.1532073d0, -7.9171005d0, -0.7328821d0, -3.9462429d0, -5.2399008d0, -3.8560651d0,&
 -3.9285067d0, -2.6946944d0,  2.9458929d0,  2.0117375d0, 3.0783027d0,  3.0900519d0 /)
 Ge_c(1:13)=(/ -2.6576158d0, -3.0617383d0, -3.7005286d0, -0.3111748d0, -3.0352412d0, -2.3111274d0, -2.9811668d0,&
 -3.2984505d0, -5.0364384d0, -2.3328737d0, -1.7915931d0, -2.9878136d0, -2.9963114d0 /)
 H_c(1:13)=(/ -2.2064598d0, -2.2959955d0, -2.7236595d0,  0.5086940d0, -0.2737536d0,  0.0588147d0,  1.1530210d0,&
 -0.2672504d0, -1.4580862d0,  1.8130987d0,  2.0726574d0, -0.5958022d0, -0.5853890d0 /)
 K_c(1:13) =(/ -2.0986575d0, -3.6484287d0, -4.7677642d0, -0.1576185d0, -1.1627118d0, -1.7490026d0, -3.9071158d0,&
 -1.8075454d0, -3.7558104d0,  0.1207703d0, -0.3079150d0, -0.0041427d0, -0.0014909d0 /)
 Li_c(1:13)=(/ -1.5367881d0, -2.6330707d0, -3.3524597d0, -2.0004987d0, -1.2799079d0, -1.0869182d0, -2.1354921d0,&
 -1.1973241d0, -3.5005749d0, -1.9443259d0, -1.7996739d0, -0.2323250d0, -0.2339710d0 /)
 Mg_c(1:13)=(/ -4.0222988d0, -5.7430494d0, -7.3564954d0, -0.7077222d0, -1.3577857d0, -2.1617047d0, -4.8626870d0, &
 -1.9636945d0, -3.7506325d0,  1.3677628d0,  0.7837478d0, -4.5082120d0, -4.4991889d0 /)
 N_c1(1:13)=(/ 11.5652328d0, 14.2180747d0, 14.6460539d0,  5.4010535d0,  3.8477980d0,  4.1012072d0,  9.4492020d0, &
 3.0744849d0,  2.3573213d0, -5.0903492d0, -4.9053483d0, -0.4638286d0, -0.4839428d0 /)
 Na_c(1:13)=(/ -1.6778857d0, -2.3506502d0, -3.1166732d0, -1.1600472d0, -1.4895182d0, -1.4044551d0, -1.9495548d0,&
 -1.4383124d0, -3.2102972d0, -1.2592011d0, -1.1906954d0,  0.5343684d0,  0.5333051d0 /)
 O_c(1:13) =(/  5.7057832d0, 10.7285547d0, 10.5750020d0,  1.0937812d0, -0.8691888d0, -0.6544707d0,  9.5742091d0,&
 -1.3914695d0, -0.8867981d0, -8.2020562d0, -8.0405401d0, -1.7324852d0, -1.7236378d0 /)
 P_c(1:13)=(/ -0.9512452d0,  2.5041169d0,  2.2244498d0,  0.5406166d0, -1.9866538d0,  0.1441910d0,  3.3056011d0, &
 -2.2275363d0, -1.1723140d0, -2.7055925d0, -1.1146807d0, 3.0197386d0,  3.0175006d0 /)
 S_c(1:13) =(/ -0.2805680d0,  5.1068651d0,  5.0501572d0,  1.1845092d0, -2.4606933d0, -0.8962871d0,  4.9665006d0, &
 -2.6323670d0,  0.1674282d0, -2.8321679d0, -1.6735147d0, -0.5102934d0, -0.5310693d0 /)
 Se_c(1:13) =(/  4.9310452d0,  8.6946369d0,  8.5622939d0,  9.2046727d0,  0.4856682d0,  2.1304963d0,  6.0164802d0,&
 0.0977277d0,  1.5081167d0,  1.9826911d0,  3.1970830d0, -2.4053098d0, -2.4351146d0 /)
 Si_c(1:13) =(/ -2.7430260d0, -1.9178037d0, -2.7944780d0,  0.8683623d0, -1.7733875d0, -0.2642282d0, -0.7373465d0,&
 -1.9925633d0, -3.0200835d0, -0.4751246d0,  0.6561212d0, 0.9062441d0,  0.9197152d0 /)
!         wB97X   wB97X-D3   wB97X-V   wB97X-D3BJ   CAM-B3LYP    B2PLYP     B2PLYPD   mPW2PLYPD wB97M-D3BJ   wB97M-V
 Al_c(14:23) =(/  0.4588746d0,  0.0592195d0, -0.3612739d0,   6.4820125d0,  -1.6767276d0,  0.6400591d0,&
 1.7668827d0,  1.6146387d0,  3.7822670d0,  3.4940062d0 /)
 As_c(14:23) =(/  6.1773975d0,  6.9658020d0,  4.1592354d0,   7.0971634d0,   2.9720591d0,  4.8424327d0,&
 5.2540609d0,  3.8304508d0,  4.0051822d0,  4.7332082d0 /)
 B_c(14:23) =(/ -1.1995886d0, -1.1816508d0, -2.4179360d0,  -1.4943219d0,   0.4523249d0,  2.3783396d0,&
 2.3409939d0,  2.1535382d0,  2.2911286d0,  1.9479766d0 /)
 Be_c(14:23)  =(/  0.5226354d0, -0.3521238d0, -1.3003518d0,   0.2537958d0,   4.6700635d0,  5.6852966d0,&
 5.4376752d0,  5.4322213d0,  0.8813667d0, -0.4135797d0 /)
 Br_c(14:23) =(/ -2.9269317d0, -2.5876042d0, -1.3237229d0,  -0.0689727d0,  -4.5870474d0, -2.2596387d0,&
 -1.1474491d0, -2.2146955d0, -3.5001266d0, -2.1651432d0 /)
 C_c1(14:23)  =(/ -0.5243175d0, -0.4294871d0, -0.7361042d0,  -0.6910016d0,  -0.8149613d0,  0.9185418d0,&
 1.7879663d0,  0.9114513d0,  0.9003842d0,  1.7454060d0 /)
 Ca_c(14:23) =(/  5.7116832d0,  7.0500129d0,  5.9451875d0,  11.5192650d0,   5.2130832d0,  2.9989448d0,&
 3.3076781d0,  2.4480862d0,  2.1029818d0,  1.4844071d0 /)
 Cl_c(14:23)  =(/ -0.6813521d0, -0.7949261d0,  0.3462114d0,   1.0585488d0,  -3.8579575d0, -2.3910205d0,&
 -1.2341582d0, -2.0804184d0, -0.4232878d0,  0.7801183d0 /)
 F_c(14:23)  =(/  0.4864616d0, -0.1098209d0,  1.1532348d0,  -0.4522315d0,   0.8414764d0, -1.1573175d0,&
 -0.9960988d0, -1.0936408d0, -0.5908371d0,  0.7969339d0 /)
 Ga_c(14:23) =(/  0.3724798d0,  0.4032924d0, -0.0399287d0,   4.4094711d0,  -1.6658742d0,  1.4289439d0,&
 1.0442676d0,  1.8274882d0,  0.0518101d0,  0.8064740d0 /)
 Ge_c(14:23) =(/  0.9831766d0,  0.7554815d0, -0.9722038d0,   1.9714199d0,  -3.3741830d0,  2.2601697d0,&
 2.7588930d0,  1.6785527d0, -0.5482813d0,  0.1508677d0 /)
 H_c(14:23) =(/  0.2003524d0,  0.1443169d0,  0.4740354d0,   0.2365448d0,   0.5620299d0, -0.8444004d0, &
 -0.7164307d0, -0.1917812d0,  0.2020973d0,  0.1454705d0 /)
 K_c(14:23)  =(/  5.1130507d0,  2.9925251d0,  0.3327724d0,   2.7727987d0,  -0.2089223d0, -1.2334909d0,&
 -1.4128188d0, -1.4066850d0,  0.000000d0,  0.000000d0 /) ! Check
 Li_c(14:23) =(/ -1.7472696d0, -3.0076263d0, -2.4040561d0,  -1.3157931d0,  -1.1325479d0, -1.1457916d0,&
 -1.0144640d0, -0.7311489d0, 0.000000d0, 0.000000d0 /) ! Check
 Mg_c(14:23) =(/  1.4775857d0,  1.4992345d0, -0.8924393d0,   3.0898108d0,  -0.5485915d0,  0.4149921d0,&
 0.1262933d0,  0.1956737d0, -3.4257692d0, -4.7780246d0 /)
 N_c1(14:23) =(/  0.8404956d0,  0.7532098d0,  2.3942569d0,   1.8628232d0,   2.1754130d0,  3.0743470d0,&
 3.1731674d0,  1.4337607d0,  2.7341051d0,  3.3571695d0 /)
 Na_c(14:23) =(/  1.9561315d0, -0.2990952d0, -1.5593627d0,   0.2948736d0,  -1.4068612d0, -1.0119107d0,&
 -0.8857266d0, -0.5443164d0, 0.000000d0, 0.000000d0 /) ! Check
 O_c(14:23) =(/ -1.0649210d0, -0.8957580d0, -0.2543851d0,  -1.3702948d0,  -0.5224412d0,  0.1050999d0, &
 0.1878499d0, -1.3582764d0, -1.2746272d0, -0.2595423d0 /)
 P_c(14:23) =(/ -0.6606942d0, -0.4616755d0,  1.6679396d0,   5.2157321d0,  -3.3447829d0, -0.1796761d0, &
 0.9691539d0, -0.5702035d0,  4.2811561d0,  5.0305039d0 /)
 S_c(14:23) =(/ -1.3480934d0, -0.7057451d0,  0.1377852d0,   1.5025113d0,  -3.7710648d0, -0.9773748d0,&
 -0.1838131d0, -1.6378107d0,  1.0500527d0,  1.8184647d0 /)
 Se_c(14:23) =(/  0.4521921d0,  1.1212589d0,  1.7366681d0,   3.3540383d0,  -1.9252209d0,  5.0654518d0,& 
 5.9354307d0,  3.3858061d0, -0.7102763d0,  0.6846981d0 /)
 Si_c(14:23) =(/ -0.6888994d0, -0.5683087d0, -0.8371569d0,   2.9255819d0,  -2.0714559d0,  0.2632182d0,&
 1.1432797d0,  0.1733243d0,  3.2113726d0,  3.3790274d0 /)

       open(unit=101,file=trim(f2))
       do imol = 1, Nmol
          n_Ga   = 0  
          n_Cl   = 0
          n_As   = 0
          n_H    = 0
          n_Se   = 0
          n_F    = 0
          n_Be   = 0
          n_O    = 0
          n_K    = 0
          n_Br   = 0
          n_Mg   = 0
          n_Li   = 0
          n_B    = 0
          n_Al   = 0
          n_Ge   = 0
          n_Ca   = 0
          n_Na   = 0
          n_P    = 0
          n_Si   = 0
          n_S    = 0
          n_N    = 0
          n_C    = 0
       if (trim(ac) .eq. 'Y') then
       read(101,*) Nat
       read(101,*) title(imol)
       allocate(ele(1:Nat), R(1:Nat,1:3))
       do iat = 1, Nat
       read(101,*) ele(iat), R(iat,1:3)
       call To_upper(ele(iat))
       if (trim(ele(iat)) .eq. 'Ga') n_Ga = n_Ga + 1
       if (trim(ele(iat)) .eq. 'Cl') n_Cl = n_Cl + 1
       if (trim(ele(iat)) .eq. 'As') n_As = n_As + 1
       if (trim(ele(iat)) .eq.  'H')  n_H = n_H + 1
       if (trim(ele(iat)) .eq. 'Se') n_Se = n_Se + 1
       if (trim(ele(iat)) .eq.  'F')  n_F = n_F + 1
       if (trim(ele(iat)) .eq. 'Be') n_Be = n_Be + 1
       if (trim(ele(iat)) .eq.  'O')  n_O = n_O + 1
       if (trim(ele(iat)) .eq.  'K')  n_K = n_K + 1
       if (trim(ele(iat)) .eq. 'Br') n_Br = n_Br + 1
       if (trim(ele(iat)) .eq. 'Mg') n_Mg = n_Mg + 1
       if (trim(ele(iat)) .eq. 'Li') n_Li = n_Li + 1
       if (trim(ele(iat)) .eq.  'B')  n_B = n_B + 1
       if (trim(ele(iat)) .eq. 'Al') n_Al = n_Al + 1
       if (trim(ele(iat)) .eq. 'Ge') n_Ge = n_Ge + 1
       if (trim(ele(iat)) .eq. 'Ca') n_Ca = n_Ca + 1
       if (trim(ele(iat)) .eq. 'Na') n_Na = n_Na + 1
       if (trim(ele(iat)) .eq.  'P')  n_P = n_P + 1
       if (trim(ele(iat)) .eq. 'Si') n_Si = n_Si + 1
       if (trim(ele(iat)) .eq.  'S')  n_S = n_S + 1
       if (trim(ele(iat)) .eq.  'N')  n_N = n_N + 1
       if (trim(ele(iat)) .eq.  'C')  n_C = n_C + 1
       enddo
       HOF = 0.0d0
       do i = 1, nmethod
       HOF(i) = dfloat(n_Ga)*Ga_c(method(i)) + dfloat(n_Cl)*Cl_c(method(i)) + &
                dfloat(n_As)*As_c(method(i)) + dfloat(n_H)*H_c(method(i)) + &
                dfloat(n_Se)*Se_c(method(i)) + dfloat(n_F)*F_c(method(i)) + &
                dfloat(n_Be)*Be_c(method(i)) + dfloat(n_O)*O_c(method(i)) + &
                  dfloat(n_K)*K_c(method(i)) + dfloat(n_Br)*Br_c(method(i)) + &
                dfloat(n_Mg)*Mg_c(method(i)) + dfloat(n_Li)*Li_c(method(i)) + &
                  dfloat(n_B)*B_c(method(i)) +  dfloat(n_Al)*Al_c(method(i)) + &
                dfloat(n_Ge)*Ge_c(method(i)) + dfloat(n_Ca)*Ca_c(method(i)) + &
                dfloat(n_Na)*Na_c(method(i)) + dfloat(n_P)*P_c(method(i)) + &
                dfloat(n_Si)*Si_c(method(i)) + dfloat(n_S)*S_c(method(i)) + &
                 dfloat(n_N)*N_c1(method(i)) +  dfloat(n_C)*C_c1(method(i))  
       enddo
       write(*,'(a, 23f17.2)') trim(title(imol)), (HOF(i), i=1,nmethod)
       write(102,'(a, 23f17.2)') trim(title(imol)), (HOF(i), i=1,nmethod)
       deallocate(ele,R)
       endif

       if (trim(ac) .eq. 'N') then
          read(101,*) Nat
          read(101,*)
       allocate(ele(1:Nat), R(1:Nat,1:3))
       do iat = 1, Nat
       read(101,*) ele(iat), R(iat,1:3)
       if (trim(ele(iat)) .eq. 'Ga') n_Ga = n_Ga + 1
       if (trim(ele(iat)) .eq. 'Cl') n_Cl = n_Cl + 1
       if (trim(ele(iat)) .eq. 'As') n_As = n_As + 1
       if (trim(ele(iat)) .eq.  'H')  n_H = n_H + 1
       if (trim(ele(iat)) .eq. 'Se') n_Se = n_Se + 1
       if (trim(ele(iat)) .eq.  'F')  n_F = n_F + 1
       if (trim(ele(iat)) .eq. 'Be') n_Be = n_Be + 1
       if (trim(ele(iat)) .eq.  'O')  n_O = n_O + 1
       if (trim(ele(iat)) .eq.  'K')  n_K = n_K + 1
       if (trim(ele(iat)) .eq. 'Br') n_Br = n_Br + 1
       if (trim(ele(iat)) .eq. 'Mg') n_Mg = n_Mg + 1
       if (trim(ele(iat)) .eq. 'Li') n_Li = n_Li + 1
       if (trim(ele(iat)) .eq.  'B')  n_B = n_B + 1
       if (trim(ele(iat)) .eq. 'Al') n_Al = n_Al + 1
       if (trim(ele(iat)) .eq. 'Ge') n_Ge = n_Ge + 1
       if (trim(ele(iat)) .eq. 'Ca') n_Ca = n_Ca + 1
       if (trim(ele(iat)) .eq. 'Na') n_Na = n_Na + 1
       if (trim(ele(iat)) .eq.  'P')  n_P = n_P + 1
       if (trim(ele(iat)) .eq. 'Si') n_Si = n_Si + 1
       if (trim(ele(iat)) .eq.  'S')  n_S = n_S + 1
       if (trim(ele(iat)) .eq.  'N')  n_N = n_N + 1
       if (trim(ele(iat)) .eq.  'C')  n_C = n_C + 1
       enddo
       HOF = 0.0d0
       do i = 1, nmethod
       HOF(i) = dfloat(n_Ga)*Ga_c(method(i)) + dfloat(n_Cl)*Cl_c(method(i)) + &
                dfloat(n_As)*As_c(method(i)) + dfloat(n_H)*H_c(method(i)) + &
                dfloat(n_Se)*Se_c(method(i)) + dfloat(n_F)*F_c(method(i)) + &
                dfloat(n_Be)*Be_c(method(i)) + dfloat(n_O)*O_c(method(i)) + &
                  dfloat(n_K)*K_c(method(i)) + dfloat(n_Br)*Br_c(method(i)) + &
                dfloat(n_Mg)*Mg_c(method(i)) + dfloat(n_Li)*Li_c(method(i)) + &
                  dfloat(n_B)*B_c(method(i)) +  dfloat(n_Al)*Al_c(method(i)) + &
                dfloat(n_Ge)*Ge_c(method(i)) + dfloat(n_Ca)*Ca_c(method(i)) + &
                dfloat(n_Na)*Na_c(method(i)) + dfloat(n_P)*P_c(method(i)) + &
                dfloat(n_Si)*Si_c(method(i)) + dfloat(n_S)*S_c(method(i)) + &
                 dfloat(n_N)*N_c1(method(i)) +  dfloat(n_C)*C_c1(method(i))
       enddo
       write(t1,'(a,i5.5)') 'Molecule_', imol
       write(*,'(a, 23f12.2)') trim(t1), (HOF(i), i=1,nmethod)
       write(102,'(a, 23f12.2)') trim(t1), (HOF(i), i=1,nmethod)
       deallocate(ele,R)
       endif
      
       enddo
       write(*,*) 'Results printed to ', trim(acc)
       deallocate(method,title,HOF)
       close(102)
      end program atomic_bias

subroutine To_upper(str) 
   ! https://rosettacode.org/wiki/String_case#Fortran
   ! Original author -- Clive Page
   character(*), intent(inout) :: str
   integer :: i

   do i = 1, len(str)
     select case(str(1:1))
     case("a":"z")
        str(i:i) = achar(iachar(str(i:i))-32)
     end select
   end do
end subroutine To_upper
