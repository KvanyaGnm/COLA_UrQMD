        module last_event_storage
            implicit none
            include 'coms.f'
            include 'comres.f'
            include 'comhis.f'
            include 'colltab.f'
            include 'freezeout.f'
            private
            
            public :: set_evolution_essentials, 
     &          get_evolution_essentials,
     &          init_evolution_essentials

            type :: event_vars
            integer :: npart=0, nbar=0, nmes=0, apt=0, uid_cnt=0,
     &          ctag=0, dectag=0, NHardRes=0, NSoftRes=0,
     &          NDecRes=0, NBlColl=0,
     &          NElColl=0, origin(nmax)
            
            real*8 :: frr0(nmax), frrx(nmax), frry(nmax), frrz(nmax),
     &          frp0(nmax), frpx(nmax), frpy(nmax), frpz(nmax),
     &          dectime(nmax), tform(nmax), xtotfac(nmax), 
     &          ffermpx(nmax), ffermpy(nmax), ffermpz(nmax), 
     &          hisacttime(wwmax), acttime

            integer :: ncoll(nmax), lstcoll(nmax), uid(nmax)


            end type event_vars

            type(event_vars) :: instance_
            !$OMP THREADPRIVATE(instance_)

            contains

            subroutine set_evolution_essentials()
                instance_%npart     = npart
                instance_%nbar      = nbar
                instance_%nmes      = nmes
                instance_%apt       = apt
                instance_%uid_cnt   = uid_cnt
                instance_%ctag      = ctag
                instance_%dectag    = dectag
                instance_%NHardRes  = NHardRes
                instance_%NSoftRes  = NSoftRes
                instance_%NDecRes   = NDecRes
                instance_%NBlColl   = NBlColl
                instance_%NElColl   = NElColl
                instance_%acttime   = acttime
                
                ! Integer Arrays
                instance_%origin    = origin
                instance_%ncoll     = ncoll
                instance_%lstcoll   = lstcoll
                instance_%uid       = uid
                
                ! Real*8 Arrays
                instance_%frr0      = frr0
                instance_%frrx      = frrx
                instance_%frry      = frry
                instance_%frrz      = frrz
                instance_%frp0      = frp0
                instance_%frpx      = frpx
                instance_%frpy      = frpy
                instance_%frpz      = frpz
                instance_%dectime    = dectime
                instance_%tform     = tform
                instance_%xtotfac   = xtotfac
                instance_%ffermpx   = ffermpx
                instance_%ffermpy   = ffermpy
                instance_%ffermpz   = ffermpz
                instance_%hisacttime = hisacttime
            end subroutine set_evolution_essentials


            subroutine get_evolution_essentials()
                npart     = instance_%npart
                nbar      = instance_%nbar
                nmes      = instance_%nmes
                apt       = instance_%apt
                uid_cnt   = instance_%uid_cnt
                ctag      = instance_%ctag
                dectag    = instance_%dectag
                NHardRes  = instance_%NHardRes
                NSoftRes  = instance_%NSoftRes
                NDecRes   = instance_%NDecRes
                NBlColl   = instance_%NBlColl
                NElColl   = instance_%NElColl
                acttime   = instance_%acttime
                ! Integer Arrays
                origin    = instance_%origin
                ncoll     = instance_%ncoll
                lstcoll   = instance_%lstcoll
                uid       = instance_%uid
                
                ! Real*8 Arrays
                frr0      = instance_%frr0
                frrx      = instance_%frrx
                frry      = instance_%frry
                frrz      = instance_%frrz
                frp0      = instance_%frp0
                frpx      = instance_%frpx
                frpy      = instance_%frpy
                frpz      = instance_%frpz
                dectime    = instance_%dectime
                tform     = instance_%tform
                xtotfac   = instance_%xtotfac
                ffermpx   = instance_%ffermpx
                ffermpy   = instance_%ffermpy
                ffermpz   = instance_%ffermpz
                hisacttime = instance_%hisacttime
            end subroutine get_evolution_essentials


            subroutine init_evolution_essentials()
                instance_%npart     = 0
                instance_%nbar      = 0
                instance_%nmes      = 0
                instance_%apt       = 0
                instance_%uid_cnt   = 0
                instance_%ctag      = 0
                instance_%dectag    = 0
                instance_%NHardRes  = 0
                instance_%NSoftRes  = 0
                instance_%NDecRes   = 0
                instance_%NBlColl   = 0
                instance_%NElColl   = 0
                instance_%acttime   = 0.0
                ! Clear Integer Arrays
                instance_%origin    = 0
                instance_%ncoll     = 0
                instance_%lstcoll   = 0
                instance_%uid       = 0
                
                ! Clear Real*8 Arrays
                instance_%frr0      = 0.0D0
                instance_%frrx      = 0.0D0
                instance_%frry      = 0.0D0
                instance_%frrz      = 0.0D0
                instance_%frp0      = 0.0D0
                instance_%frpx      = 0.0D0
                instance_%frpy      = 0.0D0
                instance_%frpz      = 0.0D0
                instance_%dectime    = 0.0D0
                instance_%tform     = 0.0D0
                instance_%xtotfac   = 0.0D0
                instance_%ffermpx   = 0.0D0
                instance_%ffermpy   = 0.0D0
                instance_%ffermpz   = 0.0D0
                instance_%hisacttime = 0.0D0

                
                call get_evolution_essentials
            end subroutine init_evolution_essentials

        end module last_event_storage