        module last_event_storage
            implicit none
            include 'coms.f'
            include 'comres.f'
            include 'comhis.f'
            include 'colltab.f'
            include 'freezeout.f'
            include 'options.f'
            private
            
            public :: set_evolution_essentials, 
     &          get_evolution_essentials,
     &          init_evolution_essentials

            type :: event_vars
            integer :: npart=0, nbar=0, nmes=0, apt=0, uid_cnt=0,
     &          ctag=0, dectag=0, NHardRes=0, NSoftRes=0,
     &          NDecRes=0, NBlColl=0,
     &          NElColl=0, origin(nmax)
            logical :: firstseed, fixedseed
            
            real*8 :: frr0(nmax), frrx(nmax), frry(nmax), frrz(nmax),
     &          frp0(nmax), frpx(nmax), frpy(nmax), frpz(nmax),
     &          dectime(nmax), tform(nmax), xtotfac(nmax), 
     &          ffermpx(nmax), ffermpy(nmax), ffermpz(nmax), 
     &          hisacttime(wwmax), acttime

            integer :: ncoll(nmax), lstcoll(nmax), uid(nmax)


            end type event_vars

            type(event_vars) :: instance_
            !$omp threadprivate(instance_)
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
                instance_%origin(1:npart)    = origin(1:npart)
                instance_%ncoll(1:npart)     = ncoll(1:npart)
                instance_%lstcoll(1:npart)   = lstcoll(1:npart)
                instance_%uid(1:npart)       = uid(1:npart)

                instance_%firstseed          = firstseed
                instance_%fixedseed          = fixedseed

                ! Real*8 Arrays
                instance_%frr0(1:npart)      = frr0(1:npart)
                instance_%frrx(1:npart)      = frrx(1:npart)
                instance_%frry(1:npart)      = frry(1:npart)
                instance_%frrz(1:npart)      = frrz(1:npart)
                instance_%frp0(1:npart)      = frp0(1:npart)
                instance_%frpx(1:npart)      = frpx(1:npart)
                instance_%frpy(1:npart)      = frpy(1:npart)
                instance_%frpz(1:npart)      = frpz(1:npart)
                instance_%dectime(1:npart)    = dectime(1:npart)
                instance_%tform(1:npart)     = tform(1:npart)
                instance_%xtotfac(1:npart)   = xtotfac(1:npart)
                instance_%ffermpx(1:npart)   = ffermpx(1:npart)
                instance_%ffermpy(1:npart)   = ffermpy(1:npart)
                instance_%ffermpz(1:npart)   = ffermpz(1:npart)
                instance_%hisacttime(1:npart) = hisacttime(1:npart)
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
                origin(1:npart)    = instance_%origin(1:npart)
                ncoll(1:npart)     = instance_%ncoll(1:npart)
                lstcoll(1:npart)   = instance_%lstcoll(1:npart)
                uid(1:npart)       = instance_%uid(1:npart)
                
                firstseed          = instance_%firstseed
                fixedseed          = instance_%fixedseed

                ! Real*8 Arrays
                frr0(1:npart)      = instance_%frr0(1:npart)
                frrx(1:npart)      = instance_%frrx(1:npart)
                frry(1:npart)      = instance_%frry(1:npart)
                frrz(1:npart)      = instance_%frrz(1:npart)
                frp0(1:npart)      = instance_%frp0(1:npart)
                frpx(1:npart)      = instance_%frpx(1:npart)
                frpy(1:npart)      = instance_%frpy(1:npart)
                frpz(1:npart)      = instance_%frpz(1:npart)
                dectime(1:npart)    = instance_%dectime(1:npart)
                tform(1:npart)     = instance_%tform(1:npart)
                xtotfac(1:npart)   = instance_%xtotfac(1:npart)
                ffermpx(1:npart)   = instance_%ffermpx(1:npart)
                ffermpy(1:npart)   = instance_%ffermpy(1:npart)
                ffermpz(1:npart)   = instance_%ffermpz(1:npart)
                hisacttime(1:npart) = instance_%hisacttime(1:npart)
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
                
                instance_%firstseed = .true.
                instance_%fixedseed = .true.
                
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