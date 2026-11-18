! "ana_vmix.h": vertical mixing for 2-D ridge runs (VMIX_PRESCRIBED_Huaiyu).
!
! Akv/Akt are the parent 3-D run's time- and along-ridge-mean values,
! read from the grid file by get_grid (grid.F, Akv_pres/Akt_pres). This
! routine is called ONCE at initialization (main.F, after get_init);
! with LMD_MIXING undefined nothing recomputes Akv/Akt afterwards, so the
! prescribed fields stay fixed for the whole run (and are re-applied on
! every restart, since get_grid and ana_vmix run again then).
! Akt is a single diffusivity shared by all tracers in this ROMS version.
#ifdef VMIX_PRESCRIBED_Huaiyu
      Akv = Akv_pres
      Akt = Akt_pres
      if (mynode==0) then
        write(*,'(6x,A,2ES11.3)') 'ana_vmix :: prescribed Akv min/max',
     &         minval(Akv_pres(1:nx,1:ny,:)), maxval(Akv_pres(1:nx,1:ny,:))
        write(*,'(6x,A,2ES11.3)') 'ana_vmix :: prescribed Akt min/max',
     &         minval(Akt_pres(1:nx,1:ny,:)), maxval(Akt_pres(1:nx,1:ny,:))
      endif
#endif
