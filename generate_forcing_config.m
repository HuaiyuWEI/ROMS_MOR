function generate_forcing_config(interp_frc, filename)
    if nargin < 2
        filename = 'flux_frc.opt';
    end

    fid = fopen(filename, 'w');

    fprintf(fid, '      ! ****************************************************************\n');
    fprintf(fid, '      ! User inputs\n\n');

    % Interpolation flag
    fprintf(fid, '      integer     :: interp_frc = %d   ! interpolate forcing from coarser input grid (=1) or not (=0). factor 2 only for now\n\n', interp_frc);

    % Surface momentum stress
    fprintf(fid, '      type (ncforce) :: nc_sustr  = ncforce(vname=''sustr'', tname=''frc_time'')  ! sustr - surface u-momentum stress flux (input data in N/m^2)\n');
    fprintf(fid, '      type (ncforce) :: nc_svstr  = ncforce(vname=''svstr'', tname=''frc_time'')  ! svstr - surface v-momentum stress flux (input data in N/m^2)\n\n');

    % Surface heat and freshwater fluxes
    fprintf(fid, '      type (ncforce) :: nc_shflx  = ncforce(vname=''shflux'', tname=''frc_time'')  ! stflx(itemp) - surface heat flux (input data in W/m^2)\n');
    fprintf(fid, '      type (ncforce) :: nc_swflux = ncforce(vname=''swflux'', tname=''frc_time'')  ! stflx(isalt) - surface freshwater flux (input data in cm/day). Might want to use #if def SALINITY?\n');
    fprintf(fid, '      type (ncforce) :: nc_swrad  = ncforce(vname=''swrad'', tname=''frc_time'')  ! swrad - surface short-wave radiation flux (input data in W/m^2)\n\n');

    fprintf(fid, '      ! End of user inputs\n');
    fprintf(fid, '      ! ****************************************************************\n');

    fclose(fid);
    fprintf('Forcing config file written to "%s"\n', filename);
end
