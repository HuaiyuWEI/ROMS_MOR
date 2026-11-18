classdef topo
    % Toolbox for analyzing ROMS simulations
    % organized by Huaiyu
    % V4 July 2026
    %
    % Method groups:
    %   vertical grid : zlevs3
    %   staggering    : rho2u, rho2v, rho2uvp, u2rho, v2rho,
    %                   rho2w, rho2w2d, rho2w4d, w2rho, w2rho2d
    %   velocities    : rebuildWvel  (w from omega and the u slope terms)
    %   masks         : uvp_masks
    %   derivatives   : grad_2d, diffxi, diffxi3d, diffeta, diffeta3d
    %   misc          : findoutput, do_bottom_avg2D
    methods(Static)

        function [z,Cs] = zlevs3(h,zeta,theta_s,theta_b,hc,N,type,scoord)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  function z = zlevs3(h,zeta,theta_s,theta_b,hc,N,type,scoord)
            %
            %  Compute ROMS rho- or w-point depths using the zlevs3
            %  stretching convention.
            %
            %  type   : 'r' for rho points, 'w' for w points
            %  scoord : 'old1994', 'new2006', or 'new2008'
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if nargin < 8
                error('Not enough input arguments')
            end

            type = char(type);
            scoord = char(scoord);

            [M, L] = size(h);

            % Set S-curves in domain [-1 < sc < 0] at vertical W- and
            % RHO-points.
            if strcmp(type, 'w')
                sc = ((0:N) - N) / N;
                N = N + 1;
            elseif strcmp(type, 'r')
                sc = ((1:N) - N - 0.5) / N;
            else
                error('type must be ''r'' or ''w''.')
            end

            if strcmp(scoord, 'new2008')
                % New s-coordinate allowing smooth bottom refinement.
                % This follows the zlevs3 CSF convention, including the
                % theta_s = 0 fallback csrf = -sc.^2.
                if theta_s > 0.0
                    csrf = (1.0 - cosh(theta_s * sc)) ./ ...
                           (cosh(theta_s) - 1.0);
                else
                    csrf = -(sc.^2);
                end
                sc1 = csrf + 1.0;

                if theta_b > 0.0
                    Cs = (exp(theta_b * sc1) - 1.0) ./ ...
                         (exp(theta_b) - 1.0) - 1.0;
                else
                    Cs = csrf;
                end
            else
                % For old1994 and new2006 s-coordinate.
                if theta_s == 0
                    error(['theta_s = 0 is only supported here for ', ...
                           'scoord = ''new2008''.'])
                end
                cff1 = 1.0 ./ sinh(theta_s);
                cff2 = 0.5 ./ tanh(0.5 * theta_s);
                Cs = (1.0 - theta_b) * cff1 * sinh(theta_s * sc) ...
                    + theta_b * (cff2 * tanh(theta_s * (sc + 0.5)) - 0.5);
            end

            z = zeros(N, M, L);
            if strcmp(scoord, 'old1994')
                disp('--- using old s-coord')
                hinv = 1.0 ./ h;
                cff = hc * (sc - Cs);
                cff1 = Cs;
                for k = 1:N
                    z0 = cff(k) + cff1(k) * h;
                    z(k,:,:) = z0 + zeta .* (1.0 + z0 .* hinv);
                end
            else
                if strcmp(scoord, 'new2006')
                    disp('--- using new s-coord (2006)')
                end
                hinv = 1.0 ./ (h + hc);
                cff = hc * sc;
                cff1 = Cs;
                for k = 1:N
                    z(k,:,:) = zeta + (zeta + h) .* ...
                        (cff(k) + cff1(k) * h) .* hinv;
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Staggered-grid interpolation in the horizontal.
        % All of these accept 2D (x,y), 3D (x,y,z) or 4D (x,y,z,t) input;
        % the averaging acts on the x (dim 1) or y (dim 2) axis and any
        % trailing dimensions pass through unchanged.

        function var_u = rho2u(var_rho)
            var_u = 0.5 * (var_rho(2:end,:,:,:) + var_rho(1:end-1,:,:,:));
        end

        function var_v = rho2v(var_rho)
            var_v = 0.5 * (var_rho(:,2:end,:,:) + var_rho(:,1:end-1,:,:));
        end

        function [ufield,vfield,pfield] = rho2uvp(rfield)
            vfield = 0.5*(rfield(:,1:end-1,:,:) + rfield(:,2:end,:,:));
            ufield = 0.5*(rfield(1:end-1,:,:,:) + rfield(2:end,:,:,:));
            pfield = 0.5*(ufield(:,1:end-1,:,:) + ufield(:,2:end,:,:));
        end

        function var_rho = u2rho(var_u)
            % Interior rho points are u-point averages; the first/last
            % x-rows are copied from their neighbours.
            sz = size(var_u);
            sz(1) = sz(1) + 1;
            var_rho = zeros(sz);
            var_rho(2:end-1,:,:,:) = 0.5 * (var_u(1:end-1,:,:,:) + var_u(2:end,:,:,:));
            var_rho(1,:,:,:) = var_rho(2,:,:,:);
            var_rho(end,:,:,:) = var_rho(end-1,:,:,:);
        end

        function var_rho = v2rho(var_v)
            % Interior rho points are v-point averages; the first/last
            % y-columns are copied from their neighbours.
            sz = size(var_v);
            sz(2) = sz(2) + 1;
            var_rho = zeros(sz);
            var_rho(:,2:end-1,:,:) = 0.5 * (var_v(:,1:end-1,:,:) + var_v(:,2:end,:,:));
            var_rho(:,1,:,:) = var_rho(:,2,:,:);
            var_rho(:,end,:,:) = var_rho(:,end-1,:,:);
        end

        % Dimension-suffixed aliases of the generic implementations above,
        % which handle any input dimensionality.
        function out = rho2u_2d(in), out = topo.rho2u(in); end
        function out = rho2v_2d(in), out = topo.rho2v(in); end
        function out = u2rho_2d(in), out = topo.u2rho(in); end
        function out = v2rho_2d(in), out = topo.v2rho(in); end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Vertical rho <-> w interpolation.

        function var_w = rho2w(var_r, z_r, z_w)
            % Interpolate a rho-level field onto w levels (linear in z,
            % linear extrapolation to the bottom/surface w points).
            % var_r may be 3D (x,y,z) or 4D (x,y,z,t); z_r/z_w are 3D and
            % broadcast over the trailing time dimension.
            var_w = zeros([size(z_w,1), size(z_w,2), size(z_w,3), size(var_r,4)]);

            % bottom layer computation
            var_w(:,:,1,:) = var_r(:,:,1,:) + (var_r(:,:,2,:) - var_r(:,:,1,:)) .* ...
                (z_w(:,:,1) - z_r(:,:,1)) ./ (z_r(:,:,2) - z_r(:,:,1));

            % surface layer computation
            var_w(:,:,end,:) = var_r(:,:,end,:) + (var_r(:,:,end,:) - var_r(:,:,end-1,:)) .* ...
                (z_w(:,:,end) - z_r(:,:,end)) ./ (z_r(:,:,end) - z_r(:,:,end-1));

            % Middle layers computation
            var_w(:,:,2:end-1,:) = var_r(:,:,1:end-1,:) .* (z_r(:,:,2:end) - z_w(:,:,2:end-1)) + ...
                var_r(:,:,2:end,:)   .* (z_w(:,:,2:end-1) - z_r(:,:,1:end-1));
            var_w(:,:,2:end-1,:) = var_w(:,:,2:end-1,:) ./ (z_r(:,:,2:end) - z_r(:,:,1:end-1));
        end

        function var_w = rho2w4d(var_r, z_r, z_w)
            % Alias of rho2w, which handles 4D input directly.
            var_w = topo.rho2w(var_r, z_r, z_w);
        end

        function var_w = rho2w2d(var_r, z_r, z_w)
            % 2D (x,z) version of rho2w.
            var_w = zeros(size(z_w));

            % bottom layer computation
            var_w(:,1) = var_r(:,1) + (var_r(:,2) - var_r(:,1)) .* ...
                (z_w(:,1) - z_r(:,1)) ./ (z_r(:,2) - z_r(:,1));

            % surface layer computation
            var_w(:,end) = var_r(:,end) + (var_r(:,end) - var_r(:,end-1)) .* ...
                (z_w(:,end) - z_r(:,end)) ./ (z_r(:,end) - z_r(:,end-1));

            % Middle layers computation
            var_w(:,2:end-1) = var_r(:,1:end-1) .* (z_r(:,2:end) - z_w(:,2:end-1)) + ...
                var_r(:,2:end)   .* (z_w(:,2:end-1) - z_r(:,1:end-1));
            var_w(:,2:end-1) = var_w(:,2:end-1) ./ (z_r(:,2:end) - z_r(:,1:end-1));
        end

        function var_rho = w2rho(var_w)
            % Interpolate a w-level field onto rho levels (ROMS 4th-order
            % weights). var_w may be 3D (x,y,z) or 4D (x,y,z,t).
            N = size(var_w, 3);
            sz = size(var_w);
            sz(3) = N - 1;
            var_rho = zeros(sz);

            % Compute for the middle layers
            for iz = 2:N-2
                var_rho(:,:,iz,:) = 0.5625 * (var_w(:,:,iz+1,:) + var_w(:,:,iz,:)) ...
                    - 0.0625 * (var_w(:,:,iz+2,:) + var_w(:,:,iz-1,:));
            end

            % Compute for the first and last layers separately
            var_rho(:,:,1,:) = -0.125 * var_w(:,:,3,:) + 0.75 * var_w(:,:,2,:) + 0.375 * var_w(:,:,1,:);
            var_rho(:,:,end,:) = -0.125 * var_w(:,:,N-2,:) + 0.75 * var_w(:,:,N-1,:) + 0.375 * var_w(:,:,N,:);
        end

        function var_rho = w2rho2d(var_w)
            % 2D (x,z) version of w2rho.
            [M,  N] = size(var_w);

            var_rho = zeros(M,  N-1);

            % Compute for the middle layers
            for iz = 2:N-2
                var_rho(:,iz) = 0.5625 * (var_w(:,iz+1) + var_w(:,iz)) - 0.0625 * (var_w(:,iz+2) + var_w(:,iz-1));
            end

            % Compute for the first and last layers separately
            var_rho(:,1) = -0.125 * var_w(:,3) + 0.75 * var_w(:,2) + 0.375 * var_w(:,1);
            var_rho(:,end) = -0.125 * var_w(:,N-2) + 0.75 * var_w(:,N-1) + 0.375 * var_w(:,N);
        end

        function w_rho = rebuildWvel(omega_w, uvel_u, z2D_r, dx)
            % Rebuild the vertical velocity at rho points on a y-mean
            % (x,z) section, following ROMS src/wvlcty.F:
            %
            %   w = <omega interpolated from the w faces to rho points>
            %       + 0.5*u(i-1/2)*(z_r(i)  -z_r(i-1))/dx
            %       + 0.5*u(i+1/2)*(z_r(i+1)-z_r(i)  )/dx
            %
            % (uniform grid, so wvlcty's pm(i)+pm(i+-1) reduces to 2/dx).
            % The two one-sided u contributions are summed, as in wvlcty.
            %
            % When it is needed: in ROMS builds whose basic_output.F passes
            % the whole `work` array (allocated 0:nz) to wvlcty in the
            % averaging path, the 'w' of the _avg files is shifted by one
            % level: the file stores w(k+1) at level k and never stores
            % w(1). omega, u, v and temp are written by direct assignment
            % and are not affected, so rebuilding from them recovers every
            % level. The development branch at commit a0e2896 passes
            % work(:,:,1:nz) in both the history and the averaging path.
            %
            % NOTE: the face-to-rho step is a plain midpoint average, NOT
            % topo.w2rho2d -- wvlcty uses 2nd-order midpoints, and the
            % 4th-order weights in w2rho2d would not reproduce the model.
            %
            %   omega_w : (Nx, Nlayer+1) w-face omega in m/s, with
            %             omega_w(:,1) = 0 at the seabed, matching wvlcty's
            %             implicit zero bottom interface
            %   uvel_u  : (Nx-1, Nlayer) u at u points
            %   z2D_r   : (Nx, Nlayer) rho-point depths
            %
            % The v*dz/deta terms are dropped: on a y-mean section of a
            % y-periodic domain with h uniform in y, z_r carries no
            % along-ridge structure beyond the free surface itself. The two
            % boundary columns are left at the omega-only value.
            [Nx, Nlayer] = size(z2D_r);

            % w faces -> rho points; the bottom cell is effectively halved
            % because the seabed interface is zero
            w_rho = 0.5*(omega_w(:, 2:end) + omega_w(:, 1:end-1));

            % u * dz/dxi, as the sum of the two one-sided face contributions
            slope = zeros(Nx, Nlayer);
            slope(2:end-1, :) = ...
                  0.5*uvel_u(1:end-1, :).*(z2D_r(2:end-1, :) - z2D_r(1:end-2, :))./dx ...
                + 0.5*uvel_u(2:end,   :).*(z2D_r(3:end,   :) - z2D_r(2:end-1, :))./dx;

            w_rho = w_rho + slope;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This function computes the Land/Sea masks on U-, V-, and PSI-points
        % from the mask on RHO-points.
        function [umask,vmask,pmask]=uvp_masks(rmask)
            [Lp,Mp]=size(rmask);
            L=Lp-1;
            M=Mp-1;

            %  Land/Sea mask on U-points.

            umask(1:L,1:Mp)=rmask(2:Lp,1:Mp).*rmask(1:L,1:Mp);

            %  Land/Sea mask on V-points.

            vmask(1:Lp,1:M)=rmask(1:Lp,2:Mp).*rmask(1:Lp,1:M);

            %  Land/Sea mask on PSI-points. The old computation of PSI-mask is replaced
            %  with the newer scheme.  Notice the shift because Matlab doesnot support
            %  zero index arrays.

            pmask=nan([L M]);

            for jr=2:Mp
                jp=jr-1;
                for ir=2:Lp
                    ip=ir-1;
                    if ((rmask(ir-1,jr  ) > 0.5) &&                     ...
                            (rmask(ir  ,jr  ) > 0.5) &&                 ...
                            (rmask(ir-1,jr-1) > 0.5) &&                 ...
                            (rmask(ir  ,jr-1) > 0.5))
                        pmask(ip,jp)=1;
                    elseif ((rmask(ir-1,jr  ) < 0.5) &&                 ...
                            (rmask(ir  ,jr  ) > 0.5) &&                 ...
                            (rmask(ir-1,jr-1) > 0.5) &&                 ...
                            (rmask(ir  ,jr-1) > 0.5))
                        pmask(ip,jp)=1;
                    elseif ((rmask(ir-1,jr  ) > 0.5) &&                 ...
                            (rmask(ir  ,jr  ) < 0.5) &&                 ...
                            (rmask(ir-1,jr-1) > 0.5) &&                 ...
                            (rmask(ir  ,jr-1) > 0.5))
                        pmask(ip,jp)=1;
                    elseif ((rmask(ir-1,jr  ) > 0.5) &&                 ...
                            (rmask(ir  ,jr  ) > 0.5) &&                 ...
                            (rmask(ir-1,jr-1) < 0.5) &&                 ...
                            (rmask(ir  ,jr-1) > 0.5))
                        pmask(ip,jp)=1;
                    elseif ((rmask(ir-1,jr  ) > 0.5) &&                 ...
                            (rmask(ir  ,jr  ) > 0.5) &&                 ...
                            (rmask(ir-1,jr-1) > 0.5) &&                 ...
                            (rmask(ir  ,jr-1) < 0.5))
                        pmask(ip,jp)=1;
                    elseif ((rmask(ir-1,jr  ) > 0.5) &&                 ...
                            (rmask(ir  ,jr  ) < 0.5) &&                 ...
                            (rmask(ir-1,jr-1) > 0.5) &&                 ...
                            (rmask(ir  ,jr-1) < 0.5))
                        pmask(ip,jp)=2;
                    elseif ((rmask(ir-1,jr  ) < 0.5) &&                 ...
                            (rmask(ir  ,jr  ) > 0.5) &&                 ...
                            (rmask(ir-1,jr-1) < 0.5) &&                 ...
                            (rmask(ir  ,jr-1) > 0.5))
                        pmask(ip,jp)=2;
                    elseif ((rmask(ir-1,jr  ) > 0.5) &&                 ...
                            (rmask(ir  ,jr  ) > 0.5) &&                 ...
                            (rmask(ir-1,jr-1) < 0.5) &&                 ...
                            (rmask(ir  ,jr-1) < 0.5))
                        pmask(ip,jp)=2;
                    elseif ((rmask(ir-1,jr  ) < 0.5) &&                 ...
                            (rmask(ir  ,jr  ) < 0.5) &&                 ...
                            (rmask(ir-1,jr-1) > 0.5) &&                 ...
                            (rmask(ir  ,jr-1) > 0.5))
                        pmask(ip,jp)=2;
                    else
                        pmask(ip,jp)=0;
                    end
                end
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function for calculating horizontal gradients of a 2D variable.
        % NOTE: boundary rows/columns 1 and N are left at zero, and the
        % rows/columns 2 and N-1 use first-order one-sided differences.
        function  [var_dx,var_dy] = grad_2d(var,dx,dy)
            [Nx, Ny] = size(var);
            var_dx = zeros(Nx, Ny);
            var_dy = zeros(Nx, Ny);
            % Compute the gradient in the x-direction (var_dx)
            var_dx(3:Nx-2, :) = (var(4:Nx-1, :) - var(2:Nx-3, :)) / (2 * dx);
            var_dx(2, :) = (var(3, :) - var(2, :)) / dx;
            var_dx(Nx-1, :) = (var(Nx-1, :) - var(Nx-2, :)) / dx;
            % Compute the gradient in the y-direction (var_dy)
            var_dy(:, 3:Ny-2) = (var(:, 4:Ny-1) - var(:, 2:Ny-3)) / (2 * dy);
            var_dy(:, 2) = (var(:, 3) - var(:, 2)) / dy;
            var_dy(:, Ny-1) = (var(:, Ny-1) - var(:, Ny-2)) / dy;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function for finding the output file names
        function filelist = findoutput(exppath,outputName)
            % Return the sorted <outputName>.<digits>.nc files. Only names
            % with a single all-digit token between the prefix and .nc are
            % kept, so joined full-domain output is used and per-node
            % partition files (e.g. neptune_his.<stamp>.000.nc) that
            % ncjoin has not consumed yet are excluded. To collect the
            % rank files of one stamp instead, include the stamp in
            % outputName (e.g. 'neptune_par.<stamp>').
            filelist_all = dir(fullfile(exppath, [outputName '.*.nc']));
            if isempty(filelist_all)
                filelist = cell(0,1);
                return;
            end

            filenames = {filelist_all.name};
            keepPattern = ['^' regexptranslate('escape', outputName) '\.\d+\.nc$'];
            keep = ~cellfun('isempty', regexp(filenames, keepPattern, 'once'));

            % Build a sorted cell array of full file paths.
            filenames = sort(filenames(keep));
            filelist = fullfile(exppath, filenames(:));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%% do_bottom_avg2D
        %%%
        %%% Thickness-weighted average of a 2D (x,z) field within a given
        %%% height above the bottom, on a refined vertical grid.
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function var_bottom_avg = do_bottom_avg2D(var, depthforavg, z2D, Dz, mask)
            ffac = 10;
            [Nx, Nz] = size(var);
            if(isscalar(depthforavg))
                depthforavg = ones(Nx,1)*depthforavg;
            end
            if(isscalar(mask))
                mask = ones(Nx,1)*mask;
            end
            % Initialize the output matrix with NaNs
            var_bottom_avg = nan(Nx,1);

            % Loop through each grid point
            for ii = 1:Nx
                % Check if the point is within the valid mask
                if any(mask(ii,:) ~= 0)
                    % interp
                    DZf = zeros(Nz*ffac,1);
                    for k = 1:Nz
                        for nfac = 1:ffac
                            DZf((k-1)*ffac + nfac) = Dz(ii,k)/ffac;
                        end
                    end

                    fine_depth = -cumsum(DZf,'reverse');
                    var_interp = interp1(z2D(ii,:), var(ii,:), fine_depth, 'linear', 'extrap');

                    % Find the indices where the depth is within the specified range
                    ind_zavg = find(fine_depth <= fine_depth(1) + depthforavg(ii));
                    % Calculate the weighted average
                    var_bottom_avg(ii) = sum(var_interp(ind_zavg) .* DZf(ind_zavg)) ...
                        ./ sum(DZf(ind_zavg));
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%% diffxi(rho,dx,z2D_r,z2D_w)
        %%%
        %%% Compute horizontal derivatives using the chain rule
        %%%
        %%% (d/dx)_z = (d/dx)_sigma - [(dz/dx)_sigma]* [(d/dz)]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [rho_dxi,rho_dz]=diffxi(rho,dx,z2D_r,z2D_w)
            [Nx, Nz] = size(rho);
            %%% compute vertical gradient
            rho_w = topo.rho2w2d(rho,z2D_r,z2D_w);
            rho_dz= (rho_w(:,2:end) - rho_w(:,1:end-1))./(z2D_w(:,2:end) - z2D_w(:,1:end-1));
            %%% compute (drho/dx)_sigma
            rho_dx = zeros(Nx, Nz);
            rho_dx(2:Nx-1, :) = (rho(3:Nx, :) - rho(1:Nx-2, :)) / (2 * dx);
            rho_dx(1, :) = (-3*rho(1, :) + 4*rho(2, :) - rho(3, :)) / (2 * dx);
            rho_dx(Nx, :) = (3*rho(Nx, :) - 4*rho(Nx-1, :) + rho(Nx-2, :)) / (2 * dx);
            %%% compute (dz/dx)_sigma
            S_slope = zeros(Nx, Nz);
            S_slope(2:Nx-1, :) = (z2D_r(3:Nx, :) - z2D_r(1:Nx-2, :)) / (2 * dx);
            S_slope(1, :) = (-3*z2D_r(1, :) + 4*z2D_r(2, :) - z2D_r(3, :)) / (2 * dx);
            S_slope(Nx, :) = (3*z2D_r(Nx, :) - 4*z2D_r(Nx-1, :) + z2D_r(Nx-2, :)) / (2 * dx);
            %%% compute (d/dx)_z
            %%% (d/dx)_z = (d/dx)_sigma - [(dz/dx)_sigma]* [(d/dz)]
            rho_dx_sigma_extra= -rho_dz.*S_slope;
            rho_dxi=rho_dx+rho_dx_sigma_extra;
        end


        function [rho_dxi,rho_dz]=diffxi3d(rho,dx,z3D_r,z3D_w)
            [Nx, Ny, Nz] = size(rho);
            %%% compute vertical gradient
            rho_w = topo.rho2w(rho,z3D_r,z3D_w);
            rho_dz= (rho_w(:,:,2:end) - rho_w(:,:,1:end-1))./(z3D_w(:,:,2:end) - z3D_w(:,:,1:end-1));
            %%% compute (drho/dx)_sigma
            rho_dx = zeros(Nx, Ny, Nz);
            rho_dx(2:Nx-1, :, :) = (rho(3:Nx, :, :) - rho(1:Nx-2, :, :)) / (2 * dx);
            rho_dx(1, :, :) = (-3*rho(1, :, :) + 4*rho(2, :, :) - rho(3, :, :)) / (2 * dx);
            rho_dx(Nx, :, :) = (3*rho(Nx, :, :) - 4*rho(Nx-1, :, :) + rho(Nx-2, :, :)) / (2 * dx);
            %%% compute (dz/dx)_sigma
            S_slope = zeros(Nx, Ny,Nz);
            S_slope(2:Nx-1, :, :) = (z3D_r(3:Nx, :, :) - z3D_r(1:Nx-2, :, :)) / (2 * dx);
            S_slope(1, :, :) = (-3*z3D_r(1, :, :) + 4*z3D_r(2, :, :) - z3D_r(3, :, :)) / (2 * dx);
            S_slope(Nx, :, :) = (3*z3D_r(Nx, :, :) - 4*z3D_r(Nx-1, :, :) + z3D_r(Nx-2, :, :)) / (2 * dx);
            %%% compute (d/dx)_z
            %%% (d/dx)_z = (d/dx)_sigma - [(dz/dx)_sigma]* [(d/dz)]
            rho_dx_sigma_extra= -rho_dz.*S_slope;
            rho_dxi=rho_dx+rho_dx_sigma_extra;
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%% diffeta(rho,dy,z2D_r,z2D_w)
        %%%
        %%% Compute horizontal derivatives using the chain rule
        %%%
        %%% (d/dy)_z = (d/dy)_sigma - [(dz/dy)_sigma]* [(d/dz)]
        %%%
        %%% NOTE: y-rows 1 and Ny are left at zero (channel walls) and
        %%% rows 2 and Ny-1 use first-order one-sided differences, unlike
        %%% diffxi which fills its x-boundaries with second-order stencils.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        function [rho_deta,rho_dz] = diffeta(rho,dy,z2D_r,z2D_w)
            [Ny, Nz] = size(rho);
            %%% compute vertical gradient
            rho_w = topo.rho2w2d(rho,z2D_r,z2D_w);
            rho_dz = (rho_w(:,2:end) - rho_w(:,1:end-1)) ./ ...
                (z2D_w(:,2:end) - z2D_w(:,1:end-1));
            %%% compute (drho/dy)_sigma
            rho_dy = zeros(Ny, Nz);
            rho_dy(3:Ny-2, :) = (rho(4:Ny-1, :) - rho(2:Ny-3, :)) / (2 * dy);
            rho_dy(2, :)      = (rho(3, :) - rho(2, :)) / dy;
            rho_dy(Ny-1, :)   = (rho(Ny-1, :) - rho(Ny-2, :)) / dy;
            %%% compute (dz/dy)_sigma
            S_slope = zeros(Ny, Nz);
            S_slope(3:Ny-2, :) = (z2D_r(4:Ny-1, :) - z2D_r(2:Ny-3, :)) / (2 * dy);
            S_slope(2, :)      = (z2D_r(3, :) - z2D_r(2, :)) / dy;
            S_slope(Ny-1, :)   = (z2D_r(Ny-1, :) - z2D_r(Ny-2, :)) / dy;

            %%% compute (d/dy)_z
            rho_dy_sigma_extra = -rho_dz .* S_slope;
            rho_deta = rho_dy + rho_dy_sigma_extra;
        end

        function [rho_deta,rho_dz] = diffeta3d(rho,dy,z3D_r,z3D_w)
            [Nx, Ny, Nz] = size(rho);
            %%% compute vertical gradient
            rho_w = topo.rho2w(rho,z3D_r,z3D_w);
            rho_dz = (rho_w(:,:,2:end) - rho_w(:,:,1:end-1)) ./ ...
                (z3D_w(:,:,2:end) - z3D_w(:,:,1:end-1));
            %%% compute (drho/dy)_sigma
            rho_dy = zeros(Nx, Ny, Nz);
            rho_dy(:,3:Ny-2, :) = (rho(:,4:Ny-1, :) - rho(:,2:Ny-3, :)) / (2 * dy);
            rho_dy(:,2, :)      = (rho(:,3, :) - rho(:,2, :)) / dy;
            rho_dy(:,Ny-1, :)   = (rho(:,Ny-1, :) - rho(:,Ny-2, :)) / dy;
            %%% compute (dz/dy)_sigma
            S_slope = zeros(Nx, Ny, Nz);
            S_slope(:,3:Ny-2, :) = (z3D_r(:,4:Ny-1, :) - z3D_r(:,2:Ny-3, :)) / (2 * dy);
            S_slope(:,2, :)      = (z3D_r(:,3, :) - z3D_r(:,2, :)) / dy;
            S_slope(:,Ny-1, :)   = (z3D_r(:,Ny-1, :) - z3D_r(:,Ny-2, :)) / dy;

            %%% compute (d/dy)_z
            rho_dy_sigma_extra = -rho_dz .* S_slope;
            rho_deta = rho_dy + rho_dy_sigma_extra;
        end

    end
end
