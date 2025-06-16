classdef rt
    % Toolbox for analyzing ROMS simulations
    % organized by Huaiyu
    % V0 August 2024
    methods(Static)


        function w_w = rho2w(var_r, z_r, z_w)
            % Initialize w_w with the same size as z_w
            w_w = zeros(size(z_w));

            % First layer computation
            w_w(:,:,1) = var_r(:,:,1) + (var_r(:,:,2) - var_r(:,:,1)) .* ...
                (z_w(:,:,1) - z_r(:,:,1)) ./ (z_r(:,:,2) - z_r(:,:,1));

            % Middle layers computation
            w_w(:,:,2:end-1) = var_r(:,:,1:end-1) .* (z_r(:,:,2:end) - z_w(:,:,2:end-1)) + ...
                var_r(:,:,2:end)   .* (z_w(:,:,2:end-1) - z_r(:,:,1:end-1));
            w_w(:,:,2:end-1) = w_w(:,:,2:end-1) ./ (z_r(:,:,2:end) - z_r(:,:,1:end-1));
        end


        function w_w = rho2w4d(var_r, z_r, z_w)
            % Initialize w_w with the same size as z_w
            w_w = zeros([size(z_w),size(var_r,4)]);

            % First layer computation
            w_w(:,:,1,:) = var_r(:,:,1,:) + (var_r(:,:,2,:) - var_r(:,:,1,:)) .* ...
                (z_w(:,:,1) - z_r(:,:,1)) ./ (z_r(:,:,2) - z_r(:,:,1));

            % Middle layers computation
            w_w(:,:,2:end-1,:) = var_r(:,:,1:end-1,:) .* (z_r(:,:,2:end) - z_w(:,:,2:end-1)) + ...
                var_r(:,:,2:end,:)   .* (z_w(:,:,2:end-1) - z_r(:,:,1:end-1));
            w_w(:,:,2:end-1,:) = w_w(:,:,2:end-1,:) ./ (z_r(:,:,2:end) - z_r(:,:,1:end-1));
        end

        function var_rho = w2rho(var_w)
            [M, L, N] = size(var_w);

            % Initialize var_rho with the same size as var_w but with one less layer
            var_rho = zeros(M, L, N-1);

            % Compute for the middle layers
            for iz = 2:N-2
                var_rho(:,:,iz) = 0.5625 * (var_w(:,:,iz+1) + var_w(:,:,iz)) - 0.0625 * (var_w(:,:,iz+2) + var_w(:,:,iz-1));
            end

            % Compute for the first and last layers separately
            var_rho(:,:,1) = -0.125 * var_w(:,:,3) + 0.75 * var_w(:,:,2) + 0.375 * var_w(:,:,1);
            var_rho(:,:,end) = -0.125 * var_w(:,:,N-2) + 0.75 * var_w(:,:,N-1) + 0.375 * var_w(:,:,N);
        end




        % Main v2rho Function
        function var_rho = v2rho(var_v)
            ndim = ndims(var_v);
            if ndim == 2
                var_rho = rt.v2rho_2d(var_v);
            elseif ndim == 3
                var_rho = rt.v2rho_3d(var_v);
            else
                var_rho = rt.v2rho_4d(var_v);
            end
        end

        % 2D v2rho Function
        function var_rho = v2rho_2d(var_v)
            [Mp, L] = size(var_v);
            Lp = L + 1;
            Lm = L - 1;
            var_rho = zeros(Mp, Lp);
            var_rho(:, 2:L) = 0.5 * (var_v(:, 1:Lm) + var_v(:, 2:L));
            var_rho(:, 1) = var_rho(:, 2);
            var_rho(:, Lp) = var_rho(:, L);
        end

        % 3D v2rho Function
        function var_rho = v2rho_3d(var_v)
            [Mp, L, N] = size(var_v);
            Lp = L + 1;
            Lm = L - 1;
            var_rho = zeros(Mp, Lp, N);
            var_rho(:, 2:L, :) = 0.5 * (var_v(:, 1:Lm, :) + var_v(:, 2:L, :));
            var_rho(:, 1, :) = var_rho(:, 2, :);
            var_rho(:, Lp, :) = var_rho(:, L, :);
        end

        % 4D v2rho Function
        function var_rho = v2rho_4d(var_v)
            [Mp, L, N, Nt] = size(var_v);
            Lp = L + 1;
            Lm = L - 1;
            var_rho = zeros(Mp, Lp, N, Nt);
            var_rho(:, 2:L, :, :) = 0.5 * (var_v(:, 1:Lm, :, :) + var_v(:, 2:L, :, :));
            var_rho(:, 1, :, :) = var_rho(:, 2, :, :);
            var_rho(:, Lp, :, :) = var_rho(:, L, :, :);
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Main u2rho Function
        function var_rho = u2rho(var_u)
            ndim = ndims(var_u);
            if ndim == 2
                var_rho = rt.u2rho_2d(var_u);
            elseif ndim == 3
                var_rho = rt.u2rho_3d(var_u);
            else
                var_rho = rt.u2rho_4d(var_u);
            end
        end

        % 2D u2rho Function
        function var_rho = u2rho_2d(var_u)
            [M, Lp] = size(var_u);
            Mp = M + 1;
            Mm = M - 1;
            var_rho = zeros(Mp, Lp);
            var_rho(2:M, :) = 0.5 * (var_u(1:Mm, :) + var_u(2:M, :));
            var_rho(1, :) = var_rho(2, :);
            var_rho(Mp, :) = var_rho(M, :);
        end

        % 3D u2rho Function
        function var_rho = u2rho_3d(var_u)
            [M, Lp, N] = size(var_u);
            Mp = M + 1;
            Mm = M - 1;
            var_rho = zeros(Mp, Lp, N);
            var_rho(2:M, :, :) = 0.5 * (var_u(1:Mm, :, :) + var_u(2:M, :, :));
            var_rho(1, :, :) = var_rho(2, :, :);
            var_rho(Mp, :, :) = var_rho(M, :, :);
        end

        % 4D u2rho Function
        function var_rho = u2rho_4d(var_u)
            [M, Lp, N, Nt] = size(var_u);
            Mp = M + 1;
            Mm = M - 1;
            var_rho = zeros(Mp, Lp, N, Nt);
            var_rho(2:M, :, :, :) = 0.5 * (var_u(1:Mm, :, :, :) + var_u(2:M, :, :, :));
            var_rho(1, :, :, :) = var_rho(2, :, :, :);
            var_rho(Mp, :, :, :) = var_rho(M, :, :, :);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Main u2psi Function
        function var_psi = u2psi(var_u)
            ndim = ndims(var_u);
            if ndim == 2
                var_psi = rt.u2psi_2d(var_u);
            elseif ndim == 3
                var_psi = rt.u2psi_3d(var_u);
            else
                var_psi = rt.u2psi_4d(var_u);
            end
        end

        % 2D u2psi Function
        function var_psi = u2psi_2d(var_u)
            var_psi = 0.5 * (var_u(:, 1:end-1) + var_u(:, 2:end));
        end

        % 3D u2psi Function
        function var_psi = u2psi_3d(var_u)
            var_psi = 0.5 * (var_u(:, 1:end-1,:) + var_u(:, 2:end,:));
        end

        % 4D u2psi Function
        function var_psi = u2psi_4d(var_u)
            var_psi = 0.5 * (var_u(:, 1:end-1,:,:) + var_u(:, 2:end,:,:));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Main v2psi Function
        function var_psi = v2psi(var_v)
            ndim = ndims(var_v);
            if ndim == 2
                var_psi = rt.v2psi_2d(var_v);
            elseif ndim == 3
                var_psi = rt.v2psi_3d(var_v);
            else
                var_psi = rt.v2psi_4d(var_v);
            end
        end

        % 2D v2psi Function
        function var_psi = v2psi_2d(var_v)
            var_psi = 0.5 * (var_v(1:end-1,:) + var_v(2:end,:));
        end

        % 3D v2psi Function
        function var_psi = v2psi_3d(var_v)
            var_psi = 0.5 * (var_v(1:end-1,:,:) + var_v(2:end,:,:));
        end

        % 4D v2psi Function
        function var_psi = v2psi_4d(var_v)
            var_psi = 0.5 * (var_v(1:end-1,:,:,:) + var_v(2:end,:,:,:));
        end






        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Main rho2psi Function
        function var_psi = rho2psi(var_rho)
            ndim = ndims(var_rho);
            if ndim < 3
                var_psi = rt.rho2psi_2d(var_rho);
            elseif ndim == 3
                var_psi = rt.rho2psi_3d(var_rho);
            elseif ndim == 4
                var_psi = rt.rho2psi_4d(var_rho);
            end
        end

        % 2D rho2psi Function
        function var_psi = rho2psi_2d(var_rho)
            var_psi = 0.25 * (var_rho(2:end, 2:end) + var_rho(2:end, 1:end-1) + var_rho(1:end-1, 1:end-1) + var_rho(1:end-1, 2:end));
        end

        % 3D rho2psi Function
        function var_psi = rho2psi_3d(var_rho)
            var_psi = 0.25 * (var_rho(2:end, 2:end, :) + var_rho(2:end, 1:end-1, :) + var_rho(1:end-1, 1:end-1, :) + var_rho(1:end-1, 2:end, :));
        end

        % 4D rho2psi Function
        function var_psi = rho2psi_4d(var_rho)
            var_psi = 0.25 * (var_rho(2:end, 2:end, :, :) + var_rho(2:end, 1:end-1, :, :) + var_rho(1:end-1, 1:end-1, :, :) + var_rho(1:end-1, 2:end, :, :));
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Main rho2u Function
        function var_u = rho2u(var_rho)
            ndim = ndims(var_rho);
            if ndim == 1
                var_u = 0.5 * (var_rho(2:end) + var_rho(1:end-1));
            elseif ndim == 2
                var_u = rt.rho2u_2d(var_rho);
            elseif ndim == 3
                var_u = rt.rho2u_3d(var_rho);
            elseif ndim == 4
                var_u = rt.rho2u_4d(var_rho);
            end
        end

        % 2D rho2u Function
        function var_u = rho2u_2d(var_rho)
            var_u = 0.5 * (var_rho(2:end, :) + var_rho(1:end-1, :));
        end

        % 3D rho2u Function
        function var_u = rho2u_3d(var_rho)
            var_u = 0.5 * (var_rho(2:end, :, :) + var_rho(1:end-1, :, :));
        end

        % 4D rho2u Function
        function var_u = rho2u_4d(var_rho)
            var_u = 0.5 * (var_rho(2:end, :, :, :) + var_rho(1:end-1, :, :, :));
        end




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Main rho2v Function
        function var_v = rho2v(var_rho)
            ndim = ndims(var_rho);
            if ndim == 1
                var_v = 0.5 * (var_rho(2:end) + var_rho(1:end-1));
            elseif ndim == 2
                var_v = rt.rho2v_2d(var_rho);
            elseif ndim == 3
                var_v = rt.rho2v_3d(var_rho);
            elseif ndim == 4
                var_v = rt.rho2v_4d(var_rho);

            end
        end

        % 2D rho2v Function
        function var_v = rho2v_2d(var_rho)
            var_v = 0.5 * (var_rho(:, 2:end) + var_rho(:, 1:end-1));
        end

        % 3D rho2v Function
        function var_v = rho2v_3d(var_rho)
            var_v = 0.5 * (var_rho(:, 2:end, :) + var_rho(:, 1:end-1, :));
        end

        % 4D rho2v Function
        function var_v = rho2v_4d(var_rho)
            var_v = 0.5 * (var_rho(:, 2:end, :, :) + var_rho(:, 1:end-1, :, :));
        end




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Main rho2uvp Function
        function [ufield,vfield,pfield]=rho2uvp(rfield)
            ndim = ndims(rfield);
            if ndim == 2
                [ufield,vfield,pfield]=rt.rho2uvp_2d(rfield);
            elseif ndim == 3
                [ufield,vfield,pfield]=rt.rho2uvp_3d(rfield);
            elseif ndim == 4
                [ufield,vfield,pfield]=rt.rho2uvp_4d(rfield);
            end
        end

        % 2D rho2uvp Function
        function  [ufield,vfield,pfield]=rho2uvp_2d(rfield)
            [Nxp1,Nyp1]=size(rfield);
            Nx=Nxp1-1;
            Ny=Nyp1-1;
            %
            vfield=0.5*(rfield(:,1:Ny)+rfield(:,2:Nyp1));
            ufield=0.5*(rfield(1:Nx,:)+rfield(2:Nxp1,:));
            pfield=0.5*(ufield(:,1:Ny)+ufield(:,2:Nyp1));
        end

        % 3D rho2uvp Function
        function [ufield,vfield,pfield]=rho2uvp_3d(rfield)
            [Nxp1,Nyp1,N]=size(rfield);
            Nx=Nxp1-1;
            Ny=Nyp1-1;
            vfield=0.5*(rfield(:,1:Ny,:)+rfield(:,2:Nyp1,:));
            ufield=0.5*(rfield(1:Nx,:,:)+rfield(2:Nxp1,:,:));
            pfield=0.5*(ufield(:,1:Ny,:)+ufield(:,2:Nyp1,:));
        end

        % 4D rho2uvp Function
        function [ufield,vfield,pfield]=rho2uvp_4d(rfield)
            [Nxp1,Nyp1,~,~]=size(rfield);
            Nx=Nxp1-1;
            Ny=Nyp1-1;
            vfield=0.5*(rfield(:,1:Ny,:,:)+rfield(:,2:Nyp1,:,:));
            ufield=0.5*(rfield(1:Nx,:,:,:)+rfield(2:Nxp1,:,:,:));
            pfield=0.5*(ufield(:,1:Ny,:,:)+ufield(:,2:Nyp1,:,:));
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
            %
            %  pmask(1:L,1:M)=rmask(1:L,1:M ).*rmask(2:Lp,1:M ).*                   ...
            %                 rmask(1:L,2:Mp).*rmask(2:Lp,2:Mp);

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
        % function for calculating horizontal gradients of a 2D variable
        function  [var_dx,var_dy] = grad_2d(var,dx,dy)
            % Get the size of the input variable
            [Nx, Ny] = size(var);
            % Initialize the output variables
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
        %function for find the output file names
        function filelist = findoutput(exppath,outputName)
            % Generate the pattern to match files
            filelist_all = dir(fullfile(exppath, [outputName '.*.nc']));
            % Initialize the filelist array
            filelist = {};
            ind = 1;
            % Loop through all matching files and store their full paths
            for i = 1:numel(filelist_all)
                filelist{ind} = fullfile(exppath, filelist_all(i).name);
                ind = ind +1;
            end
            % Sort the filelist
            filelist = sort(filelist);
        end
        function ind_avg = find_avgind(diag_time)

            if(max(diag_time) > 3500)
                ind_avg = find(diag_time>=3100 & diag_time<=3600);
                ind_avg = [length(diag_time)-4:1:length(diag_time)];
            else
                ind_avg = find(diag_time>=1300 & diag_time<=1800) ;
                if(length(ind_avg)<5)
                ind_avg = [length(diag_time)-4:1:length(diag_time)];
                end
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %functions for loading varibales from netcdf files
        function [temp_u,temp_v] = load_nc_uv(nc,Var_str,fillvalue,ind_avg)

            temp_u = nc{['u_',Var_str]}(ind_avg,:,:,:);
            temp_u = permute(temp_u,[4 3 2 1]);
            temp_u(temp_u==fillvalue) = nan;
            temp_u(temp_u==0) = nan;

            temp_v = nc{['v_',Var_str]}(ind_avg,:,:,:);
            temp_v = permute(temp_v,[4 3 2 1]);
            temp_v(temp_v==fillvalue) = nan;
            temp_v(temp_v==0) = nan;

        end



        function [temp] = load_nc_3d(nc,Var_str,fillvalue,ind_avg)


            temp = nc{Var_str}(ind_avg,:,:,:);
            temp = permute(temp,[4 3 2 1]);
            temp(temp==fillvalue) = nan;
            temp(temp==0) = nan;


        end



        function [temp] = load_nc_2d(nc,Var_str,fillvalue,ind_avg)


            temp = nc{Var_str}(ind_avg,:,:);
            temp = permute(temp,[3 2 1]);
            temp(temp==fillvalue) = nan;
            temp(temp==0) = nan;


        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %function for calculating vorticity; u and v are located at U- and V-
        %points, respectively. The resulting vorticity is at Psi-point.

        function curl = calc_curl (u,v,dx,dy)

            curl = - (u(:,2:end,:)-u(:,1:end-1,:))./dy;
            curl = curl  + (v(2:end,:,:)-v(1:end-1,:,:))./dx;

        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %function for calculating divergence; u and v are located at U- and V-
        %points, respectively. The resulting divergence is at Rho-point.

        function div = calc_div (u,v,dx,dy)
            div = zeros(size(v,1),size(u,2),size(u,3));
            div(2:end-1,2:end-1,:) =  (u(2:end,2:end-1,:)-u(1:end-1,2:end-1,:))./dx ...
                + (v(2:end-1,2:end,:)-v(2:end-1,1:end-1,:))./dy;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%% calc_curl_zlev
        %%%
        %%% Convenience function to calculate curl of a vector by first
        %%% interpolating it from sigma levels to z levles
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        function [curl_var,curl_var_zint] = calc_curl_zlev(nc,varname,fillvalue,ind_avg,Dz_u,Dz_v,z3D_u,z3D_v,mask_u,mask_v,zz,dx,dy)
            [temp_u, temp_v] = rt.load_nc_uv(nc,varname,fillvalue,ind_avg);
            temp_u_zint = squeeze(sum(mean(temp_u,4),3)); temp_v_zint = squeeze(sum(mean(temp_v,4),3));
            temp_u = mean(temp_u./Dz_u,4); temp_v = mean(temp_v./Dz_v,4);
            [temp_u_zlev,temp_v_zlev] = rt.interp2z_uv(temp_u,temp_v,z3D_u,z3D_v,mask_u,mask_v,zz);
            curl_var = rt.calc_curl(temp_u_zlev,temp_v_zlev,dx,dy);
            curl_var_zint = rt.calc_curl(temp_u_zint,temp_v_zint,dx,dy);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% calc_iso_int
        %%%
        %%% Convenience function to calculate along-isobath integrals
        %%% on isobaths specified by the vector dd (positive).
        %%% x and y gridspacings are dx and dy.
        %%% Assumes that the data (Nx-1 x Ny-1 x Nlay) and the bathymetric
        %%% elevation etab (Nx-1 x Ny-1) are on vorticity points.
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function data_int = calc_iso_int (data,dx,dy,dd,etab,mask,ffac)

            %%% Grid refinement factor
            dx_f = dx/ffac;
            dy_f = dy/ffac;

            %%% Grid dimensions
            Nx = size(data,1)-1; %%% WHY hacks: number of interior rho points
            Ny = size(data,2)-1;
            Nlay = size(data,3);
            Nd = length(dd);
            Lx = Nx*dx;
            Ly = Ny*dy;

            %%% Grids for interpolation
            xx = 0:dx:Lx;
            yy = 0:dy:Ly;%%% WHY hacks: inclue point at y=Ly
            xx_f = 0:dx_f:Lx; %%% WHY hacks: inclue point at x=Lx
            yy_f = 0:dy_f:Ly; %%% WHY hacks: inclue point at y=Ly
            Nx_f = length(xx_f);
            Ny_f = length(yy_f);
            [XX,YY] = meshgrid(xx,yy);
            [XX_f,YY_f] = meshgrid(xx_f,yy_f);

            %%% WHY hacks: remove periodicity
            %%% Extend input grid to account for periodicity
            % data = [data ; data(1,:,:)];
            % etab = [etab ; etab(1,:,:)];

            %%% Interpolate data onto fine grid
            data_f = zeros(Nx_f,Ny_f,Nlay);
            for k=1:Nlay
                data_f(:,:,k) = interp2(XX,YY,data(:,:,k)',XX_f,YY_f,'linear')';
            end
            etab_f = interp2(XX,YY,etab',XX_f,YY_f,'linear')';
            mask_f = interp2(XX,YY,mask',XX_f,YY_f,'linear')';

            %%% Calculate contour lengths
            data_int = zeros(Nd,Nlay);
            for n=1:Nd
                % msk_dxc = repmat(xor(-etab_f(:,1:Ny_f-1)>dd(n),-etab_f(:,2:Ny_f)>dd(n)), [1 1 Nlay]);
                msk_dxc = repmat( and( xor(-etab_f(:,1:Ny_f-1)>dd(n),-etab_f(:,2:Ny_f)>dd(n)), ...
                    mask_f(:,1:Ny_f-1).*mask_f(:,2:Ny_f)~=0 ), [1 1 Nlay]);

                % msk_dyc = repmat( xor(-etab_f(1:Nx_f-1,:)>dd(n),-etab_f(2:Nx_f,:)>dd(n)), [1 1 Nlay]);
                msk_dyc = repmat( and(xor(-etab_f(1:Nx_f-1,:)>dd(n),-etab_f(2:Nx_f,:)>dd(n)), ...
                    mask_f(1:Nx_f-1,:).*mask_f(2:Nx_f,:)~=0), [1 1 Nlay]);
                % figure;pcolor(msk_dxc'); figure;pcolor(msk_dyc')
                data_dxc = 0.5*(data_f(:,1:Ny_f-1,:)+data_f(:,2:Ny_f,:)); %WHY hacks: remove periodicity
                data_dyc = 0.5*(data_f(1:Nx_f-1,:,:)+data_f(2:Nx_f,:,:));
                data_int(n,:) = sum(sum(dx_f.*msk_dxc.*data_dxc,1),2) + sum(sum(dy_f.*msk_dyc.*data_dyc,1),2);
            end

        end


        %%% convinient function when var is on rho point
        function   var_isoba = calc_iso_int_rho(var_rho,dx,dy,dd,coord_q,mask,interpFac,cntrlen)
            [~,~,var_q] = rt.rho2uvp_3d(var_rho);
            var_q(isnan(var_q))=0;
            var_isoba = rt.calc_iso_int(var_q,dx,dy,dd,coord_q,mask,interpFac) ./ cntrlen;
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% calc_iso_int
        %%%
        %%% Convenience function to calculate along-isobath integrals
        %%% on isobaths specified by the vector dd (positive).
        %%% x and y gridspacings are dx and dy.
        %%% Assumes that the data (Nx-1 x Ny-1 x Nlay) and the bathymetric
        %%% elevation etab (Nx-1 x Ny-1) are on vorticity points.
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [data_avg,data_remap_rho,data_bin_rho] = calc_iso_int_remap(data,dx,dy,dd,etab,mask,ffac,cntrlen)

            if(ffac~=1)
                error("ffac must equal 1 for remaping isobath-averaged quantities")
            end
            Nx = size(data,1)-1;
            Ny = size(data,2)-1
            Nlay = size(data,3);
            data_bin = nan(Nx,Ny);
            data_remap = nan(Nx,Ny,Nlay);
            %%% Grid refinement factor
            dx_f = dx;
            dy_f = dy;


            Nd = length(dd);

            %%% Interpolate data onto fine grid
            data_f = data;
            etab_f = etab;
            mask_f = mask;

            %%% Calculate contour lengths
            data_avg = zeros(Nd,Nlay);
            for n=1:Nd
                % msk_dxc = repmat(xor(-etab_f(:,1:Ny_f-1)>dd(n),-etab_f(:,2:Ny_f)>dd(n)), [1 1 Nlay]);
                msk_dxc = repmat( and( xor(-etab_f(:,1:end-1)>dd(n),-etab_f(:,2:end)>dd(n)), ...
                    mask_f(:,1:end-1).*mask_f(:,2:end)~=0 ), [1 1 Nlay]);
                % msk_dyc = repmat( xor(-etab_f(1:end-1,:)>dd(n),-etab_f(2:end,:)>dd(n)), [1 1 Nlay]);
                msk_dyc = repmat( and(xor(-etab_f(1:end-1,:)>dd(n),-etab_f(2:end,:)>dd(n)), ...
                    mask_f(1:end-1,:).*mask_f(2:end,:)~=0), [1 1 Nlay]);
                % figure;pcolor(msk_dxc'); figure;pcolor(msk_dyc')
                data_dxc = 0.5*(data_f(:,1:end-1,:)+data_f(:,2:end,:)); %WHY hacks: remove periodicity
                data_dyc = 0.5*(data_f(1:end-1,:,:)+data_f(2:end,:,:));
                data_avg(n,:) = (sum(sum(dx_f.*msk_dxc.*data_dxc,1),2) + sum(sum(dy_f.*msk_dyc.*data_dyc,1),2))./cntrlen(n);

                data_rho = repmat(reshape(data_avg(n,:),[1,1,Nlay]),[Nx,Ny,1]);
                msk_dxc = squeeze(msk_dxc(:,:,1)); msk_dyc = squeeze(msk_dyc(:,:,1));
                msk_rho = msk_dxc(1:end-1,:) + msk_dxc(2:end,:) + msk_dyc(:,1:end-1) + msk_dyc(:,2:end) >= 2;
                if(~isnan(data_bin(msk_rho)))
                    disp(['isobath ',num2str(n)])
                    error('one grid is assigned to two isobaths')
                else
                    data_bin(msk_rho)= n ;
                    msk_rho = repmat(msk_rho, [1 1 Nlay]);
                    data_remap(msk_rho)= data_rho(msk_rho);
                end
            end

            data_remap_rho = nan(Nx+2,Ny+2,Nlay);
            data_remap_rho(2:end-1,2:end-1,:) = data_remap;
            data_bin_rho = nan(Nx+2,Ny+2);
            data_bin_rho(2:end-1,2:end-1)=data_bin;
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%% calc_iso_circ
        %%%
        %%% Convenience function to calculate along-isobath circulation on
        %%% isobaths specified by the vector dd (positive).
        %%% x and y gridspacings are dx and dy.
        %%% Assumes that the curl (Nx-1 x Ny-1 x Nlay) and the bathymetric elevation
        %%% etab (Nx-1 x Ny-1) are on vorticity points.
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function circ = calc_iso_circ (curl,dx,dy,dd,etab,mask,ffac)

            %%% Grid refinement factor
            dx_f = dx/ffac;
            dy_f = dy/ffac;

            %%% Grid dimensions
            Nx = size(curl,1)-1; %%% WHY hacks: Nx now is the number of rho points -2
            Ny = size(curl,2)-1;
            Nlay = size(curl,3);
            Nd = length(dd);
            Lx = Nx*dx;
            Ly = Ny*dy;

            %%% Grids for interpolation
            xx = 0:dx:Lx;
            yy = 0:dy:Ly;%%% WHY hacks: inclue point at y=Ly
            xx_f = 0:dx_f:Lx; %%% WHY hacks: inclue point at x=Lx
            yy_f = 0:dy_f:Ly; %%% WHY hacks: inclue point at y=Ly
            Nx_f = length(xx_f);
            Ny_f = length(yy_f);
            [XX,YY] = meshgrid(xx,yy);
            [XX_f,YY_f] = meshgrid(xx_f,yy_f);

            %%% WHY hacks: remove periodicity
            %%% Extend input grid to account for periodicity
            % curl = [curl ; curl(1,:,:)];
            % etab = [etab ; etab(1,:,:)];

            %%% Interpolate curl onto fine grid
            curl_f = zeros(Nx_f,Ny_f,Nlay);
            for k=1:Nlay
                curl_f(:,:,k) = interp2(XX,YY,curl(:,:,k)',XX_f,YY_f,'linear')';
            end
            etab_f = interp2(XX,YY,etab',XX_f,YY_f,'linear')';
            mask_f = interp2(XX,YY,mask',XX_f,YY_f,'linear')';
            %%% Perform area intergrals to compute circulation
            circ = zeros(Nd,Nlay);


            msk = repmat( mask_f==0, [1 1 Nlay]);
            circ_mid = sum(sum(curl_f.*msk*dx_f*dy_f,1),2);

            for n = 1:Nd
                msk = repmat( and(etab_f<-dd(n),mask_f~=0),[1 1 Nlay]);
                circ(n,:) = circ_mid+ sum(sum(curl_f.*msk*dx_f*dy_f,1),2);
            end


            for nn = 1:Nlay
                temp = circ(2:end,nn) - circ(1:end-1,nn);
                temp = find(temp ==0,1)+1;
                circ(temp:end,nn)=0;
            end

        end

        %%% convinient function when var is on rho point




        function   circ_val = calc_iso_circ_pack(curl_var,dx,dy,dd,coord_q,mask,interpFac,cntrlen)
            Nz = size(curl_var,3);
            curl_var(isnan(curl_var)) = 0;
            curl_var(isinf(curl_var)) = 0;
            circ_val = rt.calc_iso_circ(curl_var,dx,dy,dd,coord_q,mask,interpFac) ./ repmat(cntrlen,[1 Nz]);
            circ_val(isinf(circ_val)) = nan;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%% interp2z
        %%%
        %%% Convenience function to interpolate variables on sigma levels
        %%% to z levels
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [u_zlev,v_zlev] = interp2z_uv(u,v,z_u,z_v,mask_u,mask_v,z_new)

            Nx = size(v,1);
            Ny = size(u,2);
            Nz = length(z_new);

            u_zlev = zeros(Nx-1,Ny,Nz)  ;
            for ii = 1:Nx-1
                for jj = 1 : Ny
                    if(mask_u(ii,jj)~=0)
                        u_zlev(ii,jj,:) = interp1(squeeze(z_u(ii,jj,:)),squeeze(u(ii,jj,:)),z_new);
                    else
                        u_zlev(ii,jj,:) =nan;
                    end
                end
            end

            v_zlev = zeros(Nx,Ny-1,Nz)  ;
            for ii = 1:Nx
                for jj = 1 : Ny -1
                    if(mask_v(ii,jj)~=0)
                        v_zlev(ii,jj,:) = interp1(squeeze(z_v(ii,jj,:)),squeeze(v(ii,jj,:)),z_new);
                    else
                        v_zlev(ii,jj,:) =nan;
                    end
                end
            end

        end


        function var_zlev = interp2z(var,z,mask,z_new)

            Nx = size(var,1);
            Ny = size(var,2);
            Nz = length(z_new);

            var_zlev = zeros(Nx,Ny,Nz)  ;
            for ii = 1:Nx
                for jj = 1 : Ny
                    if(mask(ii,jj)~=0)
                        var_zlev(ii,jj,:) = interp1(squeeze(z(ii,jj,:)),squeeze(var(ii,jj,:)),z_new);
                    else
                        var_zlev(ii,jj,:) =nan;
                    end
                end
            end



        end


        function var_zlev = interp2z_2d(var,z,z_new)

            Nx = size(var,1);
            Nz = length(z_new);

            var_zlev = zeros(Nx,Nz)  ;
            for ii = 1:Nx

                var_zlev(ii,:) = interp1(squeeze(z(ii,:)),squeeze(var(ii,:)),z_new);

            end



        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%% isobaProj
        %%%
        %%% Convenience function to project momentum budget to
        %%% along-isobath direction
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        function [var_alongisob, var_alongisob_zint] = isobaProj(nc, varname, fillvalue, ind_avg, Dz_u, Dz_v, ax, ay)
            % Load u and v data from netCDF file
            [temp_u, temp_v] = rt.load_nc_uv(nc, varname , fillvalue, ind_avg);

            % Calculate the vertically integrated components
            temp_u_zint = squeeze(sum(mean(temp_u, 4), 3));
            temp_v_zint = squeeze(sum(mean(temp_v, 4), 3));

            % Calculate the mean u and v components divided by their respective cell thickness
            temp_u = mean(temp_u ./ Dz_u, 4);
            temp_v = mean(temp_v ./ Dz_v, 4);

            % Project u and v onto isobath directions
            var_alongisob = rt.u2rho(temp_u) .* ax + rt.v2rho(temp_v) .* ay;
            var_alongisob_zint = rt.u2rho(temp_u_zint) .* ax + rt.v2rho(temp_v_zint) .* ay;
        end


        function [var_alongisob, var_alongisob_zint,var_crossisob, var_crossisob_zint] = isocrossProj(nc, varname, fillvalue, ind_avg, Dz_u, Dz_v, ax, ay, cx, cy)
            % Load u and v data from netCDF file
            [temp_u, temp_v] = rt.load_nc_uv(nc, varname , fillvalue, ind_avg);

            % Calculate the vertically integrated components
            temp_u_zint = squeeze(sum(mean(temp_u, 4), 3));
            temp_v_zint = squeeze(sum(mean(temp_v, 4), 3));

            % Calculate the mean u and v components divided by their respective cell thickness
            temp_u = mean(temp_u ./ Dz_u, 4);
            temp_v = mean(temp_v ./ Dz_v, 4);

            % Project u and v onto along isobath directions
            var_alongisob = rt.u2rho(temp_u) .* ax + rt.v2rho(temp_v) .* ay;
            var_alongisob_zint = rt.u2rho(temp_u_zint) .* ax + rt.v2rho(temp_v_zint) .* ay;

            % Project u and v onto corss isobath directions
            var_crossisob = rt.u2rho(temp_u) .* cx + rt.v2rho(temp_v) .* cy;
            var_crossisob_zint = rt.u2rho(temp_u_zint) .* cx + rt.v2rho(temp_v_zint) .* cy;            
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%% do_bottom_avg
        %%%
        %%% Conduct avergaed within certain depth about the bottom
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        function var_bottom_avg = do_bottom_avg(var, depthforavg, z3D, Dz, mask)
            % Function to compute bottom average of topostrophy up to a specified depth

            % Get the dimensions of the input matrices
            [Nx, Ny, ~] = size(var);

            % Initialize the output matrix with NaNs
            var_bottom_avg = nan(Nx, Ny);

            % Loop through each grid point
            for ii = 1:Nx
                for jj = 1:Ny
                    % Check if the point is within the valid mask
                    if mask(ii,jj) ~= 0
                        % Find the indices where the depth is within the specified range
                        ind_zavg = find(z3D(ii,jj,:) <= z3D(ii,jj,1) - 0.5* Dz(ii,jj,1) + depthforavg);
                        % Calculate the weighted average of topostrophy
                        var_bottom_avg(ii,jj) = sum(var(ii,jj,ind_zavg) .* Dz(ii,jj,ind_zavg), 3) ...
                            ./ sum(Dz(ii,jj,ind_zavg), 3);
                    end
                end
            end
        end
























    end
end
