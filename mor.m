classdef mor
    % Toolbox for analyzing ROMS simulations
    % organized by Huaiyu
    % V2 Nov 2024
    methods(Static)

        % function [z,Cs,sc] = zlevs4(h,zeta,theta_s,theta_b,hc,N,type,scoord_type,alpha,beta)
        function [z,Cs] = zlevs4(h,zeta,theta_s,theta_b,hc,N,type,scoord_type,alpha,beta)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  function z = zlevs4(h,zeta,theta_s,theta_b,hc,N,type,scoord_type,alpha,beta)
            %
            %  this function compute the depth of rho or w points for ROMS
            %
            %  On Input:
            %
            %    type    'r': rho point 'w': w point
            %    scoord_type     : 'old1994' (Song, 1994),
            %                 'new2006' (Sasha, 2006)
            %                 'new2008' bottom stretching included (Sasha, 2008)
            %                 'new2012' latest version
            %    alpha,beta : optional, used for 'new2008'-type s-coordinate
            %
            %  On Output:
            %
            %    z       Depths (m) of RHO- or W-points (3D matrix).
            %
            %  Further Information:
            %  http://www.brest.ird.fr/Roms_tools/
            %
            %  This file is part of ROMSTOOLS
            %
            %  ROMSTOOLS is free software; you can redistribute it and/or modify
            %  it under the terms of the GNU General Public License as published
            %  by the Free Software Foundation; either version 2 of the License,
            %  or (at your option) any later version.
            %
            %  ROMSTOOLS is distributed in the hope that it will be useful, but
            %  WITHOUT ANY WARRANTY; without even the implied warranty of
            %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
            %  GNU General Public License for more details.
            %
            %  You should have received a copy of the GNU General Public License
            %  along with this program; if not, write to the Free Software
            %  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
            %  MA  02111-1307  USA
            %
            %  Copyright (c) 2002-2006 by Pierrick Penven
            %  e-mail:Pierrick.Penven@ird.fr
            %
            %  modified by Yusuke Uchiyama, UCLA, 2008
            %  further modified by Evan Mason, UCLA, 2008
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if nargin<8
                error('Not enough input arguments')
            else
                if nargin < 9
                    %disp('------ default values taken: alpha=0, beta=1')
                    alpha=0; beta=1;
                end
                %
                %   if (strcmp(scoord_type,'new2008'))
                %       disp('--- using new s-coord (2008)')
                %   elseif (strcmp(scoord_type,'new2012'))
                %       disp('--- using new s-coord (2012)')
                %   elseif (strcmp(scoord_type,'new2006'))
                %       disp('--- using new s-coord (2006)')
                %   end
            end
            %
            [M, L] = size(h);
            %
            % Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
            %
            if (strcmp(type,'w'))
                sc = ((0:N) - N) / N;
                N = N + 1;
            else
                sc=((1:N)-N-0.5) / N;
            end

            if (strcmp(scoord_type,'new2008') || strcmp(scoord_type,'new2012'))
                % new s-coordinate allowing smooth bottom refinement
                % alpha=-1: return to pure surface s-coord; -1 < alpha <\infty: bottom refinement
                % beta:
                if (theta_b>0);
                    csrf = (1.0-cosh(theta_s*sc))./(cosh(theta_s)-1.0);
                    if (strcmp(scoord_type,'new2012'))
                        Cs = (exp(theta_b*csrf)-1)./(1.-exp(-theta_b));
                    else % new2008
                        x = sc+1.0;
                        wgt  = x.^alpha./beta.*(alpha+beta-alpha.*x.^beta);
                        cbot = sinh(theta_b*x)./sinh(theta_b)-1.0;
                        Cs   = wgt.*csrf+(1.0-wgt).*cbot;
                    end
                else
                    Cs   = (1.0-cosh(theta_s*sc))./(cosh(theta_s)-1.0);
                end
            else
                % for 'old' and 'new2006' s-coordinate
                cff1 = 1./sinh(theta_s);
                cff2 = 0.5/tanh(0.5*theta_s);
                Cs = (1.-theta_b) * cff1 * sinh(theta_s * sc)...
                    + theta_b * (cff2 * tanh(theta_s * (sc + 0.5)) - 0.5);
            end
            %
            % Create S-coordinate system: based on model topography h(i,j),
            % fast-time-averaged free-surface field and vertical coordinate
            % transformation metrics compute evolving depths of of the three-
            % dimensional model grid.
            %
            z=zeros(N,M,L);
            if (strcmp(scoord_type,'old1994'))
                disp('--- using old s-coord')
                if ( hc > 12 || hc < 5 )
                    error ([mfilename ':grdchk'], '\n\t FATAL (%s) : For old s-coord, current tcline(hc) value (%g) is wrong.\n',mfilename,hc)
                end
                hinv=1./h;
                cff=hc*(sc-Cs);
                cff1=Cs;
                %cff2=sc+1;
                for k=1:N
                    z0=cff(k)+cff1(k)*h;
                    z(k,:,:)=z0+zeta.*(1.+z0.*hinv);
                end
            else
                %if (strcmp(scoord_type,'new2006'))
                %disp(['--- using new s-coord (', scoord_type, ')'])
                %end
                %  if ( hc > 125 || hc < 110 )
                %     error ([mfilename ':grdchk'], '\n\t FATAL (%s) : For new s-coord, current tcline(hc) value (%g) is wrong.\n',mfilename,hc)
                %  end
                hinv=1./(h+hc);
                cff=hc*sc;
                cff1=Cs;
                for k=1:N
                    z(k,:,:)=zeta+(zeta+h).*(cff(k)+cff1(k)*h).*hinv;
                end
            end

        end

        function var_w = rho2w(var_r, z_r, z_w)
            % Initialize w_w with the same size as z_w
            var_w = zeros(size(z_w));

            % bottom layer computation
            var_w(:,:,1) = var_r(:,:,1) + (var_r(:,:,2) - var_r(:,:,1)) .* ...
                (z_w(:,:,1) - z_r(:,:,1)) ./ (z_r(:,:,2) - z_r(:,:,1));

            % surface layer computation
            var_w(:,:,end) = var_r(:,:,end) + (var_r(:,:,end) - var_r(:,:,end-1)) .* ...
                (z_w(:,:,end) - z_r(:,:,end)) ./ (z_r(:,:,end) - z_r(:,:,end-1));

            % Middle layers computation
            var_w(:,:,2:end-1) = var_r(:,:,1:end-1) .* (z_r(:,:,2:end) - z_w(:,:,2:end-1)) + ...
                var_r(:,:,2:end)   .* (z_w(:,:,2:end-1) - z_r(:,:,1:end-1));
            var_w(:,:,2:end-1) = var_w(:,:,2:end-1) ./ (z_r(:,:,2:end) - z_r(:,:,1:end-1));
        end

        function var_w = rho2w2d(var_r, z_r, z_w)
            % Initialize w_w with the same size as z_w
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

        function var_w = rho2w4d(var_r, z_r, z_w)
            % Initialize w_w with the same size as z_w
            var_w = zeros([size(z_w),size(var_r,4)]);

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

        function var_rho = w2rho2d(var_w)
            [M,  N] = size(var_w);

            % Initialize var_rho with the same size as var_w but with one less layer
            var_rho = zeros(M,  N-1);

            % Compute for the middle layers
            for iz = 2:N-2
                var_rho(:,iz) = 0.5625 * (var_w(:,iz+1) + var_w(:,iz)) - 0.0625 * (var_w(:,iz+2) + var_w(:,iz-1));
            end

            % Compute for the first and last layers separately
            var_rho(:,1) = -0.125 * var_w(:,3) + 0.75 * var_w(:,2) + 0.375 * var_w(:,1);
            var_rho(:,end) = -0.125 * var_w(:,N-2) + 0.75 * var_w(:,N-1) + 0.375 * var_w(:,N);
        end
        function var_rho = w2rho4d(var_w)
            [M, L, N,TT] = size(var_w);

            % Initialize var_rho with the same size as var_w but with one less layer
            var_rho = zeros(M, L, N-1,TT);

            % Compute for the middle layers
            for iz = 2:N-2
                var_rho(:,:,iz,:) = 0.5625 * (var_w(:,:,iz+1,:) + var_w(:,:,iz,:)) - 0.0625 * (var_w(:,:,iz+2,:) + var_w(:,:,iz-1,:));
            end

            % Compute for the first and last layers separately
            var_rho(:,:,1,:) = -0.125 * var_w(:,:,3,:) + 0.75 * var_w(:,:,2,:) + 0.375 * var_w(:,:,1,:);
            var_rho(:,:,end,:) = -0.125 * var_w(:,:,N-2,:) + 0.75 * var_w(:,:,N-1,:) + 0.375 * var_w(:,:,N,:);
        end

        % function var_rho = v2rho_2d_yz(var_v)
        %     [L , Z] = size(var_v);
        %     Lp = L + 1;
        %     Lm = L - 1;
        %     var_rho = zeros(Lp,Z);
        %     var_rho(2:L,:) = 0.5 * (var_v(1:Lm,:) + var_v(2:L,:));
        %     var_rho(1,:) = var_rho(2,:);
        %     var_rho(Lp,:) = var_rho(L,:);
        % end



        % Main v2rho Function
        function var_rho = v2rho(var_v)
            ndim = ndims(var_v);
            if ndim == 2
                var_rho = mor.v2rho_2d(var_v);
            elseif ndim == 3
                var_rho = mor.v2rho_3d(var_v);
            else
                var_rho = mor.v2rho_4d(var_v);
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
                var_rho = mor.u2rho_2d(var_u);
            elseif ndim == 3
                var_rho = mor.u2rho_3d(var_u);
            else
                var_rho = mor.u2rho_4d(var_u);
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
                var_psi = mor.u2psi_2d(var_u);
            elseif ndim == 3
                var_psi = mor.u2psi_3d(var_u);
            else
                var_psi = mor.u2psi_4d(var_u);
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
                var_psi = mor.v2psi_2d(var_v);
            elseif ndim == 3
                var_psi = mor.v2psi_3d(var_v);
            else
                var_psi = mor.v2psi_4d(var_v);
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
                var_psi = mor.rho2psi_2d(var_rho);
            elseif ndim == 3
                var_psi = mor.rho2psi_3d(var_rho);
            elseif ndim == 4
                var_psi = mor.rho2psi_4d(var_rho);
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
                var_u = mor.rho2u_2d(var_rho);
            elseif ndim == 3
                var_u = mor.rho2u_3d(var_rho);
            elseif ndim == 4
                var_u = mor.rho2u_4d(var_rho);
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

        % 2D rho2v Function
        % function var_v = rho2v_2dyz(var_rho)
        %     var_v = 0.5 * (var_rho(2:end,:) + var_rho(1:end-1,:));
        % end

        function var_v = rho2v(var_rho)
            ndim = ndims(var_rho);
            if ndim == 1
                var_v = 0.5 * (var_rho(2:end) + var_rho(1:end-1));
            elseif ndim == 2
                var_v = mor.rho2v_2d(var_rho);
            elseif ndim == 3
                var_v = mor.rho2v_3d(var_rho);
            elseif ndim == 4
                var_v = mor.rho2v_4d(var_rho);

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
                [ufield,vfield,pfield]=mor.rho2uvp_2d(rfield);
            elseif ndim == 3
                [ufield,vfield,pfield]=mor.rho2uvp_3d(rfield);
            elseif ndim == 4
                [ufield,vfield,pfield]=mor.rho2uvp_4d(rfield);
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
            [Nxp1,Nyp1,~]=size(rfield);
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
            if isempty(filelist_all)
                filelist = cell(0,1);
                return;
            end

            % Build a sorted cell array of full file paths.
            filenames = sort({filelist_all.name});
            filelist = fullfile(exppath, filenames(:));
        end

        % function ind_avg = find_avgind(diag_time)
        %
        %     if(max(diag_time) > 3500)
        %         ind_avg = find(diag_time>=3100 & diag_time<=3600);
        %         ind_avg = [length(diag_time)-4:1:length(diag_time)];
        %     else
        %         ind_avg = find(diag_time>=1300 & diag_time<=1800) ;
        %         if(length(ind_avg)<5)
        %         ind_avg = [length(diag_time)-4:1:length(diag_time)];
        %         end
        %     end
        % end


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
            if ismatrix(temp)
                temp = permute(temp,[2 1]);
            elseif (ndims(temp)==3)
                temp = permute(temp,[3 2 1]);
            else
                temp = permute(temp,[4 3 2 1]);
            end
            temp(temp==fillvalue) = nan;
            temp(temp==0) = nan;


        end



        function [temp] = load_nc_2d(nc,Var_str,fillvalue,ind_avg)


            temp = nc{Var_str}(ind_avg,:,:);
            if ismatrix(temp)
                temp = temp';
            else
                temp = permute(temp,[3 2 1]);
            end
            temp(temp==fillvalue) = nan;
            temp(temp==0) = nan;
        end




        function temp_all = load_nc_2d_cat(filelist,Var_str,fillvalue,ind_avg_list)
            nfile = numel(filelist);
            if nfile ~= numel(ind_avg_list)
                error('filelist and ind_avg_list must have the same number of entries.');
            end

            if nfile == 0
                temp_all = [];
                return;
            end

            temp_chunks = cell(1,nfile);
            for ifile = 1:nfile
                nc = netcdf(filelist{ifile});
                ind_avg = ind_avg_list{ifile};
                temp = nc{Var_str}(ind_avg,:,:);
                if ismatrix(temp)
                    temp = temp';
                else
                    temp = permute(temp,[3 2 1]);
                end
                temp(temp==fillvalue) = nan;
                temp(temp==0) = nan;
                temp_chunks{ifile} = temp;
            end
            temp_all = cat(3,temp_chunks{:});
        end


        function temp_all = load_nc_3d_cat(filelist,Var_str,fillvalue,ind_avg_list)
            nfile = numel(filelist);
            if nfile ~= numel(ind_avg_list)
                error('filelist and ind_avg_list must have the same number of entries.');
            end

            if nfile == 0
                temp_all = [];
                return;
            end

            temp_chunks = cell(1,nfile);
            for ifile = 1:nfile
                nc = netcdf(filelist{ifile});
                ind_avg = ind_avg_list{ifile};
                temp = nc{Var_str}(ind_avg,:,:,:);
                if ismatrix(temp)
                    temp = permute(temp,[2 1]);
                elseif (ndims(temp)==3)
                    temp = permute(temp,[3 2 1]);
                else
                    temp = permute(temp,[4 3 2 1]);
                end
                temp(temp==fillvalue) = nan;
                temp(temp==0) = nan;
                temp_chunks{ifile} = temp;
            end
            temp_all = cat(4,temp_chunks{:});
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
            [temp_u, temp_v] = mor.load_nc_uv(nc,varname,fillvalue,ind_avg);
            temp_u_zint = squeeze(sum(mean(temp_u,4),3)); temp_v_zint = squeeze(sum(mean(temp_v,4),3));
            temp_u = mean(temp_u./Dz_u,4); temp_v = mean(temp_v./Dz_v,4);
            [temp_u_zlev,temp_v_zlev] = mor.interp2z_uv(temp_u,temp_v,z3D_u,z3D_v,mask_u,mask_v,zz);
            curl_var = mor.calc_curl(temp_u_zlev,temp_v_zlev,dx,dy);
            curl_var_zint = mor.calc_curl(temp_u_zint,temp_v_zint,dx,dy);
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
            yy = 0:dy:Ly;
            xx_f = 0:dx_f:Lx-dx_f; %%% N.B. don't inclue point at x=Lx to avoid double counting on periodic boundary condition
            yy_f = 0:dy_f:Ly-dy; %%% N.B. we ignore y in [Ly-dy,Ly] because it doesn't contribute to along-slope flows
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
                % pcolor(msk_dxc')
                % msk_dyc = repmat( xor(-etab_f(1:Nx_f-1,:)>dd(n),-etab_f(2:Nx_f,:)>dd(n)), [1 1 Nlay]);
                msk_dyc = repmat( and(xor(-etab_f(1:Nx_f-1,:)>dd(n),-etab_f(2:Nx_f,:)>dd(n)), ...
                    mask_f(1:Nx_f-1,:).*mask_f(2:Nx_f,:)~=0), [1 1 Nlay]);

                % pcolor(msk_dyc')

                % figure;pcolor(msk_dxc'); figure;pcolor(msk_dyc')
                data_dxc = 0.5*(data_f(:,1:Ny_f-1,:)+data_f(:,2:Ny_f,:)); %WHY hacks: remove periodicity
                data_dyc = 0.5*(data_f(1:Nx_f-1,:,:)+data_f(2:Nx_f,:,:));
                data_int(n,:) = sum(sum(dx_f.*msk_dxc.*data_dxc,1),2) + sum(sum(dy_f.*msk_dyc.*data_dyc,1),2);
            end



        end


        %%% convinient function when var is on rho point
        function   var_isoba = calc_iso_int_rho(var_rho,dx,dy,dd,coord_q,mask,interpFac,cntrlen)
            [~,~,var_q] = mor.rho2uvp_3d(var_rho);
            var_q(isnan(var_q))=0;
            var_isoba = mor.calc_iso_int(var_q,dx,dy,dd,coord_q,mask,interpFac) ./ cntrlen;
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
            Ny = size(data,2)-1;
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
                if any(~isnan(data_bin(msk_rho)))
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
            xx_f = 0:dx_f:Lx-dx_f; %%% WHY hacks: inclue point at x=Lx
            yy_f = 0:dy_f:Ly-dy; %%% WHY hacks: inclue point at y=Ly
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
            for n = 1:Nd
                msk = repmat(and(etab_f<-dd(n),mask_f~=0),[1 1 Nlay]);
                circ(n,:) = sum(sum(curl_f.*msk*dx_f*dy_f,1),2);
            end



            % circ = zeros(Nd,Nlay);


            % msk = repmat( mask_f==0, [1 1 Nlay]);
            % circ_mid = sum(sum(curl_f.*msk*dx_f*dy_f,1),2);
            %
            %
            % for n = 1:Nd
            %     msk = repmat( and(etab_f<-dd(n),mask_f~=0),[1 1 Nlay]);
            %     circ(n,:) = circ_mid+ sum(sum(curl_f.*msk*dx_f*dy_f,1),2);
            % end



            % for n = 1:Nd
            %     msk = repmat( and(etab_f<-dd(n),mask_f~=0),[1 1 Nlay]);
            %     circ(n,:) = circ_mid+ sum(sum(curl_f.*msk*dx_f*dy_f,1),2);
            % end


            for nn = 1:Nlay
                temp = circ(2:end,nn) - circ(1:end-1,nn);
                temp = find(temp ==0,1,'last')-1;
                circ(1:temp,nn)=0;
            end

        end

        %%% convinient function when var is on rho point




        function   circ_val = calc_iso_circ_pack(curl_var,dx,dy,dd,coord_q,mask,interpFac,cntrlen)
            Nz = size(curl_var,3);
            curl_var(isnan(curl_var)) = 0;
            curl_var(isinf(curl_var)) = 0;
            circ_val = mor.calc_iso_circ(curl_var,dx,dy,dd,coord_q,mask,interpFac) ./ repmat(cntrlen,[1 Nz]);
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
            [temp_u, temp_v] = mor.load_nc_uv(nc, varname , fillvalue, ind_avg);

            % Calculate the vertically integrated components
            temp_u_zint = squeeze(sum(mean(temp_u, 4), 3));
            temp_v_zint = squeeze(sum(mean(temp_v, 4), 3));

            % Calculate the mean u and v components divided by their respective cell thickness
            temp_u = mean(temp_u ./ Dz_u, 4);
            temp_v = mean(temp_v ./ Dz_v, 4);

            % Project u and v onto isobath directions
            var_alongisob = mor.u2rho(temp_u) .* ax + mor.v2rho(temp_v) .* ay;
            var_alongisob_zint = mor.u2rho(temp_u_zint) .* ax + mor.v2rho(temp_v_zint) .* ay;
        end

        function [var_alongisob, var_alongisob_zint,var_crossisob, var_crossisob_zint] = isocrossProj(nc, varname, fillvalue, ind_avg, Dz_u, Dz_v, ax, ay, cx, cy)
            % Load u and v data from netCDF file
            [temp_u, temp_v] = mor.load_nc_uv(nc, varname , fillvalue, ind_avg);

            % Calculate the vertically integrated components
            temp_u_zint = squeeze(sum(mean(temp_u, 4), 3));
            temp_v_zint = squeeze(sum(mean(temp_v, 4), 3));

            % Calculate the mean u and v components divided by their respective cell thickness
            temp_u = mean(temp_u ./ Dz_u, 4);
            temp_v = mean(temp_v ./ Dz_v, 4);

            % Project u and v onto along isobath directions
            var_alongisob = mor.u2rho(temp_u) .* ax + mor.v2rho(temp_v) .* ay;
            var_alongisob_zint = mor.u2rho(temp_u_zint) .* ax + mor.v2rho(temp_v_zint) .* ay;

            % Project u and v onto corss isobath directions
            var_crossisob = mor.u2rho(temp_u) .* cx + mor.v2rho(temp_v) .* cy;
            var_crossisob_zint = mor.u2rho(temp_u_zint) .* cx + mor.v2rho(temp_v_zint) .* cy;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%% do_bottom_avg
        %%%
        %%% Conduct avergaed within certain depth about the bottom
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function var_bottom_avg = do_bottom_avg2D(var, depthforavg, z2D, Dz, mask)
            % Function to compute bottom average of topostrophy up to a specified depth
            ffac = 10;
            % Get the dimensions of the input matrices
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
                    % Calculate the weighted average of topostrophy
                    var_bottom_avg(ii) = sum(var_interp(ind_zavg) .* DZf(ind_zavg)) ...
                        ./ sum(DZf(ind_zavg));
                end
            end
        end



        function var_bottom_avg = do_bottom_avg_4d(var, depthforavg, z3D, Dz, mask)
            % Function to compute bottom average of a certain quantity up to a specified depth

            % Get the dimensions of the input matrices
            [Nx, Ny, ~, Nt] = size(var);

            % Initialize the output matrix with NaNs
            var_bottom_avg = nan(Nx, Ny, Nt);

            % Loop through each grid point
            for ii = 1:Nx
                for jj = 1:Ny
                    % Check if the point is within the valid mask
                    if mask(ii,jj) ~= 0
                        % Find the indices where the depth is within the specified range
                        ind_zavg = find(z3D(ii,jj,:) <= z3D(ii,jj,1) - 0.5* Dz(ii,jj,1) + depthforavg);
                        % Calculate the weighted average of topostrophy
                        var_bottom_avg(ii,jj,:) = sum(var(ii,jj,ind_zavg,:) .* Dz(ii,jj,ind_zavg), 3) ...
                            ./ sum(Dz(ii,jj,ind_zavg), 3);
                    end
                end
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%% do_bottom_avg
        %%%
        %%% Conduct avergaed within certain depth about the bottom
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function var_bottom_avg = do_bottom_int2D(var, depthforavg, z2D, Dz, mask)
            % Function to compute bottom average of topostrophy up to a specified depth
            ffac = 10;
            % Get the dimensions of the input matrices
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
                    % Calculate the weighted average of topostrophy
                    var_bottom_avg(ii) = sum(var_interp(ind_zavg) .* DZf(ind_zavg));
                end
            end
        end

        function var_tavg = computeAverage(matName, variableName, t_start, t_end)
            % computeAverage calculates the mean of a specified variable
            % within a .mat file over a given range of the third dimension.
            %
            % Parameters:
            %   matName      - String specifying the .mat file name.
            %   variableName - String specifying the variable to load.
            %   t_start      - Starting index for the third dimension.
            %   t_end        - Ending index for the third dimension.
            %
            % Returns:
            %   usq_u_tavg - The averaged 2D matrix.

            matName = char(matName);
            variableName = char(variableName);

            variableInfo = whos('-file', matName, variableName);
            if isempty(variableInfo)
                error('Variable "%s" not found in %s.', variableName, matName);
            end

            % Validate the dimensions of the variable
            if numel(variableInfo.size) < 3
                error('Variable "%s" does not have at least three dimensions.', variableName);
            end

            % Ensure the time indices are within bounds
            sz = variableInfo.size(3);
            if t_start < 1 || t_end > sz || t_start > t_end
                error('Invalid time indices: t_start=%d, t_end=%d, but size in 3rd dimension is %d.', ...
                    t_start, t_end, sz);
            end

            % Use MATFILE when possible to avoid loading the full variable into memory.
            try
                dataObj = matfile(matName);
                dataSlice = dataObj.(variableName)(:, :, t_start:t_end);
            catch
                data = load(matName, variableName);
                dataSlice = data.(variableName)(:, :, t_start:t_end);
            end

            % Compute the mean over the specified range in the third dimension
            var_tavg = mean(dataSlice, 3);
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
            rho_w = mor.rho2w2d(rho,z2D_r,z2D_w);
            rho_dz= (rho_w(:,2:end) - rho_w(:,1:end-1))./(z2D_w(:,2:end) - z2D_w(:,1:end-1));
            %%% compute (drho/dx)_sigma
            rho_dx = zeros(Nx, Nz);
            rho_dx(3:Nx-2, :) = (rho(4:Nx-1, :) - rho(2:Nx-3, :)) / (2 * dx);
            rho_dx(2, :) = (rho(3, :) - rho(2, :)) / dx;
            rho_dx(Nx-1, :) = (rho(Nx-1, :) - rho(Nx-2, :)) / dx;
            %%% compute (dz/dx)_sigma
            S_slope = zeros(Nx, Nz);
            S_slope(3:Nx-2, :) = (z2D_r(4:Nx-1, :) - z2D_r(2:Nx-3, :)) / (2 * dx);
            S_slope(2, :) = (z2D_r(3, :) - z2D_r(2, :)) / dx;
            S_slope(Nx-1, :) = (z2D_r(Nx-1, :) - z2D_r(Nx-2, :)) / dx;
            %%% compute (d/dx)_z
            %%% (d/dx)_z = (d/dx)_sigma - [(dz/dx)_sigma]* [(d/dz)]
            rho_dx_sigma_extra= -rho_dz.*S_slope;
            rho_dxi=rho_dx+rho_dx_sigma_extra;
        end


        function [rho_dxi,rho_dz]=diffxi3d(rho,dx,z3D_r,z3D_w)
            [Nx, Ny, Nz] = size(rho);
            %%% compute vertical gradient
            rho_w = mor.rho2w(rho,z3D_r,z3D_w);
            rho_dz= (rho_w(:,:,2:end) - rho_w(:,:,1:end-1))./(z3D_w(:,:,2:end) - z3D_w(:,:,1:end-1));
            %%% compute (drho/dx)_sigma
            rho_dx = zeros(Nx, Ny, Nz);
            rho_dx(3:Nx-2, :, :) = (rho(4:Nx-1, :, :) - rho(2:Nx-3, :, :)) / (2 * dx);
            rho_dx(2, :, :) = (rho(3, :, :) - rho(2, :, :)) / dx;
            rho_dx(Nx-1, :, :) = (rho(Nx-1, :, :) - rho(Nx-2, :, :)) / dx;
            %%% compute (dz/dx)_sigma
            S_slope = zeros(Nx, Ny,Nz);
            S_slope(3:Nx-2, :, :) = (z3D_r(4:Nx-1, :, :) - z3D_r(2:Nx-3, :, :)) / (2 * dx);
            S_slope(2, :, :) = (z3D_r(3, :, :) - z3D_r(2, :, :)) / dx;
            S_slope(Nx-1, :, :) = (z3D_r(Nx-1, :, :) - z3D_r(Nx-2, :, :)) / dx;
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        function [rho_deta,rho_dz] = diffeta(rho,dy,z2D_r,z2D_w)
            [Ny, Nz] = size(rho);
            %%% compute vertical gradient
            rho_w = mor.rho2w2d(rho,z2D_r,z2D_w);
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
            rho_w = mor.rho2w(rho,z3D_r,z3D_w);
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


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%% vort3D(u,v,dx,dy,z3D_r,z3D_w)
        %%%
        %%% Compute relative vorticity on a 3D grid
        %%%
        %%% Inputs:
        %%%   u,v     : velocity components (Nx,Ny,Nz)
        %%%   dx,dy   : grid spacing in x and y
        %%%   z3D_r   : depth at rho points (Nx,Ny,Nz)
        %%%   z3D_w   : depth at w points   (Nx,Ny,Nz+1)
        %%%
        %%% Output:
        %%%   vort    : relative vorticity (Nx,Ny,Nz)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function vort = vort3D(u,v,dx,dy,z3D_r,z3D_w)

            [Nx,Ny,Nz] = size(u);

            %---------------------------------------------------------
            % dv/dx at constant z
            % Loop over j (y-index) to apply diffxi in x-z planes
            dv_dx = zeros(Nx,Ny,Nz);
            for j = 1:Ny
                [dv_dx(:,j,:), ~] = mor.diffxi(squeeze(v(:,j,:)), dx, ...
                    squeeze(z3D_r(:,j,:)), ...
                    squeeze(z3D_w(:,j,:)));
            end

            %---------------------------------------------------------
            % du/dy at constant z
            % Loop over i (x-index) to apply diffeta in y-z planes
            du_dy = zeros(Nx,Ny,Nz);
            for i = 1:Nx
                [du_dy(i,:,:), ~] = mor.diffeta(squeeze(u(i,:,:)), dy, ...
                    squeeze(z3D_r(i,:,:)), ...
                    squeeze(z3D_w(i,:,:)));
            end

            %---------------------------------------------------------
            % Relative vorticity
            vort = dv_dx - du_dy;

        end






        function  Var_cravg = CRavg(var,Nx)

            if(length(var)~=Nx)

                Var_cravg = 0.5*(var(Nx/2-1:-1:1,:)+var(Nx/2:end,:));
            else


                if(size(var)==1)
                    Var_cravg = 0.5*(var(Nx/2:-1:2)+var(Nx/2+1:end-1));
                else
                    Var_cravg = 0.5*(var(Nx/2:-1:2,:)+var(Nx/2+1:end-1,:));
                end

            end
        end






    end
end
