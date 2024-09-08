function [absorption_no, reflect_no, trans_no] = monte_carlo_double_layer(h1, h2, wl, scat_prob1, scat_prob2, ext_tot1, ext_tot2, g1, g2, QY_modified, start_wl, number_wl, inv_cdf, cos_teta_prime, sur_reflection, n_medium, k_medium, n_subs, k_subs)
    absorption_no = 0;
    reflect_no = zeros(number_wl, 1);
    trans_no = zeros(number_wl, 1);
    wl_index = wl - start_wl + 1;
    reflect_no(wl_index) = sur_reflection;
    energy = 1 - sur_reflection;
   
    % Initialize position and direction
    x = 0; y = 0; z = 0;
    phi = 2 * pi * rand;
    cos_phi = cos(phi);
    sin_phi = sqrt(1 - cos_phi^2);
    sin_teta_prime = sqrt(1 - cos_teta_prime^2);
    s_x = sin_teta_prime * sin_phi;
    s_y = sin_teta_prime * cos_phi;
    s_z = cos_teta_prime;
    alive = 1;
    l_beta = -log(rand()) / ext_tot1(wl_index); % extinction length in top layer
    current_layer = 1; % Start in the top layer
   
    while alive
        if current_layer == 1
            scat_prob = scat_prob1;
            ext_tot = ext_tot1;
            g = g1;
        else
            scat_prob = scat_prob2;
            ext_tot = ext_tot2;
            g = g2;
        end
       
        % Calculate distance to the next boundary
        if current_layer == 1 % Top layer
            if s_z > 0 % Photon moving towards bottom boundary
                l_w = (h1 - z) / s_z;
            else % Photon moving towards top boundary
                l_w = -z / s_z;
            end
        else % Bottom layer
            if s_z > 0 % Photon moving towards bottom boundary of the second layer
                l_w = (h1 + h2 - z) / s_z;
            else % Photon moving towards top boundary of the second layer
                l_w = (h1 - z) / s_z;
            end
        end
       
        % Determine whether the photon will hit a boundary or be scattered/absorbed first
        if l_w < l_beta
            min_index = 1;
            min_l = l_w;
        else
            min_index = 2;
            min_l = l_beta;
        end
       
        % Update the photon's position
        x = x + min_l * s_x;
        y = y + min_l * s_y;
        z = z + min_l * s_z;
       
     if min_index == 1 % Photon hits a boundary
        if current_layer == 1 % Top layer boundary interaction
            if s_z > 0 % Photon moving down
                alive2 =snell(s_z, n_medium(wl_index), k_medium(wl_index), n_medium(wl_index), k_medium(wl_index));
                    if alive2 ==0  %if medium2 is same as medium1, alive will always be zero. Including it for the general case when medium2 neq to medium1
                        current_layer = 2;
                        z = h1; % Set position at the top of the bottom layer
                        l_beta = -log(rand()) / ext_tot2(wl_index); % Recalculate extinction length for bottom layer
                    else
                        s_z = -s_z; % Reflect back within the top layer
                        l_beta = l_beta - l_w;
                    end
           else % Photon moving up
                alive = snell(s_z, n_medium(wl_index), k_medium(wl_index), n_subs(wl_index), k_subs(wl_index));
                     if alive==0  
                     reflect_no(wl_index)=reflect_no(wl_index)+energy;
                     else
                     s_z = -s_z; % Reflect back within the top layer    
                     l_beta = l_beta - l_w;
                    end
            end
       else %current layer is the second layer
          
            if s_z > 0 % Photon moving down
       
                alive = snell(s_z, n_medium(wl_index), k_medium(wl_index), n_subs(wl_index), k_subs(wl_index));
                    if alive ==0  
                       %disp('photon is in layer 2 and will transmit');
                        trans_no(wl_index)=trans_no(wl_index)+energy;
                    else
                        %disp('photon is in layer 2 and will reflect from bottom layer');
                        s_z = -s_z; % Reflect back within the second layer
                        l_beta = l_beta - l_w;
                    end
             else % Photon moving up
                alive3 = snell(s_z, n_medium(wl_index), k_medium(wl_index), n_medium(wl_index), k_medium(wl_index));
                     if alive3==0 % in this case alive will always be zero since the medium is same on both sides
                        current_layer = 1;
                        z = h1; % Set position at the top of the bottom layer
                        l_beta = -log(rand()) / ext_tot1(wl_index); % Recalculate extinction length for bottom layer
                    else
                        s_z = -s_z; % Reflect back within the second layer
                        l_beta = l_beta - l_w;
                    end
            end
         end
        else % Photon interacts within the layer (scattering or absorption)
            if rand() < scat_prob(wl_index)
                % Scattering event
                [s_x, s_y, s_z] = scatter(g(wl_index), s_x, s_y, s_z);
                l_beta = -log(rand()) / ext_tot(wl_index);
%             elseif rand() < scat_prob2(wl_index)
%                 [s_x, s_y, s_z] = scatter(g2(wl_index), s_x, s_y, s_z);
%                 l_beta = -log(rand()) / ext_tot2(wl_index);
            else
                % Absorption or fluorescence emission event
                if current_layer == 1 && rand() > QY_modified(wl_index)
                    % Absorption in the top layer without re-emission
                    alive = 0;
                    absorption_no = absorption_no + energy;
                elseif current_layer == 1
                    % Fluorescence emission in the top layer
                    angles = random_angles();
                    s_x = angles(1);
                    s_y = angles(2);
                    s_z = angles(3);
                    wl_new = round(inv_cdf(ceil(rand() * 10000)));
                    absorption_no = absorption_no + energy * (1 - wl / wl_new);
                    energy = energy * wl / wl_new;
                    wl = wl_new;
                    wl_index = wl - start_wl + 1;
                    l_beta = -log(rand()) / ext_tot(wl_index);
                    % Ensure fluorescence emission is scattered by the bottom layer
                    %current_layer = 2;
                else
                    % Absorption in the bottom layer without fluorescence
                    alive = 0;
                    absorption_no = absorption_no + energy;
                end
            end
      end
    end
end

