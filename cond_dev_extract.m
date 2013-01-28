%------------------------------------------------------------------------------------------------
% Author: Hassan Rezaee, Email: h.rezaee@ut.ac.ir
% Paper: "Multiple-Point Geostatistical Simulation Using the Bunch-Pasting Direct Sampling Method"
%------------------------------------------------------------------------------------------------

% This function extracts the conditioning data from the simulation grid
% which contains both primary conditioning and previously-simulated data.

%---------------------------------------------------------------------------------------------
% Input parameters:
%   window_type: search window type within which the data are extracted

%   point: this is the point under simulation. The window of pre-defined
%          shape is located over this point

%   conditioning_data: self-explained!

%   R: search dimensions. in based on the window type, these
%                dimensions may refer to different parameter (see Bunch_DS.m)

%   sim_dim: simulation grid dimensions along X, Y and Z. The algorithm
%            calls this automatically
%---------------------------------------------------------------------------------------------

% Outputs:
%   data_event_cond: this is in the form [x_coord y_coord z_coord variable]

%   informed_nodes_id: to better be addressed the IDs of all conditioning 
%                      data extracted from the simulation grid and used for
%                      further processing

function [data_event_cond informed_nodes_id] = cond_dev_extract(window_type , point , conditioning_data , R , sim_dim)

sim_x1 = point(1,1);
sim_y1 = point(1,2);
sim_z1 = point(1,3);

a=length(R);

switch a
    case 3
        search_radius_x = R(1);
        search_radius_y = R(2);
        search_radius_z = R(3);
    case 2
        outer_search_radius = R(1);
        inner_search_radius = R(2);
    case 1
        search_radius = R;
end

    switch window_type
        case 1   % 1.  Extracting conditioning data within a CIRCLE            
            dist=sqrt((sim_x1-conditioning_data(:,1)).^2+(sim_y1-conditioning_data(:,2)).^2+(sim_z1-conditioning_data(:,3)).^2);   
            informed_nodes_id = find(and(dist <= floor(search_radius/2) , dist > 0));
            data_event_cond = conditioning_data(informed_nodes_id,1:5);

        case 2   % 2.  Extracting conditioning data within a SQUARE           
            if sim_dim > 1
                informed_nodes_id = find(and(and(and(conditioning_data(:,1)>(sim_x1-floor(search_radius_x/2)),conditioning_data(:,1)<(sim_x1+floor(search_radius_x/2))),...
                and(conditioning_data(:,2)>(sim_y1-floor(search_radius_y/2)),conditioning_data(:,2)<(sim_y1+floor(search_radius_y/2)))),...
                and(conditioning_data(:,3)>(sim_z1-floor(search_radius_z/2)),conditioning_data(:,3)<(sim_z1+floor(search_radius_z/2)))));
            elseif sim_dim == 1
                informed_nodes_id = find(and(and(conditioning_data(:,1)>(sim_x1-floor(search_radius_x/2)),conditioning_data(:,1)<(sim_x1+floor(search_radius_x/2))),...
                and(conditioning_data(:,2)>(sim_y1-floor(search_radius_y/2)),conditioning_data(:,2)<(sim_y1+floor(search_radius_y/2)))));
            end
                data_event_cond = conditioning_data(informed_nodes_id,1:5);
        
        case 3   % 3.  Extracting from a DISC-SHAPED area
            dist_disc = sqrt((sim_x1-conditioning_data(:,1)).^2+(sim_y1-conditioning_data(:,2)).^2+(sim_z1-conditioning_data(:,3)).^2);
            informed_nodes_id = find(and(dist_disc >= inner_search_radius , dist_disc <= outer_search_radius));
            data_event_cond = conditioning_data(informed_nodes_id,:);
    end
    