%------------------------------------------------------------------------------------------------
% Author: Hassan Rezaee, Email: h.rezaee@ut.ac.ir
% Paper: "Multiple-Point Geostatistical Simulation Using the Bunch-Pasting Direct Sampling Method"
%------------------------------------------------------------------------------------------------

% This is the main function of Bunch_DS simulation algorith. Both 
% Conditional/Non-conditional Simulation of Categorical/Continuous
% Variables Using Bunch-Pasting/Pixel-Based DS Method can be done by this 
% code. The function is run from another Matlab script named
% running_the_main_function.m. All the input parameters are explained 
% below and also in the User Guide provided in the code package.

% USE: [Final_Realizations post_maps] = Bunch_DS(ti_name,Condition_mode,conditioning_weight,cond_data_name,sim_dim,Search_Radii,fract_of_ti_to_scan,...
%                            dist_threshold,max_no_nodes_in_data_event,window_type,bunch_size,sim_path_mode,ti_path_mode,...
%                            No_Realizations,post_proc_mode,simulation_show_mode);

%------------------------------------------------------------------------------------------------
% INPUT:
%   ti_name                     : the name of the ti. ti should be saved in a txt file in the directory
%                                   please read the manual on how to prepare the ti data file

%   Condition_mode              : consider the following codes to do:
%                                   0. non-Conditional Simulation
%                                   1. conditional Simulation.

%   conditioning_weight         : the weight of conditioning data when calculating the distance between Nx and Ny

%   cond_data                   : If the conditional simulation is selected then the 
%                                   you should enter the name of the text
%                                   file of the conditioning data. This file should be located in the directory of 
%                                   the matlab code. This data file should be in the form [x_coord y_coord z_coord variable]
%                                   please read the manual on how to prepare the cond data file

%   sim_dim                     : number of simulation grid nodes along x, y and z directions
%                                   for example [10 10 1] is 2D

%   Search_Radii                : this should be in the form of [X Y Z] to determine the search radius in each direction.
%                                   Please read the manual on how these values are assinged to the 
%                                   search radii in different window types.

%   fract_of_ti_to_scan         : the fraction of ti to be scanned; a number is to be chosen between [0-1]
%                                   for example 0.1 or 0.5

%   dist_threshold              : distance threshold between dev_ti and dev_cond, default = 0.05

%   max_no_nodes_in_data_event  : maximum number of nodes in data event
%                                   which is formed by both conditioning and previously-simulated data
%                                   a number of 30-50 would deliver realistic results
%                                   higher numbers exponentially decreases the computational efficiency of the algorithm
                                 
%   window_type                 : consider the following codes to specify the search window type
%                                   1. circular
%                                   2. squared
%                                   3. disc Shape

%   bunch_size                  : if bunch_size > 0 the simulation mode is through bunch pasting manner
%                                   otherwise (bunch_size=0), pixel-based (traditional) DS is done
%                                   for example if "bunch_size = 1", then a number of (2*bunch_size+1)^2 = 9 nodes 
%                                   in 2D and (2*bunch_size+1)^3 = 27 nodes in 3D  are to be pasted
%                                   If you want to run the original DS (Martiethoz, 2010) consider bunch_size = 0;

%   sim_path_mode               : use the following code to specify the simulation path
%                                   1. random
%                                   2. unilateral

%   ti_path_mode                : use the following codes to determine the way through which ti is scanned:
%                                   1. random
%                                   3. unilateral with starting point

%   No_Realizations             : number of realizations to be produced

%   post_proc_mode              : if No_Realizations > 1 then use the following codes to carry out post-processing steps:
%                                   0. no post-processing step is taken after the simulation
%                                   1. E-type, Inter Quartile Range (IQR) and Quantile (0.5) maps are calculated after the simulation

%   simulation_show_mode        : to illustrate the realization while simulation consider 1 otherwise any other code
%                                 use this in case you want to check  the simulation procedure other do not use it
%                                 since it overloads the code.

%-----------------------------------------------------------------------------------------------------------
% OUPUT: 
%   Final_Realizations          : final realizations produced by this method in the form 
%                                 [x_coord y_coord z_coord realz#1 realz#2 ,... realz#n], n = No_Realizations
%                                 This is the input format for SGeMS (Remy et al., 2008).

%   post_maps                   : post-processed maps of E-type, IQR and Quantile values. 

%   parfile                     : for every simulation you do the input parameters are saved in a text file 
%                                 under the name of parfile.txt. You may change the name of this file each 
%                                 time at the last line of this script.

%-----------------------------------------------------------------------------------------------------------
% This program uses the following subroutines:
%   dlmcell.m                    : write cell array to text file (Roland Pfister, Version: 01.06.2010)
%   cond_dev_extract.m           : extracts the data event from the neighborhood of the node under simulation
%   sim_post_proc.m              : calculates the post-processing maps of E-Type, Interquartile and Quantile maps

%-----------------------------------------------------------------------------------------------------------
% An Example of Input Parameters for Simulation is provided in
% running_the_main_function.m matlab script.

%-----------------------------------------------------------------------------------------------------------
function [Final_Realizations post_maps] = Bunch_DS(ti_name,Condition_mode,conditioning_weight,cond_data_name,sim_dim,Search_Radii,fract_of_ti_to_scan,...
                           dist_threshold,max_no_nodes_in_data_event,window_type,bunch_size,sim_path_mode,ti_path_mode,...
                           No_Realizations,post_proc_mode,simulation_show_mode)

fid = fopen([ti_name '.txt']);
ti_grid = textscan(fid, '%f %f %f %f');
fclose(fid);
ti_grid=cell2mat(ti_grid);

ti_size(1:3) = (max(ti_grid(:,1:3))-min(ti_grid(:,1:3)))+1;

if and(ti_size(3)==1,sim_dim(3)>1)
       errordlg('You are not allowed to use a 2D TI for 3D simulation, please reduce sim_dim(3) to 1');
       pause(3)
       error('You are not allowed to use a 2D TI for 3D simulation, please reduce sim_dim(3) to 1');
end

if and(post_proc_mode == 1 , No_Realizations==1)
   errordlg('The number of realizations (input name: "No_Realizations" ) should be more than 1 to produce post-processing maps');
   pause(3)
   error('The number of realizations (input name: "No_Realizations" ) should be more than 1 to produce post-processing maps');
end

%% Search Radii Specifics
search_radius_x = Search_Radii(1);
search_radius_y = Search_Radii(2);
search_radius_z = Search_Radii(3);

switch window_type
    case 1
        if sim_dim(3)<=1
            circ_search_radius = min([search_radius_x ; search_radius_y]); 
        elseif sim_dim(3)>1
            circ_search_radius = min([search_radius_x ; search_radius_y ; search_radius_z]); 
        end
        search_radii = circ_search_radius;

    case 2
            sqr_search_radius_x = search_radius_x;
            sqr_search_radius_y = search_radius_y;
            sqr_search_radius_z = search_radius_z;
            search_radii = [sqr_search_radius_x sqr_search_radius_y sqr_search_radius_z];

    case 3
        if sim_dim(3)<=1
            disc_inner_search_radius = min([search_radius_x ; search_radius_y]); 
            disc_outer_search_radius = max([search_radius_x ; search_radius_y]);
        elseif sim_dim(3)>1
            disc_inner_search_radius = min([search_radius_x ; search_radius_y ; search_radius_z]); 
            disc_outer_search_radius = max([search_radius_x ; search_radius_y ; search_radius_z]);
        end
        if disc_outer_search_radius == disc_inner_search_radius;
            disc_inner_search_radius = 0.5*disc_outer_search_radius;
        end
        search_radii = [disc_outer_search_radius ; disc_inner_search_radius];
end

% check for the type of variable; it is supposed that the number of facies in categorical TI is less than 10!
[B]=sort(ti_grid(:,4),1);
No_facies=0;
for i=2:length(B)
    if (B(i-1)~=B(i))
        No_facies = No_facies+1;
    end
end

%% Preparing the TI
% the grid of the training image [x_coord y_coord z_coord variable] should be introduced to the algorithm
ti=reshape(ti_grid(:,4),ti_size(1),ti_size(2),ti_size(3));
if ti_size(3)==1
    imagesc(ti);axis xy;axis equal;xlim([1 ti_size(2)]);ylim([1 ti_size(1)]);colorbar;title('Training Image')
elseif ti_size(3)>1
    imagesc(ti(:,:,fix(mean(ti_grid(:,3)))));axis xy;axis equal;xlim([1 ti_size(2)]);ylim([1 ti_size(1)]);
    colorbar;title('Training Image')
end
if No_facies<=2;colormap('gray');end

% hw = waitbar(0,'Please wait...');

if No_facies>10
    prompt={'Is this your TI? Is it Continuous? (Y/N)'};
elseif No_facies<=10
    prompt={['Is this your TI (a section of it if 3D) ? Is it Categorical with: ' num2str(No_facies) ' facies ? (Y/N)']};
end

name='TI Check';
numlines=1;
defaultanswer={'Y','Cat'};
answer=inputdlg(prompt,name,numlines,defaultanswer);

if strcmp(answer(1),{'Y'}) == 1
    close all
    msgbox('OK! The simulation will start in few seconds');
    pause(3)
elseif strcmp(answer,{'N'}) == 1
    msgbox(['Sorry! First, please read the manual on how to prepare the training'...
           ' image format. It should be noted that the No of facies within'...
           ' a given TI is 10 to consider it as a categorical one, if you'...
           ' wish to consider the categorical ti with more than 10 facies, please advise here']);
    pause(5)
    prompt={'Categorical with more than 10 facies'};
    name='Do you wish to consider yout ti as a categorical one still ? (Y/N)';
    numlines=1;
    defaultanswer={'Y'};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    
    if strcmp(answer(1),{'Y'}) == 1
        Distance_type = 1; %#ok<NASGU>
        msgbox('OK! The simulation will start in few seconds');
        close all
    elseif strcmp(answer(1),{'N'}) == 1
        msgbox(['We guess the problem is due to the original ti'...
                            ' file. Please check the coordinates of your ti.' ...
                            ' Please read the manual on how to prepare the ti file.']);
        close all
        error(['We guess the problem is due to the original ti'...
                            ' file. Please check the coordinates of your ti.' ...
                            ' Please read the manual on how to prepare the ti file.']);
    end
end

allHandle = allchild(0);
allTag = get(allHandle, 'Tag');
isMsgbox = strncmp(allTag, 'Msgbox_', 7);
delete(allHandle(isMsgbox));


if No_facies > 10
    disp('Continuous Variable is selected')
    pause(2)
    Distance_type=2;
else
    disp('Categorical Variable is selected')
    pause(2)
    Distance_type=1;
end

% comparing the search radii with ti dimensions.
if search_radius_x > ti_size(1)
    disp('The search radius cannot exceed the ti size in that direction')
    disp(['In your case search_radius_x is changed from ' num2str(search_radius_x) ' to ' num2str(ti_size(1)+1)])
    search_radius_x=ti_size(1)+1;
    pause(0.3)
end
if search_radius_y>ti_size(2)
    disp('The search radius cannot exceed the ti size in that direction')
    disp(['In your case search_radius_y is changed from ' num2str(search_radius_y) ' to ' num2str(ti_size(2)+1)])
    search_radius_y=ti_size(2)+1;
    pause(0.3)
end
if search_radius_z>ti_size(3)
    disp('The search radius cannot exceed the ti size in that direction')
    disp(['In your case search_radius_z is changed from ' num2str(search_radius_z) ' to ' num2str(ti_size(3)+1)])
    search_radius_z=ti_size(3)+1;
    pause(0.3)
end

% comparing the search radii with simulation grid dimensions.
if search_radius_x > sim_dim(1)
    disp('The search radius cannot exceed the simulation grid in that direction')
    disp(['In your case search_radius_x is changed from ' num2str(search_radius_x) ' to ' num2str(floor(sim_dim(1)/2)+1)])
    search_radius_x=floor(sim_dim(1)/2)+1;
    pause(0.3)
end
if search_radius_y > sim_dim(2)
    disp('The search radius cannot exceed the simulation grid in that direction')
    disp(['In your case search_radius_y is changed from ' num2str(search_radius_y) ' to ' num2str(floor(sim_dim(2)/2)+1)])
    search_radius_y=floor(sim_dim(2)/2)+1;
    pause(0.3)
end
if search_radius_z > sim_dim(3)
    disp('The search radius cannot exceed the simulation grid in that direction')
    disp(['In your case search_radius_z is changed from ' num2str(search_radius_z) ' to ' num2str(floor(sim_dim(3)/2)+1)])
    search_radius_z=floor(sim_dim(3)/2)+1;
    pause(0.3)
end
if length(find((sim_dim)==1))>1
   search_radius_x=ceil(sim_dim(1)/2)+1; 
   search_radius_y=ceil(sim_dim(2)/2)+1; 
   search_radius_z=ceil(sim_dim(3)/2)+1;
   if simulation_show_mode == 1
       msgbox('Simulation grid dimensions are too limited to illustrate it while doing the simulation');
       pause(3)
       allHandle = allchild(0);
       allTag = get(allHandle, 'Tag');
       isMsgbox = strncmp(allTag, 'Msgbox_', 7);
       delete(allHandle(isMsgbox));
   end
end
if and(sim_dim(3)>1,search_radius_z<=1)
    search_radius_z=2;
end
%% creating empty simulation grid

% Creating an outer shell around Simulation Grid
[y z x]=meshgrid(1:sim_dim(2)+2*floor(search_radius_y/2) ,...
                 1:sim_dim(3)+2*floor(search_radius_z/2) ,...
                 1:sim_dim(1)+2*floor(search_radius_x/2));
             
xx=reshape(x,[],1);
yy=reshape(y,[],1);
zz=reshape(z,[],1);

sim_grid(:,1)=xx;
sim_grid(:,2)=yy;
sim_grid(:,3)=zz;
sim_grid(:,4)=min(ti_grid(:,4))-2;

x_active = and(sim_grid(:,1)>=((floor(search_radius_x/2))) , sim_grid(:,1) < max(sim_grid(:,1))-(floor(search_radius_x/2)));
y_active = and(sim_grid(:,2)>=((floor(search_radius_y/2))),sim_grid(:,2) < max(sim_grid(:,2))-(floor(search_radius_y/2)));
z_active = and(sim_grid(:,3)>=((floor(search_radius_z/2))),sim_grid(:,3) < max(sim_grid(:,3))-(floor(search_radius_z/2)));

if sim_dim(3) == 1
    z_active=ones(size(x_active,1),1);
elseif sim_dim(2) == 1
    y_active=ones(size(y_active,1),1);
elseif sim_dim(1) == 1
    x_active=ones(size(x_active,1),1);
end

id=and(and(x_active==1 , y_active==1) , z_active==1);
sim_grid(id==1,4)=nan;
SIM_GRID = sim_grid;
% extra variables to show the contribution of different states that a node can be simulated
c_empty_dv=0;
c_count_thre=0;
c_distance=0;
T_one=zeros(No_Realizations,1);

%% Realizations over sim_grid are to be produced
Realization=zeros(size(sim_grid,1),No_Realizations);
for real = 1:No_Realizations
    
    active = sim_grid(and(and((x_active==1) , y_active==1) , z_active==1), : );
    active(:,4) = nan;
    active(:,5)=0;   
    diff=sum(isnan(sim_grid(:,4)))-sim_dim(1)*sim_dim(2)*sim_dim(3);

% Conditional simulation
if Condition_mode == 1
    
    % First to read the conditioning data file
    if real == 1       
       fid = fopen([cond_data_name '.txt']);
       cond_data = textscan(fid, '%f %f %f %f');
       fclose(fid);
       cond_data=cell2mat(cond_data);      
       Hard=cond_data;   
    elseif real > 1
       cond_data=Hard;
    end
    
    % Due to the outer shell around sim_grid, hard data are re-located.
    cond_data(:,1) = cond_data(:,1)+floor(search_radius_x/2);
    cond_data(:,2) = cond_data(:,2)+floor(search_radius_y/2);
    cond_data(:,3) = cond_data(:,3)+floor(search_radius_z/2);
    
    cond_data(:,5) = conditioning_weight;
    
    % just the "cond data" within the simulation grid should be identified
    a = and(and(cond_data(:,1)<sim_dim(1)+floor(search_radius_x/2),...
                cond_data(:,2)<sim_dim(2)+floor(search_radius_y/2)),...
                cond_data(:,3)<=sim_dim(3)+floor(search_radius_z/2));

    cond_data = cond_data(a,:);
    Hard_within=cond_data;
    sim_grid(:,5)=1;
    
    % Inserting cond data in simulation grid
    for i=1:size(cond_data,1)
        a = and(and(sim_grid(:,1)==cond_data(i,1),sim_grid(:,2)==cond_data(i,2)),sim_grid(:,3)==cond_data(i,3));
        sim_grid(a,4)=cond_data(i,4);
        sim_grid(a,5)=conditioning_weight;
    end
    
    % Nodes are visited along a pre-determined path: random or unilateral
    if sim_path_mode == 2
        sim_path = 1:size(active,1);
    elseif sim_path_mode == 1
        sim_path = randperm(size(active,1));
    end


    % Eliminating hard data point from the simulation path
    for i=1:length(sim_path)
        a=find(and(and(active(i,1)==cond_data(:,1),active(i,2)==cond_data(:,2)),active(i,3)==cond_data(:,3)), 1);       
        if ~isempty(a);sim_path(i)=-1;end        
    end
    
    sim_path(sim_path == -1) = [];
    hard_data = cond_data;
    
% Non-conditional simulation  
elseif Condition_mode==0    
    cond_data = [-9999 -9999 -9999 0 1];
    sim_grid(:,5)=1;
    if sim_path_mode == 2
        sim_path = 1:size(active,1);
    elseif sim_path_mode == 1
        sim_path = randperm(size(active,1));
    end
end

tic

while 1==1

% for every simulation node do
for i = 1:length(sim_path)

    % The coordinates of the nodes under simulation
    sim_x = active(sim_path(i),1);
    sim_y = active(sim_path(i),2);
    sim_z = active(sim_path(i),3);
    
    sim_node_id = and(and(sim_grid(:,1)==sim_x , sim_grid(:,2)==sim_y) , sim_grid(:,3)==sim_z);
  
if isnan(sim_grid(sim_node_id==1,4)) == 1
    
    % To show the progress of simulation
    simulated_percentage=(100-((sum(isnan(sim_grid(:,4)))-diff)/(sim_dim(1)*sim_dim(2)*sim_dim(3)))*100);

    disp(['Realization No.  ' num2str(real) '  is being produced   ' num2str(simulated_percentage), '  % completed']);
%     waitbar(simulated_percentage/100,hw,sprintf('%12.9f',fix(simulated_percentage)))
       % To extract the bunch-shape space to which the bunch from TI is to be pasted  
if sim_dim(3)>1
    
    sim_bunch_counter = 0;
    points_sim = zeros((2*bunch_size+1)^3,3);
    bunch_ids_sim_coord = zeros((2*bunch_size+1)^3,3);

    for x_id_sim = -bunch_size:bunch_size
        for y_id_sim = -bunch_size:bunch_size
            for z_id_sim = -bunch_size:bunch_size
                p_x=sim_x + x_id_sim; p_y=sim_y + y_id_sim; p_z=sim_z + z_id_sim;
                
            if and(isnan(sim_grid(and(and(sim_grid(:,1) == p_x , sim_grid(:,2) == p_y) , sim_grid(:,3) == p_z),4)) == 1,...
                   and(and(and(p_x <= max(active(:,1)) , p_x >= (min(active(:,1)))),...
                           and(p_y <= max(active(:,2)) , p_y >= (min(active(:,2))))),...
                           and(p_z <= max(active(:,3)) , p_z >= (min(active(:,3))))))
                       
                sim_bunch_counter = sim_bunch_counter+1;
                points_sim(sim_bunch_counter,1) = p_x;
                points_sim(sim_bunch_counter,2) = p_y;
                points_sim(sim_bunch_counter,3) = p_z;
                bunch_ids_sim_coord(sim_bunch_counter,1) = x_id_sim;
                bunch_ids_sim_coord(sim_bunch_counter,2) = y_id_sim;
                bunch_ids_sim_coord(sim_bunch_counter,3) = z_id_sim;
            end     
            end
        end
    end
    
elseif sim_dim(3)==1
    
    sim_bunch_counter = 0;
    points_sim = zeros((2*bunch_size+1)^2,3);
    bunch_ids_sim_coord = zeros((2*bunch_size+1)^2,3);

    for x_id_sim = -bunch_size:bunch_size
        for y_id_sim = -bunch_size:bunch_size
                p_x=sim_x + x_id_sim; p_y=sim_y + y_id_sim;
                
            if and(isnan(sim_grid(and(sim_grid(:,1) == p_x , sim_grid(:,2) == p_y) ,4)) == 1,...
                   and(and(p_x <= max(active(:,1)) , p_x >= (min(active(:,1)))),...
                       and(p_y <= max(active(:,2)) , p_y >= (min(active(:,2))))))
                       
                sim_bunch_counter = sim_bunch_counter+1;
                points_sim(sim_bunch_counter,1) = p_x;
                points_sim(sim_bunch_counter,2) = p_y;
                bunch_ids_sim_coord(sim_bunch_counter,1) = x_id_sim;
                bunch_ids_sim_coord(sim_bunch_counter,2) = y_id_sim;
            end
        end
    end
    
points_sim(1:sim_bunch_counter,3)=1;
bunch_ids_sim_coord(1:sim_bunch_counter,3)=0;

end

    sim_bunch_counter_no = 1:sim_bunch_counter;
    bunch_coord_sim = points_sim;
    s = max(sim_grid(:,2))-min(sim_grid(:,2))+1;
    s1 = max(sim_grid(:,3))-min(sim_grid(:,3))+1; 

    % finds the ID of the bunch of nodes to be simulated
    bunch_ids_sim = (s).*(s1).*(points_sim(sim_bunch_counter_no,1)-1)+...
                               (points_sim(sim_bunch_counter_no,2)-1).*(s1)+...
                               (points_sim(sim_bunch_counter_no,3));

    % conditioing data from both primary cond data and previously simulated data are extracted
    [data_event_cond informed_nodes_id] = cond_dev_extract(window_type,[sim_x sim_y sim_z],...
                                          cond_data,search_radii,sim_dim(3));
    c = size(data_event_cond,1);
    
    % maximum number of nodes in search window. In circular case the n closest nodes are considered for further calculations
    if c > max_no_nodes_in_data_event
        max_no_nodes_in_data_event_no=1:max_no_nodes_in_data_event;
       dist_to_central_node = sqrt((data_event_cond(:,1)-sim_x).^2+(data_event_cond(:,2)-sim_y).^2+(data_event_cond(:,3)-sim_z).^2);
       [~,IX] = sort(dist_to_central_node);
       data_event_cond = cond_data(informed_nodes_id(IX(max_no_nodes_in_data_event_no)),1:5);
    end

    c = size(data_event_cond,1);
    c_no=1:c;
    
    % data-event is defined based on the lag vectors h
    h_x = data_event_cond(:,1)-sim_x;  
    h_y = data_event_cond(:,2)-sim_y;   
    h_z = data_event_cond(:,3)-sim_z;

%% Define a search window in the TI grid by using the dimensions of the data event.
    h_x_pos = h_x(h_x>0); h_x_neg = h_x(h_x<0);
    h_y_pos = h_y(h_y>0); h_y_neg = h_y(h_y<0);
    h_z_pos = h_z(h_z>0); h_z_neg = h_z(h_z<0);

    if isempty(h_x_pos)==1; h_x_pos=0;end
    if isempty(h_x_neg)==1; h_x_neg=0;end
    if isempty(h_y_pos)==1; h_y_pos=0;end
    if isempty(h_y_neg)==1; h_y_neg=0;end
    if isempty(h_z_pos)==1; h_z_pos=0;end
    if isempty(h_z_neg)==1; h_z_neg=0;end

    x_up = max(h_x_pos); x_down = abs(min(h_x_neg));
    y_up = max(h_y_pos); y_down = abs(min(h_y_neg));
    z_up = max(h_z_pos); z_down = abs(min(h_z_neg));

    if isempty(x_up)==1; x_up=0;end
    if isempty(x_down)==1; x_down=0;end
    if isempty(y_up)==1; y_up=0;end
    if isempty(y_down)==1; y_down=0;end
    if isempty(z_up)==1; z_up=0;end
    if isempty(z_down)==1; z_down=0;end
    
    ti_active_id_x = and(ti_grid(:,1) >= (min(ti_grid(:,1))+x_down) , ti_grid(:,1) <= (max(ti_grid(:,1))-x_up));
    ti_active_id_y = and(ti_grid(:,2) >= (min(ti_grid(:,2))+y_down) , ti_grid(:,2) <= (max(ti_grid(:,2))-y_up));
    ti_active_id_z = and(ti_grid(:,3) >= (min(ti_grid(:,3))+z_down) , ti_grid(:,3) <= (max(ti_grid(:,3))-z_up));
 
    if ti_size(3)==1; ti_active_id_z=ones(size(ti_active_id_z,1),1);end
    
    ti_active_id = and(and(ti_active_id_x==1 , ti_active_id_y==1) , ti_active_id_z==1);
    
    % The sacanable area of the TI
    ti_active = ti_grid(ti_active_id,:);

    ti_active_width_x = max(ti_active(:,1))-min(ti_active(:,1))+1; 
    ti_active_width_y = max(ti_active(:,2))-min(ti_active(:,2))+1;
    ti_active_width_z = max(ti_active(:,3))-min(ti_active(:,3))+1;
    
    min_x = min(ti_active(:,1));
    min_y = min(ti_active(:,2));
    min_z = min(ti_active(:,3));

%% 2. scan the TI for matching pattern
   path_size = ti_active_width_x*ti_active_width_y*ti_active_width_z;
    if ti_path_mode == 2
        ti_path=1:path_size;
        path=ti_path;
        seed=randperm(path_size,1);
        ti_path(1:(length(path)-seed(1)))=path(seed(1)+1:end);
        ti_path((length(path)-seed(1))+1:end)=path(1:seed(1));    
    elseif ti_path_mode == 1
        ti_path = randperm(path_size);
    end

    % matching criteria is reached when Distance<=distance threshold

    count_thre = floor(size(ti_path,2)*fract_of_ti_to_scan)-1;
    counter = 0;
    mindist = inf;
            
        % Sacnning the TI starts here
        while 1==1

              counter = counter+1;
              ti_node_id = ti_path(counter);
              point_ti = ti_active(ti_node_id,1:3);
              
              data_event_ti_coord_x = ti_active(ti_node_id,1)+h_x;
              data_event_ti_coord_y = ti_active(ti_node_id,2)+h_y;
              data_event_ti_coord_z = ti_active(ti_node_id,3)+h_z;
              
              if ti_size(3) == 1
                  data_event_ti_coord_z = ones(size(data_event_ti_coord_x,1),1);
              end
              
              % Extracting a bunch of nodes around the central node in the ti data event              
              Points_x = point_ti(1,1) + bunch_ids_sim_coord(sim_bunch_counter_no,1);
              Points_y = point_ti(1,2) + bunch_ids_sim_coord(sim_bunch_counter_no,2);
              Points_z = point_ti(1,3) + bunch_ids_sim_coord(sim_bunch_counter_no,3);
           
              % the ID of the nodes of the dev_ti which corresponds the data event from simulation grid
               
              nodes_id_ti = (ti_active_width_x).*(ti_active_width_y).*(Points_z-(min_z))+...
                                                 (Points_y-(min_y)).*(ti_active_width_x)+(Points_x)-(min_x)+1;     

              bunch_ids_ti = nodes_id_ti(and(nodes_id_ti <= size(ti_active,1) , nodes_id_ti>0),1);
              ti_bunch_counter=size(bunch_ids_ti,1);

              if ti_bunch_counter==0;break;end
    
              min_size = 1:min(ti_bunch_counter , sim_bunch_counter);
              ti_bunch = ti_active(bunch_ids_ti,4);

              % the case when data event from simulation grid is empty: a random pattern is pasted
              if c < 1
                 c_empty_dv = c_empty_dv+1;
                 sim_grid(bunch_ids_sim(min_size),4) = ti_bunch;
                 
                 nodes_id_sim = (ti_active_width_y).*(ti_active_width_z).*(sim_grid(bunch_ids_sim,1)-(floor(search_radius_x/2)+1))+...
                                                (sim_grid(bunch_ids_sim,2)-(floor(search_radius_y/2)+1)).*(ti_active_width_z)+...
                                                (sim_grid(bunch_ids_sim,3)-(floor(search_radius_z/2)+1))+1;
                 nodes_id_sim = nodes_id_sim(and((nodes_id_sim)>0 , nodes_id_sim<=size(active,1)),1);
                 active(nodes_id_sim,5) = -1;                
                 break
              end
              
              % based on the lag vectors defined above the data_event_ti is made then
              data_event_ti = zeros(c,1);
              for ii=1:c
                  data_event_ti(ii,1)=ti(data_event_ti_coord_x(ii,1),data_event_ti_coord_y(ii,1),data_event_ti_coord_z(ii,1));
              end
              
              % evaluating the distance between the data event in the simulation and that of ti              
              if Distance_type == 2
%                   d_max = max(abs(data_event_ti(c_no,1)))-min(abs(data_event_ti(c_no,1)));
                  d_max = max(abs(data_event_ti(c_no,1)-data_event_cond(c_no,4)));
                  if d_max ~= 0
                      Distance = (sum(abs((data_event_ti(c_no,1)-data_event_cond(c_no,4)).*data_event_cond(c_no,5))./d_max))/(sum(data_event_cond(c_no,5)));
                  elseif d_max == 0
                     Distance = 0;
                  end
              elseif Distance_type == 1
                      Distance = mean(data_event_cond(c_no,4)~=data_event_ti(c_no,1));
              end
              
              % The minimun distance calculated so fat is updated              
              if Distance < mindist
                 mindist = Distance;
                 best_bunch = ti_bunch;
                 best_bunch_no = 1:size(best_bunch,1);
              end

              % the case when the pre-specified fraction of TI is scanned and the matching pattern is not found yet
              if counter > count_thre
                 c_count_thre = c_count_thre+1;
                 id = (ti_active_width_y).*(bunch_coord_sim(best_bunch_no,1)-(floor(search_radius_x/2)+1))+...
                           (bunch_coord_sim(best_bunch_no,2)-(floor(search_radius_y/2)));                
                 id = id(and(id <= size(active,1),id>0));
                 active(id,5) = -1;
                 sim_grid(bunch_ids_sim(best_bunch_no),4) = best_bunch;
              break
    
              % the case when the matching pattern is found that satisfies the distance threshold
              elseif Distance <= dist_threshold
                 c_distance = c_distance+1;               
                 id2 = (ti_active_width_y).*(ti_active_width_z).*(bunch_coord_sim(best_bunch_no,1)-(floor(search_radius_x/2)+1))+...
                                       (bunch_coord_sim(best_bunch_no,2)-(floor(search_radius_y/2)+1)).*(ti_active_width_z)+...
                                       (bunch_coord_sim(best_bunch_no,3)-(floor(search_radius_z/2)+1))+1;                
                 id2 = id2(and(id2 <= size(active,1),id2>0));
                 active(id2,5) = -1;
                 sim_grid(bunch_ids_sim(best_bunch_no),4) = ti_bunch;
              break
              end

        end

    % with the addittion of simulated nodes the conditioning data is inflated each time
    cond_data = sim_grid(and(isnan(sim_grid(:,4))==0,sim_grid(:,4)~=-2),:);

    % The illustration of realization while doing the simulation
    if simulation_show_mode == 1
        if length(find((sim_dim)==1))==1
        if sim_dim(3) == 1
            id_level=sim_grid(:,3) == sim_z;       
            sim = reshape(sim_grid(id_level==1,4),sim_dim(2)+2*floor(search_radius_y/2),...
                                                  sim_dim(1)+2*floor(search_radius_x/2));
            shg;imagesc(sim);axis xy;axis equal;xlim([1 sim_dim(1)+2*floor(search_radius_x/2)]);
            ylim([1 sim_dim(2)+2*floor(search_radius_y/2)]);
            title('Bunch DS simulation-XY Plane')
        end
        if sim_dim(2) == 1
            id_level=sim_grid(:,2) == sim_y;       
            sim = reshape(sim_grid(id_level==1,4),sim_dim(3)+2*floor(search_radius_z/2),...
                                                  sim_dim(1)+2*floor(search_radius_x/2));
            shg;imagesc(sim);axis xy;axis equal;xlim([1 sim_dim(1)+2*floor(search_radius_x/2)]);
            ylim([1 sim_dim(3)+2*floor(search_radius_z/2)]);
            title('Bunch DS simulation-XZ Plane')
        end
        if sim_dim(1) == 1
            id_level=sim_grid(:,1) == sim_x;       
            sim = reshape(sim_grid(id_level==1,4),sim_dim(3)+2*floor(search_radius_z/2),...
                                                  sim_dim(2)+2*floor(search_radius_y/2));
            shg;imagesc(sim);axis xy;axis equal;xlim([1 sim_dim(2)+2*floor(search_radius_y/2)]);
            ylim([1 sim_dim(3)+2*floor(search_radius_z/2)]);
            title('Bunch DS simulation-YZ Plane')
        end
        end
    end
end
end
active = active(and(active(:,5)~=-1,active(:,1)~=0),:);
    
    if size(active,1)>=1
       sim_path = randperm(size(active,1));
    end

% checks for the criteria that determines the end of the simulation 
if or( sum(active(:,5)~=-1)<=diff , sum((isnan(sim_grid(:,4)))==0) )
   T_one(real,1)=toc;
   toc
   disp(['--------------------------','END OF SIMULATION USING BUNCH-DS-Rezaee et al., 2013'])
   break
end

end

Realization(:,real)=sim_grid(:,4);
sim_grid = SIM_GRID;
end
T_total=sum(T_one);

%% Illustration of outputs

nan_id = (Realization(:,1))~=-2;
Realizations(:,1:3)=sim_grid(nan_id==1,[1 2 3]);

for r=1:No_Realizations
    Realizations(:,r+3) = Realization(nan_id==1,r);
end

Realizations(:,1)=Realizations(:,1)-(min(Realizations(:,1))-1);
Realizations(:,2)=Realizations(:,2)-(min(Realizations(:,2))-1);
Realizations(:,3)=Realizations(:,3)-(min(Realizations(:,3))-1);

% Re-format the output realization to the GSlib format
sim1=zeros(sim_dim(1)*sim_dim(2)*sim_dim(3),No_Realizations+3);
for r=1:No_Realizations
counter=0;
for k=1:sim_dim(3)
    for j=1:sim_dim(2)
        for i=1:sim_dim(1)
            counter=counter+1;
            id=and(and(Realizations(:,1)==i,Realizations(:,2)==j),Realizations(:,3)==k);
            sim1(counter,1:3) = [i,j,k];
            sim1(counter,r+3) = Realizations(id==1,r+3);
        end
    end
end
end
Final_Realizations = sim1;
% the illustration of TI
if min(sim_dim)==1
    
       subplot(1,2,1)
       
       if ti_size(3)==1
           imagesc(ti);axis xy;axis equal;xlim([1 ti_size(2)]);ylim([1 ti_size(1)]);colorbar;title('Training Image')
       elseif ti_size(3)>1
           imagesc(ti(:,:,fix(mean(ti_grid(:,3)))));axis xy;axis equal;xlim([1 ti_size(2)]);ylim([1 ti_size(1)]);
           colorbar;title('Training Image')
       end
       if No_facies<=2;colormap('gray');end

       subplot(1,2,2)
       if sim_dim(3) == 1
           imagesc(reshape(sim1(:,4),sim_dim(1),sim_dim(2))');axis xy;axis equal;xlim([1 sim_dim(1)+0.1]);ylim([1 sim_dim(2)+0.1]);colorbar
       end
       if sim_dim(2) == 1
           imagesc(reshape(sim1(:,4),sim_dim(1),sim_dim(3))');axis xy;axis equal;xlim([1 sim_dim(1)+0.1]);ylim([1 sim_dim(3)+0.1]);colorbar
       end
       if sim_dim(1) == 1
           imagesc(reshape(sim1(:,4),sim_dim(2),sim_dim(3))');axis xy;axis equal;xlim([1 sim_dim(2)+0.1]);ylim([1 sim_dim(3)+0.1]);colorbar
       end
       title('Simulation')
       
       if Condition_mode == 1 
           if size(Hard_within,1)>1
           hold on 
           hard_data = Hard;
           hard_data(:,1:3)=hard_data(:,1:3)+1;
           scatter(hard_data(:,1),hard_data(:,2) , 'o' , 'MarkerEdgeColor','r');
           legend('hard data','Location','NorthWest')
           end
       end

       if No_facies<=2;colormap('gray');end

end

if min(sim_dim)>1
    
       subplot(1,3,1)
       id_level=sim1(:,1) == fix(mean(sim1(:,1)));
       sim_YZ_section = reshape(sim1(id_level==1,4),sim_dim(3),sim_dim(2));
       imagesc(sim_YZ_section);axis xy;axis equal;xlim([1 sim_dim(2)]);
       ylim([1 sim_dim(3)]);
       title(['Realz#1 at YZ Plane,  X = ' num2str(fix(mean(sim1(:,1))))])
       if No_facies<=2;colormap('gray');end
       
       subplot(1,3,2)
       id_level=sim1(:,2) == fix(mean(sim1(:,2)));
       sim_XZ_section = reshape(sim1(id_level==1,4),sim_dim(3),sim_dim(1));
       imagesc(sim_XZ_section);axis xy;axis equal;xlim([1 sim_dim(1)]);
       ylim([1 sim_dim(3)]);
       title(['Realz#1 at XZ Plane,  Y = ' num2str(fix(mean(sim1(:,1))))])
       if No_facies<=2;colormap('gray');end
       
       subplot(1,3,3)
       id_level=sim1(:,3) == fix(mean(sim1(:,3)));
       sim_XY_section = reshape(sim1(id_level==1,4),sim_dim(1),sim_dim(2));
       imagesc(sim_XY_section);axis xy;axis equal;xlim([1 sim_dim(2)]);
       ylim([1 sim_dim(1)]);
       title(['Realz#1 at XY Plane,  Z = ' num2str(fix(mean(sim1(:,1))))])
       if No_facies<=2;colormap('gray');end      
end

disp(' ')
disp([num2str((c_empty_dv/(c_empty_dv+c_count_thre+c_distance)*100)) ' % of nodes are simulated by pasting a random pattern from TI'])
disp([num2str((c_count_thre/(c_empty_dv+c_count_thre+c_distance)*100)) ' % of nodes are simulated by drawing the pattern of minimum distance found so far'])
disp([num2str((c_distance/(c_empty_dv+c_count_thre+c_distance)*100)) ' % of nodes are simulated by pasting the matching pattern satisfying the distance threshold'])

stats(1:3,1) = [(c_empty_dv/(c_empty_dv+c_count_thre+c_distance)*100)...
              (c_count_thre/(c_empty_dv+c_count_thre+c_distance)*100)...
              (c_distance/(c_empty_dv+c_count_thre+c_distance)*100)];

if and(stats(1)>10 , stats(3)>10)
    explode = [1 1 1];
    figure
    pie(stats,explode,{'random pattern',...
        'minimum distance',...
        'matching pattern'})
end
%% Carrying out the post processing steps on the realizations produced

if and(post_proc_mode == 1 , No_Realizations>1)

   post = sim_post_proc(sim1(:,4:end));
   post_maps(:,1:3) = sim1(:,1:3);
   post_maps(:,4:6) = post;
   
if min(sim_dim)==1
    if sim_dim(3) ==1
        etype = reshape(post(:,1), sim_dim(1) , sim_dim(2));
        figure
        subplot(1,3,1)
        imagesc(etype);axis xy;axis equal;xlim([1 sim_dim(2)]);ylim([1 sim_dim(1)]);title('E-type')
           if Condition_mode == 1
            if size(Hard_within,1)>1
                hold on
                scatter(hard_data(:,1),hard_data(:,2),'o')
                legend('hard data','Location','NorthWest')
            end
           end 
        
        subplot(1,3,2)
        IQR = reshape(post (:,2), sim_dim(1) , sim_dim(2));
        imagesc(IQR);axis xy;axis equal;xlim([1 sim_dim(2)]);ylim([1 sim_dim(1)]);title('IQR')
        if Condition_mode == 1
            if size(Hard_within,1)>1
                hold on
                scatter(hard_data(:,1),hard_data(:,2),'o')
                legend('hard data','Location','NorthWest')
            end
        end       
        subplot(1,3,3)
        Quantile = reshape(post (:,3), sim_dim(1) , sim_dim(2));
        imagesc(Quantile);axis xy;axis equal;xlim([1 sim_dim(2)]);ylim([1 sim_dim(1)]);title('Quantile (probability 0.5)')
        if Condition_mode == 1
            if size(Hard_within,1)>1
                hold on
                scatter(hard_data(:,1),hard_data(:,2),'o')
                legend('hard data','Location','NorthWest')
            end
        end
     end
   
if sim_dim(2) ==1
    etype = reshape(post(:,1), sim_dim(1) , sim_dim(3));
    figure
    subplot(1,3,1)
    imagesc(etype);axis xy;axis equal;xlim([1 sim_dim(3)]);ylim([1 sim_dim(1)]);title('E-type')
    subplot(1,3,2)
    IQR = reshape(post (:,2), sim_dim(1) , sim_dim(3));
    imagesc(IQR);axis xy;axis equal;xlim([1 sim_dim(3)]);ylim([1 sim_dim(1)]);title('IQR')
    if Condition_mode == 1
        if size(Hard_within,1)>1
            hold on
            scatter(hard_data(:,1),hard_data(:,3),'o')
            legend('hard data','Location','NorthWest')
        end
    end       
    subplot(1,3,3)
    Quantile = reshape(post (:,3), sim_dim(1) , sim_dim(3));
    imagesc(Quantile);axis xy;axis equal;xlim([1 sim_dim(3)]);ylim([1 sim_dim(1)]);title('Quantile (probability 0.5)')
    if Condition_mode == 1
        if size(Hard_within,1)>1
            hold on
            scatter(hard_data(:,1),hard_data(:,3),'o')
            legend('hard data','Location','NorthWest')
        end
    end
end
    
if sim_dim(1) ==1
    etype = reshape(post(:,1), sim_dim(2) , sim_dim(3));
    figure
    subplot(1,3,1)
    imagesc(etype);axis xy;axis equal;xlim([1 sim_dim(3)]);ylim([1 sim_dim(2)]);title('E-type')
    subplot(1,3,2)
    IQR = reshape(post (:,2), sim_dim(2) , sim_dim(3));
    imagesc(IQR);axis xy;axis equal;xlim([1 sim_dim(3)]);ylim([1 sim_dim(2)]);title('IQR')
    if Condition_mode == 1
        if size(Hard_within,1)>1
            hold on
            scatter(hard_data(:,2),hard_data(:,3),'o')
            legend('hard data','Location','NorthWest')
        end
    end       
    subplot(1,3,3)
    Quantile = reshape(post (:,3), sim_dim(2) , sim_dim(3));
    imagesc(Quantile);axis xy;axis equal;xlim([1 sim_dim(3)]);ylim([1 sim_dim(2)]);title('Quantile (probability 0.5)')
    if Condition_mode == 1
        if size(Hard_within,1)>1
            hold on
            scatter(hard_data(:,2),hard_data(:,3),'o')
            legend('hard data','Location','NorthWest')
        end
    end
end    
    
elseif min(sim_dim)>1
       subplot(1,3,1)
       id_level=sim1(:,1) == fix(mean(sim1(:,1)));
       etype_YZ_section = reshape(post(id_level==1,1),sim_dim(3),sim_dim(2));
       imagesc(etype_YZ_section);axis xy;axis equal;xlim([1 sim_dim(2)]);
       ylim([1 sim_dim(3)]);colorbar
       title(['E-type at YZ Plane,  X = ' num2str(fix(mean(sim1(:,1))))])
       
       subplot(1,3,2)
       id_level=sim1(:,2) == fix(mean(sim1(:,2)));
       iqr_XZ_section = reshape(post(id_level==1,2),sim_dim(3),sim_dim(1));
       imagesc(iqr_XZ_section);axis xy;axis equal;xlim([1 sim_dim(1)]);
       ylim([1 sim_dim(3)]);colorbar
       title(['IQR at XZ Plane,  Y = ' num2str(fix(mean(sim1(:,1))))])
       
       subplot(1,3,3)
       id_level=sim1(:,3) == fix(mean(sim1(:,3)));
       quantile_XY_section = reshape(post(id_level==1,3),sim_dim(1),sim_dim(2));
       imagesc(quantile_XY_section);axis xy;axis equal;xlim([1 sim_dim(2)]);
       ylim([1 sim_dim(1)]);colorbar
       title(['Quantile at XY Plane,  Z = ' num2str(fix(mean(sim1(:,1))))])
end
elseif post_proc_mode ~= 1
    post_maps=nan;
end

%% To save the input parameter file
parfile=cell(18,4);
    parfile(1,1) = {'TI name = '};
    parfile(1,2) = {ti_name};
    parfile(2,1) = {'Conditioning mode = '};
    if Condition_mode == 1
        parfile(2,2) = {'Conditional Simulation'};
    elseif Condition_mode == 0
        parfile(2,2) = {'Non-Conditional Simulation'};
    end    
    parfile(3,1) = {'conditioning weight  = '};
    if Condition_mode == 1
        parfile(3,2) = num2cell(conditioning_weight);
    elseif Condition_mode == 0
        parfile(3,2) = {'Non-Conditional Simulation'};
    end   
    parfile(4,1) = {'Conditioning data name = '};
    parfile(4,2) = {cond_data_name};    
    parfile(5,1) = {'Simulation Grid Dimension (x,y) = '};
    parfile(5,2:4) = num2cell(sim_dim);    
    parfile(6,1) = {'Search Radii (X , Y , Z) = '};
    parfile(6,2) = num2cell(search_radius_x);
    parfile(6,3) = num2cell(search_radius_y);
    parfile(6,4) = num2cell(search_radius_z);    
    parfile(7,1) = {'Fraction of ti which is scanned = '};
    parfile(7,2) = num2cell(fract_of_ti_to_scan);    
    parfile(8,1) = {'Distance type = '};
    if Distance_type==1
        parfile(8,2) = {'Categorical variable'};
    elseif Distance_type==2
        parfile(8,2) = {'Continuous variable'};
    end    
    parfile(9,1) = {'Distance threshold = '};
    parfile(9,2) = num2cell(dist_threshold);    
    parfile(10,1) = {'Maximum number of nodes in dataevent = '};
    parfile(10,2) = num2cell(max_no_nodes_in_data_event);    
    parfile(11,1) = {'Window type = '};
    if window_type==1
        parfile(11,2) = {'Circular'};
    elseif window_type==2 
        parfile(11,2) = {'Squared'};
    elseif window_type == 3
        parfile(11,2) = {'Disc Shape'};
    end    
    parfile(12,1) = {'Bunch size = '};
    if min(sim_dim)>1
        parfile(12,2) = num2cell((2*bunch_size+1)^3);
    elseif min(sim_dim)==1
        parfile(12,2) = num2cell((2*bunch_size+1)^2);
    end  
    parfile(13,1) = {'Simulation path mode = '};
    if sim_path_mode == 1
        parfile(13,2) = {'Nodes are visited "Randomly"'};
    elseif sim_path_mode == 2
        parfile(13,2) = {'Unilateral or Raster path is adopted'};
    end
    
    parfile(14,1) = {'The TI scanned through = '};
    if ti_path_mode == 2
        parfile(14,2) = {'Unilateral path with starting point'};
    elseif ti_path_mode == 1
        parfile(14,2) = {'Random path'};
    end    
    
    parfile(15,1) = {'Number of realizations that are produced = '};
    parfile(15,2) = num2cell(No_Realizations);     
    parfile(16,1) = {'Post-processing maps = '};
    if post_proc_mode == 1
        parfile(16,2) = {'E-type, IQR and Quantile maps are produced"'};
    elseif post_proc_mode ~= 1
        parfile(16,2) = {'No post-processing step is taken'};
    end    
    parfile(17,1) = {'---------------------------------------------------'};    
    parfile(18,1) = {'Simulation time for one realization is (sec): ' };
    parfile(18,2) = num2cell(T_one(1));    
    if No_Realizations>1
        parfile(19,1) = {'Simulation time for all realizations is (sec): ' };
        parfile(19,2) = num2cell(T_total);
    end
dlmcell('Parfile.txt',parfile) 
