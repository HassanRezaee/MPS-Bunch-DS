%------------------------------------------------------------------------------------------------
% Author: Hassan Rezaee, Email: h.rezaee@ut.ac.ir
% Paper: "Multiple-Point Geostatistical Simulation Using the Bunch-Pasting Direct Sampling Method"
%------------------------------------------------------------------------------------------------

% Thank you for your interest in Bunch DS!

% To reach a full understanding of the input parameters please read User Guide where
% the main function Bunch_DS.m is introduced and its inputs and outputs are explained further. 
% The present code is written in 7.14.0.739 (R2012a). Using former versions
% may result in errors while running the simulation. Please beware of the
% minor differences between Matlab versions scripting points.

clc
clear all
close all

%% input parameters

ti_name = 'TI_A';   %  Three examples are provided: "TI_A", "TI_D", "TI_E" (the last one is 3D)
Condition_mode = 1;   % 1. conditional simulation 0. non-conditional simulation
conditioning_weight = 10;   % the conditioning data wieght
cond_data_name = 'TI_A_Sparse_Cond';   % The samples provided are "TI_A_Cond", "TI_A_Sparse_Cond", "TI_D_Cond", "TI_E_Cond"
sim_dim = [64 64 1];   % Along X, Y and Z respectively
Search_Radii = [12 12 1];   % Along X, Y and Z respectively
fract_of_ti_to_scan = 0.7;
dist_threshold = 0.01;
max_no_nodes_in_data_event = 50;
window_type = 1;   % 1. Circular 2. Squared 3. Disc-shape
bunch_size = 2;   % In 2D the number of nodes in bunch = (2*bunch_size + 1)^2 and in 3D = (2*bunch_size + 1)^3
ti_path_mode = 2;   % 1. random path 2. unilateral with starting point
sim_path_mode = 2;   % 1. random path 2. unilateral
No_Realizations = 1;
post_proc_mode = 0;   % 0. No post-processing 1. post processing steps to calculate Etype-IQR and Quantile maps
simulation_show_mode = 1;   % 1. Show the realization while doing the simulation, any other value. does not show

%% the main function
[Final_Realizations post_maps] = Bunch_DS(ti_name,Condition_mode,conditioning_weight,cond_data_name,sim_dim,Search_Radii,fract_of_ti_to_scan,...
                           dist_threshold,max_no_nodes_in_data_event,window_type,bunch_size,sim_path_mode,ti_path_mode,...
                           No_Realizations,post_proc_mode,simulation_show_mode);
                       