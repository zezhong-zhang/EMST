% output_multislice = il_multem(system_conf, input_multem) perform TEM simulation
% Convergent beam electron diffraction (CBED) simulation
% All parameters of the input_multem structure are explained in ilm_dflt_input_multem()
% Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>

clear; clc;
addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])

%%%%%%%%%%%%%%%%%% Load multem default parameter %%%%%%%%$$%%%%%%%%%
input_multem = ilm_dflt_input_multem();          % Load default values;

%%%%%%%%%%%%%%%%%%%%% Set system configuration %%%%%%%%%%%%%%%%%%%%%
system_conf.precision = 1;                           % eP_Float = 1, eP_double = 2
system_conf.device = 2;                              % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 1;
system_conf.gpu_device = 0;

%%%%%%%%%%%%%%%%%%%% Set simulation experiment %%%%%%%%%%%%%%%%%%%%%
% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52,
% eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
input_multem.simulation_type = 21;

%%%%%%%%%%%%%% Electron-Specimen interaction model %%%%%%%%%%%%%%%%%
input_multem.interaction_model = 1;              % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multem.potential_type = 6;                 % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

%%%%%%%%%%%%%%%%%%%%%%% Potential slicing %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.potential_slicing = 1;              % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4

%%%%%%%%%%%%%%% Electron-Phonon interaction model %%%%%%%%%%%%%%%%%%
input_multem.pn_model = 3;                       % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multem.pn_coh_contrib = 0;
input_multem.pn_single_conf = 0;                 % 1: true, 0:false (extract single configuration)
input_multem.pn_nconf = 10;                      % true: specific phonon configuration, false: number of frozen phonon configurations
input_multem.pn_dim = 110;                       % phonon dimensions (xyz)
input_multem.pn_seed = 300183;                   % Random seed(frozen phonon)

%%%%%%%%%%%%%%%%%%%% Microscope parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.E_0 = 200;                          % Acceleration Voltage (keV)
input_multem.theta = 0.0;                        % Till ilumination (�)
input_multem.phi = 0.0;                          % Till ilumination (�)

%%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.illumination_model = 4;             % 1: coherente mode, 4: Numerical integration
input_multem.temporal_spatial_incoh = 1;         % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_m = 0;                   % Vortex momentum
input_multem.cond_lens_c_10 = 1110;             % Defocus (�)
input_multem.cond_lens_c_30 = 3.3;              % Third order spherical aberration (mm)
input_multem.cond_lens_c_50 = 0.00;             % Fifth order spherical aberration (mm)
input_multem.cond_lens_c_12 = 0.0;              % Twofold astigmatism (�)
input_multem.cond_lens_phi_12 = 0.0;            % Azimuthal angle of the twofold astigmatism (�)
input_multem.cond_lens_c_23 = 0.0;              % Threefold astigmatism (�)
input_multem.cond_lens_phi_23 = 0.0;            % Azimuthal angle of the threefold astigmatism (�)
input_multem.cond_lens_inner_aper_ang = 0.0;    % Inner aperture (mrad)
input_multem.cond_lens_outer_aper_ang = 21.0;   % Outer aperture (mrad)

%%%%%%%%% defocus spread function %%%%%%%%%%%%
dsf_sigma = il_iehwgd_2_sigma(32); % from defocus spread to standard deviation
input_multem.cond_lens_ti_a = 1.0;                          % Height proportion of a normalized Gaussian [0, 1]
input_multem.cond_lens_ti_sigma = dsf_sigma;                % Standard deviation of the defocus spread for the Gaussian component (�)
input_multem.cond_lens_ti_beta = 0.0;                 		% Standard deviation of the defocus spread for the Exponential component (�)
input_multem.cond_lens_ti_npts = 4;                         % Number of integration points. It will be only used if illumination_model=4

%%%%%%%%%% source spread function %%%%%%%%%%%%
ssf_sigma = il_hwhm_2_sigma(0.45); % half width at half maximum to standard deviation
input_multem.cond_lens_si_a = 1.0;                          % Height proportion of a normalized Gaussian [0, 1]
input_multem.cond_lens_si_sigma = ssf_sigma;                % Standard deviation of the source spread function for the Gaussian component: For parallel ilumination(�^-1); otherwise (�)
input_multem.cond_lens_si_beta = 0.0;                 		% Standard deviation of the source spread function for the Exponential component: For parallel ilumination(�^-1); otherwise (�)
input_multem.cond_lens_si_rad_npts = 4;                     % Number of radial integration points. It will be only used if illumination_model=4
input_multem.cond_lens_si_azm_npts = 4;                     % Number of radial integration points. It will be only used if illumination_model=4

%%%%%%%%% zero defocus reference %%%%%%%%%%%%
input_multem.cond_lens_zero_defocus_type = 1;   % eZDT_First = 1, eZDT_User_Define = 4
input_multem.cond_lens_zero_defocus_plane = 0;

%%%%%%%%%%%%%%%%%%%%%%%% Objective lens %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% zero defocus reference %%%%%%%%%%%%
input_multem.obj_lens_zero_defocus_type = 3;    % eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
input_multem.obj_lens_zero_defocus_plane = 0;   % It will be only used if obj_lens_zero_defocus_type = eZDT_User_Define

%% Task start heare
%%%%%%%%%%%%%%%%%%%%%% Specimen information %%%%%%%%%%%%%%%%%%%%%%%
n = 8; a=1.42569; rms3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_lx...
, input_multem.spec_ly, input_multem.spec_lz] = graphene(n, a, rms3d);

ilm_show_crystal(1, input_multem.spec_atoms)

%%%%%%%%%%%%%%%%%%%%%% Specimen thickness %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.thick_type = 1;                     % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
input_multem.thick = 0;   % Array of thickes (�)

%%%%%%%%%%%%%%%%%%%%%% Potential sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_max = 21/20*32/1000; % 21 mrad, 20 pixels, 32 max pixels, 1000 from mrad to rad
lambda = wave_length(200);
pixel_size = 1/3*lambda/theta_max;

input_multem.nx = ceil(input_multem.spec_lx/pixel_size);
input_multem.ny = ceil(input_multem.spec_lx/pixel_size);
input_multem.bwl = 1;                            % Band-width limit, 1: true, 0:false

%%%%%%%%%%%%%%%%%%%%%% Probe sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sannnig area
nx_start=3*3*a;
nx_end=4*3*a;
ny_start=3*sqrt(3)*a;
ny_end=4*sqrt(3)*a;
% probe position
step = 0.2;
nx=nx_start:step:nx_end;
ny=ny_start:step:ny_end;

%%%%%%%%%%%%%%%%%%%%%% Start loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count=0; % set up the counting
cbed_stack={}; % cell arry to hold the data
dose = 1e5; % dose
mkdir('output/graphene') % folder to store the data
figure()
for j = 1: length(ny) % because the scanning direction and array x-y are swtiched
    for i= 1:length(nx)
        count=count+1;
        fprintf('Calculating...%d,%d...%.2f%% \n',j,i,count/length(nx)/length(ny)*100)
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
        input_multem.iw_type = 2;                        % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
        input_multem.iw_psi = 0;    % user define incident wave
        input_multem.iw_x = nx(i);     % x position 
        input_multem.iw_y = ny(j);     % y position
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate MULTEM %%%%%%%%%%%%%%%%%%%%%%%%%%
        clear il_MULTEM;
        output_multislice = il_multem(system_conf, input_multem);
        cbed=(output_multislice.data.m2psi_tot).*dose;
        cbed= fourier_interpolate(cbed, [64 64]);
        cbed_stack{j,i}=cbed;
        imagesc(cbed)
        pause(0.5)
        filename=strcat('output/','graphene/',num2str(i),'_',num2str(j),'.tif');
        imwrite(cbed,filename)
    end
end

