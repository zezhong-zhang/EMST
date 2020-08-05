% output_multislice = il_multem(system_conf, input_multem) perform TEM simulation
% Scanning transmission electron microscopy (STEM) simulation
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
system_conf.cpu_nthread = 6;
system_conf.gpu_device = 0;

%%%%%%%%%%%%%%%%%%%% Set simulation experiment %%%%%%%%%%%%%%%%%%%%%
% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52,
% eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
input_multem.simulation_type = 11;

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

%%%%%%%%%%%%%%%%%%%%%%% Specimen information %%%%%%%%%%%%%%%%%%%%%%%
na = 8; nb = 8; nc = 5; ncu = 2; rms3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_lx...
, input_multem.spec_ly, input_multem.spec_lz...
, a, b, c, input_multem.spec_dz] = Au001_xtl(na, nb, nc, ncu, rms3d);

%%%%%%%%%%%%%%%%%%%%%% Specimen thickness %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.thick_type = 2;                     % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
input_multem.thick = c/2:c:1000;                   % Array of thickes (�)

%%%%%%%%%%%%%%%%%%%% Microscope parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.E_0 = 300;                          % Acceleration Voltage (keV)
input_multem.theta = 0.0;                        % Till ilumination (�)
input_multem.phi = 0.0;                          % Till ilumination (�)

%%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.illumination_model = 1;             % 1: coherente mode, 4: Numerical integration
input_multem.temporal_spatial_incoh = 1;         % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_m = 0;                  % Vortex momentum
input_multem.cond_lens_c_10 = -5.32;            % Defocus (�)
input_multem.cond_lens_c_30 = 1e-03;            % Third order spherical aberration (mm)
input_multem.cond_lens_c_50 = 0.00;             % Fifth order spherical aberration (mm)
input_multem.cond_lens_c_12 = 0.0;             % Twofold astigmatism (�)
input_multem.cond_lens_phi_12 = 0.0;             % Azimuthal angle of the twofold astigmatism (�)
input_multem.cond_lens_c_23 = 0.0;             % Threefold astigmatism (�)
input_multem.cond_lens_phi_23 = 0.0;             % Azimuthal angle of the threefold astigmatism (�)
input_multem.cond_lens_inner_aper_ang = 0.0;   % Inner aperture (mrad)
input_multem.cond_lens_outer_aper_ang = 21.0;  % Outer aperture (mrad)

%%%%%%%%% defocus spread function %%%%%%%%%%%%
dsf_sigma = il_iehwgd_2_sigma(32); % from defocus spread to standard deviation
input_multem.cond_lens_ti_sigma = dsf_sigma;   % standard deviation (�)
input_multem.cond_lens_ti_npts = 5;         % # of integration points. It will be only used if illumination_model=4

%%%%%%%%%% source spread function %%%%%%%%%%%%
ssf_sigma = il_hwhm_2_sigma(0.45); % half width at half maximum to standard deviation
input_multem.cond_lens_si_sigma = ssf_sigma;  	% standard deviation: For parallel ilumination(�^-1); otherwise (�)
input_multem.cond_lens_si_rad_npts = 4;         % # of integration points. It will be only used if illumination_model=4

%%%%%%%%% zero defocus reference %%%%%%%%%%%%
input_multem.cond_lens_zero_defocus_type = 1;   % eZDT_First = 1, eZDT_User_Define = 4
input_multem.cond_lens_zero_defocus_plane = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%STEM probe sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = wave_length(300);
resolution = lambda/4/input_multem.cond_lens_outer_aper_ang*1000;
step = ceil(a/resolution);

input_multem.scanning_type = 2; % eST_Line = 1, eST_Area = 2
input_multem.scanning_periodic = 1;     % 1: true, 0:false (periodic boundary conditions)
input_multem.scanning_ns = step;
input_multem.scanning_x0 = 3*a;
input_multem.scanning_y0 = 3*b;
input_multem.scanning_xe = 4*a;
input_multem.scanning_ye = 4*b;

%%%%%%%%%%%%%%%%%%%%%% x-y sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_max = 0.2;
dx = 1/3/theta_max*lambda;

input_multem.nx = pow2(ceil(log2(input_multem.spec_lx/dx)));
input_multem.ny = pow2(ceil(log2(input_multem.spec_ly/dx)));
input_multem.bwl = 0;                            % Band-width limit, 1: true, 0:false


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Detector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.detector.type = 1;  % eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3
input_multem.detector.cir(1).inner_ang = 0;  % Inner angle(mrad)
input_multem.detector.cir(1).outer_ang = 15; % Outer angle(mrad)
input_multem.detector.cir(2).inner_ang = 50;  % Inner angle(mrad)
input_multem.detector.cir(2).outer_ang = 200; % Outer angle(mrad)

%%%%%%%%%%%%%%%%%%%%%%%%%% run multem %%%%%%%%%%%%%%%%%%%%%%%%%%
clear il_multem;
tic;
output_multislice = il_multem(system_conf, input_multem);
toc;

figure(1);
for i=1:length(output_multislice.data)
	ndet = length(input_multem.detector.cir);
    for j=1:ndet
        subplot(1, ndet, j);
        image = output_multislice.data(i).image_tot(j).image;
        image = kron(ones([3 3]),image);
        image = fourier_interpolate(image, [512 512]);
        imagesc(image);
        title(strcat('Thk = ', num2str(i), ', det = ', num2str(output_multislice.thick(j))));
        axis image;
        colormap gray;
    end
    pause(0.5);
end
