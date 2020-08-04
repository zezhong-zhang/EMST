% output_multislice = il_multem(system_conf, input_multem) perform TEM simulation
% Incident wave simulation
% All parameters of the input_multem structure are explained in ilm_dflt_input_multem()
% Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>

clear; clc;
addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])

input_multem = ilm_dflt_input_multem();         % Load default values;

system_conf.precision = 1;                     % eP_Float = 1, eP_double = 2
system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 4; 
system_conf.gpu_device = 0;

input_multem.E_0 = 300;                          % Acceleration Voltage (keV)
input_multem.theta = 0.0;
input_multem.phi = 0.0;

input_multem.spec_lx = 20;
input_multem.spec_ly = 20;

input_multem.nx = 1024; 
input_multem.ny = 1024;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.iw_type = 2;   % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.iw_psi = read_psi_0_multem(input_multem.nx, input_multem.ny);    % user define incident wave
input_multem.iw_x = 0.5*input_multem.spec_lx;    % x position 
input_multem.iw_y = 0.5*input_multem.spec_ly;    % y position

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_m = 0;                  % Vortex momentum
input_multem.cond_lens_c_10 = 0;          % Defocus (�)
input_multem.cond_lens_c_30 = 0;            % Third order spherical aberration (mm)
input_multem.cond_lens_c_50 = 0.00;             % Fifth order spherical aberration (mm)
input_multem.cond_lens_c_12 = 0.0;              % Twofold astigmatism (�)
input_multem.cond_lens_phi_12 = 0.0;          % Azimuthal angle of the twofold astigmatism (�)
input_multem.cond_lens_c_23 = 0.0;              % Threefold astigmatism (�)
input_multem.cond_lens_phi_23 = 0.0;          % Azimuthal angle of the threefold astigmatism (�)
input_multem.cond_lens_inner_aper_ang = 18;   % Inner aperture (mrad) 
input_multem.cond_lens_outer_aper_ang = 21.0;  % Outer aperture (mrad)
input_multem.cond_lens_ti_sigma = 0;                % standard deviation (�)
input_multem.cond_lens_ti_npts = 0;               % # of integration points. It will be only used if illumination_model=4
input_multem.cond_lens_si_sigma = 0.2;             % standard deviation: For parallel ilumination(�^-1); otherwise (�)
input_multem.cond_lens_si_rad_npts = 8;             % # of integration points. It will be only used if illumination_model=4
input_multem.cond_lens_zero_defocus_type = 1;  % eZDT_First = 1, eZDT_User_Define = 4
input_multem.cond_lens_zero_defocus_plane = 0;

defocus = linspace(-250,250,100);
probe=zeros(1024,1024,100);
input_multem.cond_lens_inner_aper_ang = 0;
for n = 1:100
    input_multem.cond_lens_c_10=defocus(n);
    output_incident_wave = il_incident_wave(system_conf, input_multem);
    psi_0 = output_incident_wave.psi_0;
    probe(:,:,n) = abs(psi_0).^2;
end
normal_slice=squeeze(probe(512,:,:));
subplot(1,2,1)
imagesc(linspace(-10,10,1024),defocus,normal_slice');
xlabel('distance (A)')
ylabel('Defocus (A)')

probe=zeros(1024,1024,100);
input_multem.cond_lens_inner_aper_ang = 18;
for n = 1:100
    input_multem.cond_lens_c_10=defocus(n);
    output_incident_wave = il_incident_wave(system_conf, input_multem);
    psi_0 = output_incident_wave.psi_0;
    probe(:,:,n) = abs(psi_0).^2;
end
bessle_slice=squeeze(probe(512,:,:));
subplot(1,2,2)
imagesc(linspace(-10,10,1024),defocus,bessle_slice');
xlabel('distance (A)')
ylabel('Defocus (A)')
