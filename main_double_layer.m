% Codded by Refet Ali YALCIN. You can change and distribute the code
% please keep this part. The codes come WITHOUT ANY WARRANTY 
% In case of use, following articles that uses the code can be cited:
% https://doi.org/10.1088/2053-1591/ab28b8 Improving photosynthetic efficiency using greenhouse coatings with scattering and fluorescent pigments
% https://doi.org/10.1016/j.biosystemseng.2020.02.007 Improving crop production in solar illuminated vertical farms using fluorescence coatings
clear all
close all
clc

start_wl=300; % starting wavelength in nm, must be an integer
end_wl=2500; % last wavelength of the area of interest, must be an integer

repeat_no=10000;%100000; % # of montecarlo simulations for each wavelength
h1=300*10^-6; %thickness of coating in meters
h2 = 300*10^-6;
radius1=5.5*10^-6;%5.5*10^-6; % radius of fluorescent particles in meters
radius2=250*10^-9;
%sigma=18.5*10^-9;
f_v1=0.05;  % volume fraction of phosphor particles
f_v2=0.15;
polar_angle_degree=0; %polar angle of incident in degree 

polar_angle_rad=polar_angle_degree*pi/180;
wl=(start_wl:end_wl)';%(start_wl:(end_wl-start_wl)/220:end_wl)';%for 1 nm interval, need to modify monte carlo code
number_wl=length(wl); 
lamda=wl*10^-9;

QY=0.9*ones(number_wl,1);%QY_l(lamda); % %quantum yield
n_medium1=real(nk_PDMS(lamda));%real(nk_PDMS(lamda));%PMMA_n(lamda); %
%k_medium1=imag(nk_PDMS(lamda));%PMMA_k(lamda);
k_medium1=zeros(number_wl,1); % enable for non absorbing medium case

n_medium2=real(nk_PDMS(lamda));%real(nk_PDMS(lamda));%PMMA_n(lamda); %
%k_medium=imag(nk_PDMS(lamda));%PMMA_k(lamda);
k_medium2=zeros(number_wl,1); % enable for non absorbing medium case

n_subs=ones(number_wl,1);%real(nk_Ag_new(lamda));% %substrate is air
k_subs=zeros(number_wl,1);%imag(nk_Ag_new(lamda));%

n_phosphor=real(nk_BAM(lamda));%real(nk_BAM(lamda));%real(nk_TiO2_new(lamda));%real(nk_BAM(lamda));%real(nk_coumarin1(lamda));%%n_yagce(lamda);
k_phosphor=imag(nk_BAM(lamda));%imag(nk_TiO2_new(lamda));%imag(nk_BAM(lamda));%imag(nk_coumarin1(lamda));%k_yagce(lamda);

n_particle=real(nk_TiO2_new(lamda));%real(nk_TiO2_new(lamda));%real(nk_BAM(lamda));%real(nk_coumarin1(lamda));%%n_yagce(lamda);
k_particle=imag(nk_TiO2_new(lamda));%imag(nk_TiO2_new(lamda));%imag(nk_BAM(lamda));%imag(nk_coumarin1(lamda));%k_yagce(lamda);
%n_Ag = nk_Ag_new(lamda);

%pre_process % call pre process to calculate coefficients
pre_process_double_layer

% Below part calculates reflectance and refraction at the air - medium1 interface
cos_teta_prime=zeros(length(lamda),1); %the cos of the ray after refracted from air to medium1
sur_reflection=zeros(length(lamda),1); %surface reflection at air-medium1 interfece

for i=1:length(lamda)
    cos_teta_prime(i)=cos(F_fresnel_2(n_medium1(i),k_medium1(i),polar_angle_rad));
    cos_teta=cos(polar_angle_rad);
    sin_teta=sqrt(1-cos_teta*cos_teta);
    carpan2=1/(n_medium1(i)-1i*k_medium1(i));
    sin_x2=sin_teta*carpan2;
    cos_x2=sqrt(1-sin_x2*sin_x2);
    carpan1=cos_teta/cos_x2;
    carpan3=cos_x2/cos_teta;
    E_parallel=(carpan1-carpan2)/(carpan1+carpan2);
    R_parallel=E_parallel*conj(E_parallel);
    E_orth=(carpan3-carpan2)/(carpan3+carpan2);
    R_orth=E_orth*conj(E_orth);
    reflectance=real(R_parallel+R_orth)*0.5;
    sur_reflection(i)=reflectance;
end
% end of refraction and surface reflection calculation

% set database for photons
db_absorption_no=zeros(number_wl,1);
db_reflect_no=zeros(number_wl,number_wl);
db_trans_no=zeros(number_wl,number_wl);
tic

for k=start_wl:end_wl
    absorption_no=0;
    reflect_no=zeros(number_wl,1);
    trans_no=zeros(number_wl,1);
    wl_index=k-start_wl+1;
    for i=1:repeat_no
        [absorption_no_new,reflect_no_new,trans_no_new] = monte_carlo_double_layer(h1, h2,k, scat_prob1, scat_prob2, ext_tot1, ext_tot2, g1, g2, QY_modified, start_wl, number_wl, inv_cdf, cos_teta_prime(wl_index), sur_reflection(wl_index), n_medium1, k_medium1, n_subs, k_subs);%monte_carlo_both_layer_fluorescence(h1, h2, k, scat_prob1, scat_prob2, ext_tot1, ext_tot2, g1, g2, QY_modified1, QY_modified2, start_wl, number_wl, inv_cdf1, inv_cdf2, cos_teta_prime(wl_index), sur_reflection(wl_index), n_medium1, k_medium1, n_subs, k_subs);
        absorption_no=absorption_no + absorption_no_new;
        reflect_no=reflect_no + reflect_no_new;
        trans_no=trans_no + trans_no_new;
    end
    db_reflect_no(:,wl_index)=reflect_no;
    db_absorption_no(wl_index)=absorption_no;
    db_trans_no(:,wl_index)=trans_no;
    
    clc
    disp([num2str(floor(wl_index*100/number_wl)),'% has been completed.']);
end
toc
disp(['100% has been completed.']);
post_process
