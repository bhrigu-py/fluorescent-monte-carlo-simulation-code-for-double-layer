notOctave = exist('OCTAVE_VERSION', 'builtin') == 0;


beta_non_fl=zeros(number_wl,1);
beta_fl=zeros(number_wl,1);
%abs and scat coef

r1=radius1;

%
mu_tot_arr_sigma1=zeros(length(lamda),1);
beta_sigma1=zeros(length(lamda),1);
alfa_sigma1=zeros(length(lamda),1);
scat_prob_arr_sigma1=zeros(length(lamda),1);
g_arr_sigma1=zeros(length(lamda),1);
QY_modified_sigma=zeros(length(lamda),1);

    Area1=pi*r1^2;
    V1=4*pi*r1^3/3;
    for i=1:number_wl
        x1=2*pi*r1*n_medium1(i)/lamda(i);
        m1=(n_phosphor(i)+1i*k_phosphor(i))/n_medium1(i);
        fonksiyon1=Mie(m1,x1);
        Qabs1=fonksiyon1(3);
        Qsca1=fonksiyon1(2);
        g_arr_sigma1(i)=fonksiyon1(5);
        Qsca1=Qsca1*(1.0-g_arr_sigma1(i));
        alfa_sigma1(i)=f_v1*Qsca1*Area1/V1;
        beta_fl=f_v1*Qabs1*Area1/V1;
        beta_non_fl=(1-f_v1)*4*pi*k_medium1(i)/lamda(i); %absorption by medium (non fluorescent)
        beta_sigma1(i)=beta_fl+beta_non_fl;
        QY_modified_sigma(i)=QY(i)*beta_fl/beta_sigma1(i);%QY_modified_sigma(i,z)=QY(i)*beta_fl/beta_sigma(i,z); %Modified version of QY is used since the probability of re-emitting the absorbed rays by non fluorescent part is zero.
        mu_tot_arr_sigma1(i)=alfa_sigma1(i)+beta_sigma1(i);
        scat_prob_arr_sigma1(i)=alfa_sigma1(i)/mu_tot_arr_sigma1(i);
    end
%end
ext_tot1=mu_tot_arr_sigma1;
scat_prob1=scat_prob_arr_sigma1;
QY_modified=QY_modified_sigma;
g1=g_arr_sigma1;
beta1=beta_sigma1;
alfa1=alfa_sigma1;

%%
r2=radius2;
mu_tot_arr_sigma2=zeros(length(lamda),1);
beta_sigma2=zeros(length(lamda),1);
alfa_sigma2=zeros(length(lamda),1);
scat_prob_arr_sigma2=zeros(length(lamda),1);
g_arr_sigma2=zeros(length(lamda),1);

    Area2=pi*r2^2;
    V2=4*pi*r2^3/3;
    for i=1:number_wl
        x2=2*pi*r2*n_medium2(i)/lamda(i);%2*pi*r_vector(z)*n_medium(i)/lamda(i);
        m2=(n_particle(i)+1i*k_particle(i))/n_medium2(i);
        fonksiyon2=Mie(m2,x2);
        Qabs2=fonksiyon2(3);
        Qsca2=fonksiyon2(2);
        g_arr_sigma2(i)=fonksiyon2(5);
        Qsca2=Qsca2*(1.0-g_arr_sigma2(i));
        alfa_sigma2(i)=f_v2*Qsca2*Area2/V2;
        beta_particle=f_v2*Qabs2*Area2/V2;
        beta_matrix=(1-f_v2)*4*pi*k_medium2(i)/lamda(i); 
        beta_sigma2(i)=beta_particle+beta_matrix;
        mu_tot_arr_sigma2(i)=alfa_sigma2(i)+beta_sigma2(i);
        scat_prob_arr_sigma2(i)=alfa_sigma2(i)/mu_tot_arr_sigma2(i);
    end
%end
ext_tot2=mu_tot_arr_sigma2;
scat_prob2=scat_prob_arr_sigma2;
g2=g_arr_sigma2;
beta2=beta_sigma2;
alfa2=alfa_sigma2;
%QY_modified2=zeros(length(lamda),1);%0.9*ones(number_wl,1);
%%


% pdf ve cdf

data_flo=flo_emission_data_BAM();%flo_emission_data();
flo_wl_start=data_flo(1,1);
flo_wl_end=data_flo(end,1);

data_1=data_flo(:,1);
data_2=data_flo(:,2);
wl_arr=(flo_wl_start:flo_wl_end)';
result_arr=interp1(data_1,data_2,wl_arr);
wave_flo=flo_wl_start:flo_wl_end;
wave_flo_no=length(wave_flo);
x=zeros(wave_flo_no,1);
y=zeros(wave_flo_no,1);
cdf=zeros(wave_flo_no,1);

for i=1:wave_flo_no
   x(i)=(i-1)/(wave_flo_no-1);
   y(i)=result_arr(i);
end
y=y/trapz(x,y);
for i=2:wave_flo_no
    cdf(i)=trapz(x(1:i),y(1:i));
end

cdf_fl_new=zeros(10001,1);
inv_cdf=zeros(10001,1);
for i=1:10001
   cdf_fl_new(i)=i-1;
   inv_cdf(i)=interp1(cdf,wave_flo,cdf_fl_new(i)/10000);
end
% inv_cdf2=zeros(10001,1);
% for i=1:10001
%    cdf_fl_new(i)=i-1;
%    inv_cdf2(i)=interp1(cdf,wave_flo,cdf_fl_new(i)/10000);
% end

fig1=figure(1);

plot(wl,exc_BAM(wl),'-k',data_1,data_2/max(data_2),'--k','LineWidth',2)%plot(wl,exc_yagce(wl),'-k',data_1,data_2/max(data_2),'--k','LineWidth',2)
ylabel('Normalized Intensity [a.u.]')
xlabel('Wavelength [nm]')
xlim([300 800])
legend('Excitation','Emission','Location','NorthEast')
% saveas(fig1,'abs_emis.fig')
% saveas(fig1,'abs_emis.emf')

% refractive indices
% if notOctave %yyaxis problem in octave
%     fig2=figure(2);
%     set(fig2,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
%     hold on
%     yyaxis left
%     plot(wl,n_phosphor,'-k','LineWidth',2)
%     ylabel('n')
%     ylim([1.5 3])
%     yyaxis right
%     plot(wl,k_phosphor,'--k','LineWidth',2)
%     ylabel('k')
%     hold off
%     box on
%     xlabel('Wavelength [nm]')
%     xlim([300 end_wl])
%     legend('n','k','Location','SouthEast')
% %     saveas(fig2,'ref_ind.fig')
% %     saveas(fig2,'ref_ind.emf')

%     absorption and scat coeff
    fig3=figure(3);
    set(fig3,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
    hold on
    yyaxis left
    plot(wl,alfa1,'-k',wl,beta1,':k',wl,alfa2,'k',wl,beta2,'--k','LineWidth',2)
    ylabel('Scattering and Absorption Coefficient [1/m]')
%     ylim([0  1.1*max(alfa)])
    %yyaxis right
%     plot(wl,g,'--k','LineWidth',2)
%     ylabel('Asymmetry Parameter (g)')
%     ylim([-1 1])
%     hold off
    box on
    xlabel('Wavelength [nm]')
    xlim([300 end_wl])
    legend('Scattering Coefficient_BAM','Absorption Coefficient_BAM','Scattering Coefficient_TiO2','Absorption Coefficient_TiO2','Location','SouthEast')
% %     saveas(fig3,'scat_abs_coef.fig')
% %     saveas(fig3,'scat_abs_coef.emf')

%     CDF and PDF
    fig5=figure(5);
    set(fig5,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
    hold on
    yyaxis left
    plot(flo_wl_start:flo_wl_end,y/trapz(flo_wl_start:flo_wl_end,y),'-k','LineWidth',2);
    xlabel('Wavelength [nm]')
    ylabel('PDF [1/nm]')
    ylim([0 0.02])
    xlim([flo_wl_start flo_wl_end])
    yyaxis right
    plot(flo_wl_start:flo_wl_end,cdf,'--k','LineWidth',2);
    xlabel('Wavelength [nm]')
    ylabel('CDF')
    ylim([0 1])
    xlim([flo_wl_start flo_wl_end])
    legend('Probability Density Function','Cumulative Distribution Function','Location','northwest')
    hold off
    box on
%     saveas(fig5,'pdf_cdf.fig')
%     saveas(fig5,'pdf_cdf.emf')

%inv cdf
% fig6=figure(6);
% plot(linspace(0,1,length(inv_cdf)),inv_cdf,'-k','LineWidth',2);
% ylabel('Wavelength [nm]')    
% xlabel('Random Number')    
% xlim([0 1])
% ylim([flo_wl_start flo_wl_end])
% legend('Inverse of Cumulative Distribution Function','Location','southeast')
% saveas(fig6,'invcdf.fig')
% saveas(fig6,'invcdf.emf')
    