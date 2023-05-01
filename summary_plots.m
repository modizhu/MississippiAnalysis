%% Make Plots of The concat data
clear
clc
close all
startup
% raditaion sensor data
data1 = load('RBR_flux_RMet_1Min_20220609_20220812.mat');
data2 = load('RBR_flux_RMet_1Min_20220812_20220913.mat');
data3 = load('RBR_flux_RMet_1Min_20220913_20221015.mat');
data4 = load('RBR_flux_RMet_1Min_20221015_20221122.mat');

% EC sensor data
data5 = load('total_1min_20220619_20220812.mat');
data6 = load('total_1min_20220812_20220913.mat');
data7 = load('total_1min_20220913_20221015.mat');
data8 = load('total_1min_20221015_20221122.mat');

% EC Processed data
data9 = readtable('eddypro_1_dn.csv');
data10 = readtable('eddypro_2_dn.csv');
data11 = readtable('eddypro_3_dn.csv');
data12 = readtable('eddypro_4_dn.csv');

% Modi's code dealing EC data
data13 = load('Energy_Flux_Result_1.mat');
data14 = load('Energy_Flux_Result_2.mat');
data15 = load('Energy_Flux_Result_3.mat');
data16 = load('Energy_Flux_Result_4.mat');
%% Parse the data
% Data from radiometer
Time = [data1.Time; data2.Time; data3.Time; data4.Time];
Rain = [data1.Rain; data2.Rain; data3.Rain; data4.Rain];
Rl_down = [data1.Rl_down; data2.Rl_down; data3.Rl_down; data4.Rl_down];
Rl_up = [data1.Rl_up; data2.Rl_up; data3.Rl_up; data4.Rl_up];
Rn = [data1.Rn; data2.Rn; data3.Rn; data4.Rn];
Rs_down = [data1.Rs_down; data2.Rs_down; data3.Rs_down; data4.Rs_down];
Rs_up = [data1.Rs_up; data2.Rs_up; data3.Rs_up; data4.Rs_up];
water_temp = [data1.water_temp; data2.water_temp; data3.water_temp; data4.water_temp];
windspeed = [data1.windspeed; data2.windspeed; data3.windspeed; data4.windspeed];
total_Radiation = [data1.real_data; data2.real_data; data3.real_data; data4.real_data];

RnS = Rs_down - Rs_up;
RnL = Rl_down - Rl_up; RnL(find(abs(RnL)>200)) = nan;


% "TIMESTAMP" (1),"RECORD"(2),"Rn_Avg"(3),"albedo_Avg"(4),"R_SW_in_Avg"(5),"R_SW_out_Avg"(6),"R_LW_in_Avg"(7),"R_LW_out_Avg"(8),"T_nr_Avg"(9),"R_LW_in_meas_Avg"(10),"R_LW_out_meas_Avg"(11),
% "mean_windspeed_up"(12),"mean_winddir_up"(13),"std_wind_dir_up"(14),"mean_windspeed_dn"(15),"mean_winddir_dn"(16),"std_wind_dir_dn"(17),"PAR_density_Avg"(18),"SBTempC_Avg"(19),"TargTempC_Avg"(20),
% "Water_T_1_Avg"(21),"Water_T_2_Avg"(22),"Water_T_3_Avg"(23),"Water_T_4_Avg"(24),"Water_T_5_Avg"(25),"Water_T_6_Avg"(26),"Water_T_7_Avg"(27),"Water_T_8_Avg"(28),"Water_T_9_Avg"(29),"Water_T_10_Avg"(30),
% "Water_T_11_Avg"(31),"Water_T_12_Avg"(32),"Water_T_13_Avg"(33),"Water_T_14_Avg"(34),"Rain_mm_Tot"(35),"Press_irga_up_Avg"(36),"Press_irga_dn_Avg"(37),"PTemp_C_Avg"(38),"BattV_Avg"(39)
total_Radiation = [data1.real_data; data2.real_data; data3.real_data; data4.real_data];
Albedo = total_Radiation(:, 6)./total_Radiation(:, 5);
water_temp = total_Radiation(:, 21:34);
water_temp(find(water_temp>40 | water_temp < 0)) = nan;

% AirT = total_Radiation(:, 38);

% Data from EC instruments before processing
%[RECORD(1), Ux_up(2), Uy_up(3), Uz_up(4), Ts_up(5), diag_csat_up(6), co2_up(7), h2o_up(8), Press_irga_up(9), diag_irga_up(10),
% Ux_dn(11), Uy_dn(12), Uz_dn(13), Ts_dn(14), co2_dn(15), h2o_dn(16), Press_irga_dn(17), diag_irga_dn(18)];
total_EC = [data5.total_data; data6.total_data; data7.total_data; data8.total_data];
total_EC_Time = [data5.total_time; data6.total_time; data7.total_time; data8.total_time];
% Trim prev part
% 3195 is that there is some data left by Prof. Liu in his prev
% observations
total_EC = total_EC(3195:end, :);
total_EC_Time = total_EC_Time(3195:end, :);
total_time = total_EC_Time(1):minutes(1):total_EC_Time(end);

A = find(total_EC_Time == datetime(2022, 8, 13, 5, 52, 0));
B = find(total_EC_Time == datetime(2022, 8, 31, 16, 0, 0));
C = find(total_EC_Time == datetime(2022, 10, 13, 7, 35, 0));
D = find(total_EC_Time == datetime(2022, 10, 22, 9, 29, 0));
CO2_dn = total_EC(:, 15); 
CO2_dn(find(CO2_dn > 1000 | CO2_dn < 0)) = nan;
CO2_dn(A:B) = nan; CO2_dn(C:D) = nan;
CO2_dn = CO2_dn * 29/(44*1.29);       % to umol/mol
H2O_dn = total_EC(:, 16); 
H2O_dn(find(H2O_dn > 300 | H2O_dn < 0)) = nan; 
H2O_dn(A:B) = nan; H2O_dn(C:D) = nan;
H2O_dn = H2O_dn * 29/(18*1.29);       % to mmol/mol
CO2_up = total_EC(:, 7);  
CO2_up(find(CO2_up > 1000 | CO2_up < 0)) = nan; 
CO2_up = CO2_up * 29/(44*1.29);       % to umol/mol
H2O_up = total_EC(:, 8);  
H2O_up(find(H2O_up > 300 | H2O_up < 0)) = nan;  
H2O_up = H2O_up * 29/(18*1.29);       % to mmol/mol


T_dn = total_EC(:, 14);

% j = 1;
% for i = 1 : length(total_time)
%     if ismember(total_time(i), total_EC_Time)
%         j = j + 1;
%         continue
%     else
%         total_EC = [total_EC(1:j-1, :); nan(1, size(total_EC, 2)); total_EC(j:end, :)];
%         j = j + 1;
%     end
% end
% toal_EC_Time = total_time;


% Data from Eddy Pro
total_EC_processed = [data9; data10; data11; data12];
Time_EC_processed = [];
for i = 1 : size(total_EC_processed, 1)
    curr = total_EC_processed.date(i);
    curr_hour = cell2mat(total_EC_processed.time(i)); curr_hour = curr_hour(1:2); curr_hour = str2num(curr_hour);
    curr_minute = cell2mat(total_EC_processed.time(i)); curr_minute = curr_minute(4:5); curr_minute = str2num(curr_minute);
    curr_time = curr + hours(curr_hour) + minutes(curr_minute);
    Time_EC_processed = [Time_EC_processed; curr_time];
end

% total_time = Time_EC_processed(1):minutes(30):Time_EC_processed(end);
% j = 1;
% sub_table = table(nan(1,1), NaT(1, 1));
% for k = 3: size(total_EC_processed, 2)
%     eval(strcat('sub_table.Var',num2str(k),' = [nan(1,1)];'));
% end
% sub_table.Properties.VariableNames = total_EC_processed.Properties.VariableNames;
% for i = 1 : length(total_time)
%     if ismember(total_time(i), Time_EC_processed)
%         j = j + 1;
%         continue
%     else
%         total_EC_processed = [total_EC_processed(1:j-1, :); sub_table; total_EC_processed(j:end, :)];
%         j = j + 1;
%     end
% end
% Time_EC_processed = total_time;

% total_EC_processed.LE(find(abs(total_EC_processed.LE)>1000)) = nan;
% total_EC_processed.H(find(abs(total_EC_processed.H)>1000)) = nan;
Bowen = total_EC_processed.H ./ total_EC_processed.LE;
mean_LE = mean(total_EC_processed.LE,'omitnan');
mean_H = mean(total_EC_processed.H,'omitnan');
mean_bowen = mean_H/mean_LE;


% Data from own EC codes
own_EC_data_time_30 = [data13.Time; data14.Time; data15.Time; data16.Time];
own_EC_data_H = [data13.H_dn; data14.H_dn; data15.H_dn; data16.H_dn];
own_EC_data_E = [data13.E_dn; data14.E_dn; data15.E_dn; data16.E_dn];


%% Make plots
%% Radiation sensor
figure
set(gcf,'Position',[200 400 1500 300])
hold on; grid on
Albedo_hourly_mean = get_diurnal_hourly_mean(Albedo, Time, 0, 60);
plot(linspace(0, 24, 1440), Albedo_hourly_mean)
xlim([7.5, 16.5])
ylim([0 0.2])
title('Albedo')
xlabel('Hour')

figure
set(gcf,'Position',[200 100 1500 800])
subplot(3, 1, 1)
hold on; grid on
plot(Time, Rn)
ylabel('W m^{-2}')
ylim([-200 1300])
title('(a)')
subplot(3, 1, 2)
hold on; grid on
plot(Time, RnS)
ylabel('W m^{-2}')
title('(b)')
subplot(3, 1, 3)
hold on; grid on
plot(Time, RnL)
title('(c)')
ylabel('W m^{-2}')

mean_Rn = mean(Rn, 'omitnan');
mean_RnS = mean(RnS, 'omitnan');
mean_RnL = mean(RnL, 'omitnan');

diurnal_hourly_mean_Rn = get_diurnal_hourly_mean(Rn, Time, 0, 60);

figure
set(gcf,'Position',[200 400 1500 300])
hold on; grid on
plot(total_EC_Time, T_dn-273.15)
ylabel('^oC')
ylim([0 50])

depth = [0, 5, 10, 15, 20, 27.5, 35, 42.5, 50, 57.5, 82.5, 97.5, 122.5, 147.5];
figure
set(gcf,'Position',[200 100 1200 800])
for i = 1 : size(water_temp, 2)
    subplot(5, 3, i)
    hold on; grid on
    plot(total_EC_Time, T_dn-273.15)
    plot(Time, water_temp(:,i));
    title(strcat('Layer ', num2str(depth(i)),  ' cm Water Temp vs AirT'))
    ylim([0 50])
    ylabel('^oC')
end
ylabel('^oC')

figure
set(gcf,'Position',[200 100 1200 800])
hold on; grid on;
mean_water_temp = []; min_water_temp = []; max_water_temp = [];
N = 1;
for i = 1 : size(water_temp, 2)
    sub_mean = mean(water_temp(N:end,i), 'omitnan');
    sub_min = sub_mean - min(water_temp(N:end,i));
    sub_max = max(water_temp(N:end,i)) - sub_mean;
    mean_water_temp = [mean_water_temp, sub_mean];
    min_water_temp = [min_water_temp, sub_min];
    max_water_temp = [max_water_temp, sub_max];
end
mean_water_temp(2) = nan; mean_water_temp(5) = nan;
errorbar(mean_water_temp, depth, min_water_temp, max_water_temp, 'horizontal', "-o","MarkerSize",10,...
    "MarkerEdgeColor","blue","MarkerFaceColor",'r', 'LineWidth', 2, 'CapSize', 12 );
plot([mean_water_temp(1), mean_water_temp(3:4), mean_water_temp(6:end)], [depth(1), depth(3:4), depth(6:end)], 'b', 'LineWidth',2);
set(gca, 'YDir','reverse')
xlabel('Water Temperature (^oC)')
ylabel('Depth (cm)')

figure
set(gcf,'Position',[200 400 1500 300])
hold on; grid on
plot(Time, windspeed)

%% EC sensor
figure
set(gcf,'Position',[200 100 2000 800])
subplot(1, 2, 1)
hold on; grid on
plot(total_EC_Time, total_EC(:,2), 'Color','#fd79a8');
plot(total_EC_Time, total_EC(:,3), 'Color','#fdcb6e');
plot(total_EC_Time, total_EC(:,4), 'Color','#d63031');
plot(total_EC_Time, total_EC(:,11), 'Color','#6c5ce7');
plot(total_EC_Time, total_EC(:,12), 'Color','#00cec9');
plot(total_EC_Time, total_EC(:,13), 'Color','#0984e3');
legend('Ux up', 'Uy up', 'Uz up', 'Ux dn', 'Uy dn', 'Uz dn')
ylim([-10 10])
ylabel('ms^{-1}')
title('Three Axis Wind Speed')
subplot(1, 2, 2)
hold on; grid on
plot(total_EC_Time(1:30:end), smooth(total_EC(1:30:end,2)), 'Color','#fd79a8');
plot(total_EC_Time(1:30:end), smooth(total_EC(1:30:end,3)), 'Color','#fdcb6e');
plot(total_EC_Time(1:30:end), smooth(total_EC(1:30:end,4)), 'Color','#d63031');
plot(total_EC_Time(1:30:end), smooth(total_EC(1:30:end,11)), 'Color','#6c5ce7');
plot(total_EC_Time(1:30:end), smooth(total_EC(1:30:end,12)), 'Color','#00cec9');
plot(total_EC_Time(1:30:end), smooth(total_EC(1:30:end,13)), 'Color','#0984e3');
legend('Ux up', 'Uy up', 'Uz up', 'Ux dn', 'Uy dn', 'Uz dn')
ylim([-10 10])
ylabel('ms^{-1}')
title('Smoothed Three Axis Wind Speed')

mean_u_up = mean(abs(total_EC(:,2)), 'omitnan'); mean_u_dn = mean(abs(total_EC(:,11)), 'omitnan');
mean_v_up = mean(abs(total_EC(:,3)), 'omitnan'); mean_v_dn = mean(abs(total_EC(:,12)), 'omitnan');
mean_w_up = mean(abs(total_EC(:,4)), 'omitnan'); mean_w_dn = mean(abs(total_EC(:,13)), 'omitnan');


% CO2 and H2O concentration
figure
set(gcf,'Position',[200 100 1500 300])
hold on; grid on;
plot(total_EC_Time, CO2_dn);
plot(total_EC_Time, CO2_up);
ylabel('ppm')
legend('2m', '6m')


mean_CO2_dn = mean(CO2_dn, 'omitnan');
mean_CO2_up = mean(CO2_up, 'omitnan');

figure
set(gcf,'Position',[200 100 1500 600])
subplot(2, 1, 1)
hold on; grid on;
plot(total_EC_Time, H2O_dn);
plot(total_EC_Time, H2O_up);
ylabel('mmol/mol')
ylim([0 230])
legend('2m', '6m')
title('Water Vapor Concentration at Ross Barnett Reservoir')
subplot(2, 1, 2)
hold on; grid on;
plot(total_EC_Time, H2O_dn);
plot(total_EC_Time, H2O_up);
ylabel('mmol/mol')
xlim([datetime(2022, 7, 23), datetime(2022, 8, 4)])
ylim([0 100])
title('(b)')
% legend('2m', '6m')
% subplot(2, 1, 2)
% hold on; grid on;
% plot(total_EC_Time, H2O_dn);
% plot(total_EC_Time, H2O_up);
% ylabel('mmol/mol')

mean_H2O_dn = mean(H2O_dn, 'omitnan');
mean_H2O_up = mean(H2O_up, 'omitnan');

diurnal_hourly_mean_CO2_dn = get_diurnal_hourly_mean(CO2_dn,  total_EC_Time, 0, 60);
diurnal_hourly_mean_CO2_up = get_diurnal_hourly_mean(CO2_up,  total_EC_Time, 0, 60);
diurnal_hourly_mean_H2O_dn = get_diurnal_hourly_mean(H2O_dn,  total_EC_Time, 0, 60);
diurnal_hourly_mean_H2O_up = get_diurnal_hourly_mean(H2O_up,  total_EC_Time, 0, 60);

figure
set(gcf,'Position',[200 100 1200 450])
subplot(1, 2, 1)
hold on; grid on;
scatter(diurnal_hourly_mean_Rn, diurnal_hourly_mean_CO2_dn);
ylabel('ppm')
xlabel('W m^{-2}')
title('(a)')
subplot(1, 2, 2)
hold on; grid on;
scatter(diurnal_hourly_mean_Rn, diurnal_hourly_mean_CO2_up);
ylabel('ppm')
xlabel('W m^{-2}')
title('(b)')

figure
set(gcf,'Position',[200 100 1600 800])
subplot(1, 3, 1)
hold on; grid on;
yyaxis left
plot(linspace(0, 24, 1440), diurnal_hourly_mean_CO2_dn)
ylabel('ppm')
title('(a)')
xticks(0:4:24);
xticklabels({'0','4','8','12','16','20','24'});
yyaxis right
plot(linspace(0, 24, 1440), diurnal_hourly_mean_CO2_up)
ylabel('ppm')
legend('2m','6m')
subplot(1, 2, 2)
hold on; grid on;
scatter(diurnal_hourly_mean_CO2_dn, diurnal_hourly_mean_CO2_up)
ylabel('CO2 6m (ppm)')
xlabel('CO2 2m (ppm)')
title('(b)')

%% EC Processed
figure
total_EC_processed.LE(find(total_EC_processed.LE>600 | total_EC_processed.LE < -150)) = nan;
total_EC_processed.H(find(total_EC_processed.H>400 | total_EC_processed.H < -50)) = nan;
set(gcf,'Position',[200 100 1500 300])
hold on; grid on;
plot(Time, Rn, 'Color','#fdcb6e')
plot(Time_EC_processed(1: find(Time_EC_processed == datetime(2022, 7, 14, 1, 23, 0))), total_EC_processed.LE(1: find(Time_EC_processed == datetime(2022, 7, 14, 1, 23, 0))), 'Color','#fd79a8');
plot(Time_EC_processed(1: find(Time_EC_processed == datetime(2022, 7, 14, 1, 23, 0))), total_EC_processed.H(1: find(Time_EC_processed == datetime(2022, 7, 14, 1, 23, 0))), 'Color','#6c5ce7');
plot(Time_EC_processed(find(Time_EC_processed == datetime(2022, 7, 17, 22, 23, 0)):end), total_EC_processed.LE(find(Time_EC_processed == datetime(2022, 7, 17, 22, 23, 0)):end), 'Color','#fd79a8');
plot(Time_EC_processed(find(Time_EC_processed == datetime(2022, 7, 17, 22, 23, 0)):end), total_EC_processed.H(find(Time_EC_processed == datetime(2022, 7, 17, 22, 23, 0)):end), 'Color','#6c5ce7');
xlim([Time_EC_processed(1), Time_EC_processed(end)])
legend('Rn', 'E','H')
ylabel('W m^{-2}')

figure
E_dn_hourly_mean = get_diurnal_hourly_mean(total_EC_processed.LE, Time_EC_processed, 0, 2);
H_dn_hourly_mean = get_diurnal_hourly_mean(total_EC_processed.H, Time_EC_processed, 0, 2);
hold on; grid on;
plot(0:0.5:23.5, E_dn_hourly_mean)
plot(0:0.5:23.5, H_dn_hourly_mean)
ylabel('W m^{-2}')
legend('E','H')
xticks(0:4:24);
xticklabels({'0','4','8','12','16','20','24'});

figure
set(gcf,'Position',[200 100 1500 250])
hold on; grid on;
plot(Time_EC_processed, total_EC_processed.wind_speed);
ylabel('ms^{-1}')

% Power Spectrum 
[pxx,f] = plomb(total_EC_processed.LE,Time_EC_processed,[],10,'power');
f = f*86400;
[pk,f0] = findpeaks(pxx,f,'MinPeakHeight',300);
figure
set(gcf,'Position',[200 100 1500 250])
plot(f,pxx,f0,pk,'o')
xlabel('Frequency (day^{-1})')
title('Power Spectrum and Prominent Peak of E')
hold on; grid on
xlim([0.5 5])
ylabel('W m^{-2}Hz^{-1}')

[pxx,f] = plomb(total_EC_processed.H,Time_EC_processed,[],10,'power');
f = f*86400;
[pk,f0] = findpeaks(pxx,f,'MinPeakHeight',10);
figure
set(gcf,'Position',[200 100 1500 250])
plot(f,pxx,f0,pk,'o')
xlabel('Frequency (day^{-1})')
title('Power Spectrum and Prominent Peak of H')
hold on; grid on
xlim([0.5 5])
ylabel('W m^{-2}Hz^{-1}')

figure
subplot(1, 2, 1)
hold on; grid on;
scatter(abs(total_EC_processed.w_unrot), total_EC_processed.LE)
subplot(1, 2, 2)
hold on; grid on;
scatter(abs(total_EC_processed.w_unrot), total_EC_processed.H)

%% Summary of Modi's code of processing EC E and H

figure
set(gcf,'Position',[200 100 1500 250])
hold on; grid on;
plot(Time, Rn, 'Color','b')
plot(own_EC_data_time_30(1: find(own_EC_data_time_30 == datetime(2022, 7, 14, 1, 0, 0))), own_EC_data_E(1: find(own_EC_data_time_30 == datetime(2022, 7, 14, 1, 0, 0))), 'Color','m');
plot(own_EC_data_time_30(1: find(own_EC_data_time_30 == datetime(2022, 7, 14, 1, 0, 0))), own_EC_data_H(1: find(own_EC_data_time_30 == datetime(2022, 7, 14, 1, 0, 0))), 'Color','r');
plot(own_EC_data_time_30(find(own_EC_data_time_30 == datetime(2022, 7, 17, 22, 30, 0)):end), own_EC_data_E(find(own_EC_data_time_30 == datetime(2022, 7, 17, 22, 30, 0)):end), 'Color','m');
plot(own_EC_data_time_30(find(own_EC_data_time_30 == datetime(2022, 7, 17, 22, 30, 0)):end), own_EC_data_H(find(own_EC_data_time_30 == datetime(2022, 7, 17, 22, 30, 0)):end), 'Color','r');
legend('Rn', 'E','H')
ylabel('W m^{-2}')
xlim([own_EC_data_time_30(1) own_EC_data_time_30(end)])

figure
mean_Rn_hourly = get_diurnal_hourly_mean(Rn, Time, 0, 60);
mean_own_E_dn_hourly = get_diurnal_hourly_mean(own_EC_data_E, own_EC_data_time_30, 0, 2);
mean_own_H_dn_hourly = get_diurnal_hourly_mean(own_EC_data_H, own_EC_data_time_30, 0, 2);
hold on; grid on;
yyaxis left
plot(linspace(0, 24, 1440), mean_Rn_hourly, 'b-')
plot(0:0.5:23.5, mean_own_E_dn_hourly, 'Color','#00cec9','LineStyle','-')
ylabel('R_n E W m^{-2}')
yyaxis right
plot(0:0.5:23.5, mean_own_H_dn_hourly)
ylabel('H W m^{-2}')
legend('Rn', 'E','H')
xticks(0:4:24);
xticklabels({'0','4','8','12','16','20','24'});
mean_own_E = mean(own_EC_data_E,'omitnan');
mean_own_H = mean(own_EC_data_H,'omitnan');
mean_own_bowen = mean_own_H/mean_own_E;

figure
set(gcf,'unit','normalized','position',[0.2,0.2,0.54,0.32]);
maxY1 = max(mean_Rn_hourly);
maxY2 = max(mean_own_E_dn_hourly);
maxY3 = max(mean_own_H_dn_hourly);
minY1 = min(mean_Rn_hourly);
minY2 = min(mean_own_E_dn_hourly);
minY3 = min(mean_own_H_dn_hourly);
newY2 = (mean_own_E_dn_hourly - minY2)/(maxY2 - minY2);   % 归一化
newY2 = newY2*(maxY1 - minY1) + minY1;  % 反归一化
newY3 = (mean_own_H_dn_hourly - minY3)/(maxY3 - minY3);
newY3 = newY3*(maxY1 - minY1) + minY1;
constY3 = (0 - minY3)/(maxY3 - minY3)*(maxY1 - minY1) + minY1;
h1 = axes('position', [0.1 0.1 0.5 0.8]); grid on; hold on;
set(h1, 'ycolor', 'b')
plot(linspace(0, 24, 1440), mean_Rn_hourly,'b');
plot(0:0.5:23.5, newY2,'m');
plot(0:0.5:23.5, newY3,'r');
ylim([minY1, maxY1])
xticks(0:4:24);
xticklabels({'0','4','8','12','16','20','24'});
ylabel('Rn (W m^{-2})')
legend('Rn','E','H')
h2 = axes('position', [0.63 0.1 0.005 0.8]);
plot(0:0.5:23.5, mean_own_E_dn_hourly, 'w')
set(h2, 'ycolor', 'm', 'yaxislocation', 'right', 'xtick', []);
hold on
limX2 = get(h2, 'Xlim');
limY2 = get(h2, 'Ylim');
plot([limX2(2), limX2(2)], limY2, 'Color','m');
ylim([minY2, maxY2])
hold off
box off
ylabel('E (W m^{-2})');
h3 = axes('position', [0.68 0.1 0.005 0.8]);
plot(0:0.5:23.5, mean_own_H_dn_hourly, 'w')
set(h3, 'ycolor', 'r', 'yaxislocation', 'right', 'xtick', [])
hold on
limX3 = get(h3, 'Xlim');
limY3 = get(h3, 'Ylim');
plot([limX3(2), limX3(2)], limY3, 'Color', 'r');
hold off
box off
ylim([minY3, maxY3])
ylabel('H (W m^{-2})');
set(gcf,'color','white'); 


% Water Temp and Air Temp and Solar Radiation Down
figure
set(gcf,'Position',[200 100 1200 850])
subplot(3, 1, 1)
hold on; grid on;
yyaxis left
plot(total_EC_Time(10308:end), T_dn(10308:end)-273.15, 'Color','#00cec9')
plot(Time(58381:end), water_temp(58381:end, 1),'b-');
ylabel('^oC')
yyaxis right
plot(Time(58381:end), Rn(58381:end))
ylabel('(W m^{-2})')
legend('TA','TS', 'Rn')
title('(a)')
subplot(3, 1, 2)
hold on; grid on;
ylabel('(W m^{-2})')
yyaxis left
mean_AirT_hourly = get_diurnal_hourly_mean(T_dn(10308:end) - 273.15, total_EC_Time(10308:end), 0, 60);
mean_SurfaceWaterTemp_hourly = get_diurnal_hourly_mean(water_temp(58381:end, 1), Time(58381:end), 0, 60);
plot(linspace(0, 24, 1440), mean_AirT_hourly, 'Color','#00cec9');
plot(linspace(0, 24, 1440), mean_SurfaceWaterTemp_hourly, 'b-');
ylabel('^oC')
yyaxis right
plot(linspace(0, 24, 1440), mean_Rn_hourly);
legend('TA','TS', 'Rn')
ylabel('(W m^{-2})')
xticks(0:4:24);
xticklabels({'0','4','8','12','16','20','24'});
title('(b)')
subplot(3, 1, 3)
hold on; grid on;
yyaxis left
plot(linspace(0, 24, 1440), mean_SurfaceWaterTemp_hourly - mean_AirT_hourly);
ylabel('^oC')
yyaxis right
plot(0:0.5:23.5, mean_own_H_dn_hourly)
legend('TS - TA', 'H')
ylabel('(W m^{-2})')
title('(c)')
xticks(0:4:24);
xticklabels({'0','4','8','12','16','20','24'});

mean_RnL = mean(RnL, 'omitnan');
mean_RnS = mean(RnS, 'omitnan');
mean_Q = mean_RnL - mean_own_E - mean_own_H;
mean_Rs_down = mean(Rs_down, 'omitnan');
mean_Rs_up = mean(Rs_up, 'omitnan');

save Mississippi_OBS.mat