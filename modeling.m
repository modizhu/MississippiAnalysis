% Modeling of Lake Processes
clear
clc
close all

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
T_dn = total_EC(:, 14);

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

total_EC_processed.LE(find(abs(total_EC_processed.LE)>1000)) = nan;
total_EC_processed.H(find(abs(total_EC_processed.H)>1000)) = nan;
Bowen = total_EC_processed.H ./ total_EC_processed.LE;

% Data from own EC codes
own_EC_data_time = [data13.Time; data14.Time; data15.Time; data16.Time];
own_EC_data_H = [data13.H_dn; data14.H_dn; data15.H_dn; data16.H_dn];
own_EC_data_E = [data13.E_dn; data14.E_dn; data15.E_dn; data16.E_dn];

%% Do The Modeling

Is = 1600;  % Water Thermal Inertia
z = 0.5;

Temp = total_Radiation(:,9) - 273.15;
RH = ones(1, length(Temp));
qs = Qs(Temp, RH);
[EMEP_dn, HMEP_dn, GMEP_dn, I0] = F_MEP_EHG(Rn, Temp, qs, Is, z, 1, 0);
GMEP_dn = GMEP_dn - (Rs_down-Rs_up);

figure
set(gcf,'Position',[200 100 1500 700])
subplot(2, 1, 1)
hold on;grid on;
plot(Time, EMEP_dn, 'Color','#fd79a8');
plot(Time, HMEP_dn, 'Color','#fdcb6e');
plot(Time, GMEP_dn, 'Color','#00cec9');
plot(Time, Rn, 'Color','#0984e3')
legend('MEPE','MEPH','MEPQ','Rn')
ylabel('(Wm^{-2})')
subplot(2, 1, 2)
hold on;grid on;
plot(Time, EMEP_dn, 'Color','#fd79a8');
plot(Time, HMEP_dn, 'Color','#fdcb6e');
plot(Time, GMEP_dn, 'Color','#00cec9');
plot(Time, Rn, 'Color','#0984e3')
ylabel('(Wm^{-2})')
xlim([datetime(2022, 9, 17), datetime(2022, 9, 24)])


total_EC_processed.LE(find(total_EC_processed.LE>600 | total_EC_processed.LE < -150)) = nan;
total_EC_processed.H(find(total_EC_processed.H>400 | total_EC_processed.H < -50)) = nan;

A_E = []; B_E = [];
for i = 1 : length(Time)
    if ismember(Time(i), Time_EC_processed)
        A_E = [A_E, EMEP_dn(i)];
        B_E = [B_E, total_EC_processed.LE(find(Time_EC_processed == Time(i)))];
    end
end

idx = [];
for i = 1 : length(A_E)
    if ~(isnan(A_E(i)) | isnan(B_E(i)))
        idx = [idx, i];
    end
end
P = polyfit(A_E(idx), B_E(idx), 1);
r_2_E = corrcoef(A_E(idx),B_E(idx)); r_2_E = r_2_E(1,2);
RMSE_E = sqrt(mean((A_E(idx) - B_E(idx)).^2));
NRMSE_E = goodnessOfFit(A_E(idx),B_E(idx),'NRMSE');

A_H = []; B_H = [];
for i = 1 : length(Time)
    if ismember(Time(i), Time_EC_processed)
        A_H = [A_H, HMEP_dn(i)];
        B_H = [B_H, total_EC_processed.H(find(Time_EC_processed == Time(i)))];
    end
end

idx = [];
for i = 1 : length(A_H)
    if ~(isnan(A_H(i)) | isnan(B_H(i)))
        idx = [idx, i];
    end
end
P = polyfit(A_H(idx), B_H(idx), 1);
r_2_H = corrcoef(A_H(idx),B_H(idx)); r_2_H = r_2_H(1,2);
RMSE_H = sqrt(mean((A_H(idx) - B_H(idx)).^2));
NRMSE_H = goodnessOfFit(A_H(idx),B_H(idx),'NRMSE');

figure
set(gcf,'Position',[200 100 1500 300])
subplot(2, 1, 1)
hold on; grid on;
plot(Time_EC_processed, A_E);
plot(Time_EC_processed, B_E);
ylabel('Wm^{-2}')
title('(a)')
legend('MEP','OBS')
subplot(2, 1, 2)
hold on; grid on;
plot(Time_EC_processed, A_H);
plot(Time_EC_processed, B_H);
ylabel('Wm^{-2}')
title('(b)')

mean_MEP_E = mean(A_E,'omitnan');
mean_MEP_H = mean(A_H,'omitnan');
mean_OBS_E = mean(B_E,'omitnan');
mean_OBS_H = mean(B_H,'omitnan');

% figure
% set(gcf,'Position',[200 100 1500 300])
% subplot(1, 2, 1)
% hold on; grid on;
% scatter(A_E, B_E)
% title('(a)')
% subplot(1, 2, 2)
% hold on; grid on;
% scatter(A_H, B_H)
% title('(b)')