IMU_data = IMU_read('long_walk_1.h5');
dt_list = diff(IMU_data.time);

side_list = ["Right" , "Left"];
component_list = ["SupInf", "LatMed", "PostAnt"]; %directions identified from static data

results = struct;

for s=1:2
    side_str = side_list(s);
    results.(side_str) = struct;
    disp(s);
    for c=1:3
        disp(c);
        component_idx = c;
        results.(side_str).(component_list(component_idx)) = struct;
        activity = segment_extraction(IMU_data, component_idx, side_str);
        activity = state_space_analysis(activity, dt_list, component_idx, side_str);
        LE_results = plot_results(activity, component_idx, side_str);
        results.(side_str).(component_list(component_idx)) = LE_results;
    end
end

plot_BPresults(results, component_list);


function plot_BPresults(results, component_list)
    path_img_final = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\IMU sample data\Images\";
    figure('units','normalized','outerposition',[0 0 1 1]);
    for c=1:3
        component_idx = c;
        subplot(2,3,c);
        boxplot([results.Left.(component_list(component_idx)).s, results.Right.(component_list(component_idx)).s]);
        if c==1
            ylabel("\lambda_s");
        end
        ylim([100, 250]);
        title(component_list(component_idx));
        ax = gca;
        ax.FontSize = 20;
        set(gca,'XTickLabel',{' '});
        set(findobj(gca,'type','line'),'linew',1);
    end

    for c=4:6
        component_idx = c-3;
        subplot(2,3,c);
        boxplot([results.Left.(component_list(component_idx)).b, results.Right.(component_list(component_idx)).b], ["Left", "Right"]);
        if c==4
            ylabel("\lambda_b");
        end
        ylim([-4,2]);
        ax = gca;
        ax.FontSize = 20;
        set(findobj(gca,'type','line'),'linew',1);
    end

    filename_BPresults = path_img_final + "BPresults.png";
    exportgraphics(gcf, filename_BPresults);
    close(gcf);
end

function LE_results = plot_results(activity, component_idx, side_str)
    path_img_final = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\IMU sample data\Images\" + side_str + "\" + string(component_idx) + "\";
    NumActivities = length(activity);
    
    LE_results = struct;
    LE_results.s = nan(NumActivities,1);
    LE_results.b = nan(NumActivities,1);
    
    figure('units','normalized','outerposition',[0 0 1 1]);

    subplot(2,1,1);
    for k=1:NumActivities
        LE_results.s(k) = activity(k).LE_s;
        scatter(k, activity(k).LE_s, 50, 'b', 'filled');
        hold on;
    end
    xlabel("Trial number");
    ylabel("\lambda_s");
    title("Results for \lambda_s");
    ax = gca;
    ax.FontSize = 20;

    subplot(2,1,2);
    for k=1:NumActivities
        LE_results.b(k) = activity(k).LE_b;
        scatter(k, activity(k).LE_b, 50, 'r', 'filled');
        hold on;
    end
    xlabel("Trial number");
    ylabel("\lambda_b");
    title("Results for \lambda_b");
    ax = gca;
    ax.FontSize = 20;

    filename_results = path_img_final + "Results.png";
    exportgraphics(gcf, filename_results);
    close(gcf);
end

function activity = state_space_analysis(activity, dt_list, component_idx, side_str)
    fs = floor(1./dt_list(1));
    path_img_final = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\IMU sample data\Images\" + side_str + "\" + string(component_idx) + "\";
    NumActivities = length(activity);
    for k=1:NumActivities
        disp(k);
        [XR_j,eLag_j,eDim_j] = phaseSpaceReconstruction(activity(k).data_j_interp);
        activity(k).XR_j = XR_j;
        activity(k).eLag_j = eLag_j;
        activity(k).eDim_j = eDim_j;

        [~,estep,ldiv] = lyapunovExponent(XR_j,fs,eLag_j,eDim_j,'ExpansionRange',[1,1000]);

        stride_interval = mean(activity(k).stride_time);
        stride_deltaN_unit = floor(stride_interval/dt_list(1));
        estride = estep/stride_deltaN_unit;
        estride_firstTen_bool = estride<11;
        estride_firstTen = estride(estride_firstTen_bool);
        ldiv_firstTen = ldiv(estride_firstTen_bool);

        estride_small_bool = estride < 0.5;
        estride_small = estride(estride_small_bool);
        ldiv_small = ldiv(estride_small_bool);
        p_small = polyfit(estride_small, ldiv_small, 1);
        f_small = polyval(p_small, estride_small);

        estride_big_bool = (4 < estride ) & (estride < 10);
        estride_big = estride(estride_big_bool);
        ldiv_big = ldiv(estride_big_bool);
        p_big = polyfit(estride_big, ldiv_big, 1);
        f_big = polyval(p_big, estride_big);

        activity(k).LE_s = p_small(1);
        activity(k).LE_b = p_big(1);

        figure('units','normalized','outerposition',[0 0 1 1]);
        plot(estride_firstTen, ldiv_firstTen, 'k', 'LineWidth', 1, 'DisplayName', 'Logarithmic Divergence');
        hold on;
        plot(estride_small, f_small, 'r', 'LineWidth', 1.5, 'DisplayName', '\lambda_s line');
        plot(estride_big, f_big, 'b', 'LineWidth', 1.5, 'DisplayName', '\lambda_b line');
        xlabel("Number of Strides");
        ylabel("Logarithmic Divergence");
        title("\lambda_s = " + string(round(p_small(1),2)) + " \lambda_L = " + string(round(p_big(1),2)));
        grid;
        legend('location', 'best');
        ax = gca;
        ax.FontSize = 20;

        filename_results = path_img_final + "LE_plot_" + string(k) + ".png";
        exportgraphics(gcf, filename_results);
        close(gcf);
    end

end

function activity = segment_extraction(IMU_data, component_idx, side_str)
    dt_list = diff(IMU_data.time);
    fs = floor(1./dt_list(1));
    
    segment_str = side_str + "_Shin";
    
    data_j_raw = IMU_data.(segment_str).a(:, component_idx);
    frames = 1:length(data_j_raw);
    
    N = 1;
    maxDistance = 1;
    xyPoints = horzcat(frames', data_j_raw);
    [P, inlierIdx] = fitPolynomialRANSAC(xyPoints,N,maxDistance);
    RecoveredCurve = polyval(P,frames');
    outlier_frames = frames(~inlierIdx);
    outlier_data = data_j_raw(~inlierIdx);

    threshold_break_minutes = 1;
    diff_outlier_frames_minutes = diff(outlier_frames)/(fs*60);
    [~, locs] = findpeaks(diff_outlier_frames_minutes,'MinPeakHeight',threshold_break_minutes);
    activity_start_idx = outlier_frames(locs(1:end-1)+1);
    activity_end_idx = outlier_frames(locs(2:end));

    activity = struct;
    NumActivities = length(activity_start_idx);
    LyapunovExponent_list = nan(NumActivities, 1);
    for k=1:NumActivities
        activity(k).trial_idx = k;
        activity(k).start_idx = activity_start_idx(k);
        activity(k).end_idx = activity_end_idx(k);
        [stride_start_idx, stride_end_idx] = extract_stride(IMU_data, activity_start_idx(k), activity_end_idx(k));
        activity(k).stride_start_idx = activity_start_idx(k) + stride_start_idx;
        activity(k).stride_end_idx = activity_start_idx(k) + stride_end_idx;
        activity(k).num_strides = length(stride_start_idx);
        activity(k).stride_centered_start_idx = activity(k).stride_start_idx(floor(activity(k).num_strides/2)-75:floor(activity(k).num_strides/2)+75);
        activity(k).stride_centered_end_idx = activity(k).stride_end_idx(floor(activity(k).num_strides/2)-75:floor(activity(k).num_strides/2)+75);
        activity(k).num_centered_strides = length(activity(k).stride_centered_start_idx);
        activity(k).data_j_interp = [];
        activity(k).data_j_interp_frames = [];
        activity(k).stride_time = [];
        for s=1:activity(k).num_centered_strides
            x = activity(k).stride_centered_start_idx(s):activity(k).stride_centered_end_idx(s);
            xq = linspace(activity(k).stride_centered_start_idx(s),activity(k).stride_centered_end_idx(s));
            data_s = data_j_raw(x);
            data_s_interp = interp1(x, data_s, xq)';
            activity(k).data_j_interp = [activity(k).data_j_interp; data_s_interp];
            activity(k).data_j_interp_frames = [activity(k).data_j_interp_frames; xq'];
            activity(k).stride_time = [activity(k).stride_time; (activity(1).stride_centered_end_idx(2) - activity(1).stride_centered_start_idx(2))*dt_list(1)];
        end
        activity(k).data_j = data_j_raw(activity(k).stride_centered_start_idx(1):activity(k).stride_centered_end_idx(end));
        activity(k).data_j_frames = transpose(activity(k).stride_centered_start_idx(1):activity(k).stride_centered_end_idx(end));
        activity(k).LE_s = LyapunovExponent_list(k);
        activity(k).LE_b = LyapunovExponent_list(k);

    end

    path_img_final = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\IMU sample data\Images\" + side_str + "\" + string(component_idx) + "\";

    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(frames, data_j_raw, 'b', 'DisplayName', 'Raw Data');
    hold on;
    plot(frames, RecoveredCurve, 'r', 'LineWidth', 2, 'DisplayName', 'Static Line');
    start_plot = arrayfun(@(a)xline(a, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off'), activity_start_idx);
    end_plot = arrayfun(@(a)xline(a, 'k', 'LineWidth', 1, 'HandleVisibility', 'off'), activity_end_idx);
    xlabel("Frame");
    ylabel("Movement Data");
    title("Extracting Movement Trials");
    legend;
    ax = gca;
    ax.FontSize = 20;

    filename_extractTrial = path_img_final + "ExtractTrial_plot_" + string(component_idx) + ".png";
    exportgraphics(gcf, filename_extractTrial);
    close(gcf);
end

function [stride_start_idx, stride_end_idx] = extract_stride(IMU_data, trial_start_idx, trial_end_idx)
    acc_mag = vecnorm(IMU_data.Right_Shin.a,2,2);
    acc_mag_mvt = acc_mag(trial_start_idx:trial_end_idx);
    frames = 1:length(acc_mag_mvt);

    peak_separation_threshold = 0.01*(10^4);
    [~, locs] = findpeaks(acc_mag_mvt, 'MinPeakDistance', peak_separation_threshold);
    %h = histogram(diff(locs));

    stride_start_idx = frames(locs(1:end-1));
    stride_end_idx = frames(locs(2:end));

    %figure
    %plot(frames, acc_mag_mvt, 'b');
    %hold on;
    %start_plot = arrayfun(@(a)xline(a, 'k--', 'LineWidth', 1), stride_start_idx(100:110));
    %end_plot = arrayfun(@(a)xline(a, 'k--', 'LineWidth', 1), stride_end_idx(100:110));
end

function [IMU] = IMU_read(filepath)
filename = convertStringsToChars(filepath);

data_info1 = h5info(filename,'/Sensors');
data_info2 = h5info(filename,'/Processed');

num_sensors = length(data_info1.Groups);

button = h5read(filename,'/Annotations');
IMU.button.push_time = button.Time;
IMU.button.sensor = button.SensorID;
IMU.button.label = cellstr(button.Annotation');

sensor_names = cell(num_sensors,1); 
vectorlengths = zeros(num_sensors,1);

sensor_data_test = data_info1.Groups(1).Name;
sensor_data_size = length(h5read(filename,[sensor_data_test,'/Time']));
sensor_data_empty = sensor_data_size == 0;
if sensor_data_empty
    IMU = struct;
    return
end


for ii = 1:num_sensors
    sensor_data_loc1 = data_info1.Groups(ii).Name;
    sensor_data_loc2 = data_info2.Groups(ii).Name;
    
    sensor_info = h5info(filename,[sensor_data_loc1,'/Configuration']);
    sensor_label = sensor_info.Attributes(1).Value;
    sensor_label = strrep(sensor_label,' ','_');
    
    sensor_names(ii) = {sensor_label};
    IMU.(sensor_label).a = h5read(filename,[sensor_data_loc1,'/Accelerometer'])';
    IMU.(sensor_label).w = h5read(filename,[sensor_data_loc1,'/Gyroscope'])';
    IMU.(sensor_label).m = h5read(filename,[sensor_data_loc1,'/Magnetometer'])';
    IMU.(sensor_label).q = h5read(filename,[sensor_data_loc2,'/Orientation'])';
    IMU.(sensor_label).b = h5read(filename,[sensor_data_loc1,'/Barometer']);
    IMU.(sensor_label).temp = h5read(filename,[sensor_data_loc1,'/Temperature']);
    IMU.(sensor_label).time = h5read(filename,[sensor_data_loc1,'/Time']);
    
    vectorlengths(ii) = length(IMU.(sensor_label).time);
end

[maxlength, indmax] = max(vectorlengths);

for jj = 1:num_sensors
    if vectorlengths(jj) == maxlength
        
        if length(IMU.(char(sensor_names(jj))).a) == length(IMU.(char(sensor_names(jj))).time)
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'a');
        end
        if length(IMU.(char(sensor_names(jj))).w) == length(IMU.(char(sensor_names(jj))).time)
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'w');
        end
        if length(IMU.(char(sensor_names(jj))).m) == length(IMU.(char(sensor_names(jj))).time)
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'m');
        end
        if length(IMU.(char(sensor_names(jj))).b) == length(IMU.(char(sensor_names(jj))).time)
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'b');
        end
        if length(IMU.(char(sensor_names(jj))).temp) == length(IMU.(char(sensor_names(jj))).time)      
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'temp');
        end
        
    elseif vectorlengths(jj) ~= maxlength
        index = 1:maxlength;
        
        if length(IMU.(char(sensor_names(jj))).a) == length(IMU.(char(sensor_names(jj))).time)
            IMU.(char(sensor_names(jj))).a = interp1(IMU.(char(sensor_names(jj))).time,IMU.(char(sensor_names(jj))).a,IMU.(char(sensor_names(indmax))).time,'linear','extrap');
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'a');
        end
        if length(IMU.(char(sensor_names(jj))).w) == length(IMU.(char(sensor_names(jj))).time)
            IMU.(char(sensor_names(jj))).w = interp1(IMU.(char(sensor_names(jj))).time,IMU.(char(sensor_names(jj))).w,IMU.(char(sensor_names(indmax))).time,'linear','extrap');
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'w');
        end
        if length(IMU.(char(sensor_names(jj))).m) == length(IMU.(char(sensor_names(jj))).time)
            IMU.(char(sensor_names(jj))).m = interp1(IMU.(char(sensor_names(jj))).time,IMU.(char(sensor_names(jj))).m,IMU.(char(sensor_names(indmax))).time,'linear','extrap');
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'m');
        end
        if length(IMU.(char(sensor_names(jj))).b) == length(IMU.(char(sensor_names(jj))).time)
            IMU.(char(sensor_names(jj))).b = interp1(IMU.(char(sensor_names(jj))).time,IMU.(char(sensor_names(jj))).b,IMU.(char(sensor_names(indmax))).time,'linear','extrap');
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'b');
        end
        if length(IMU.(char(sensor_names(jj))).temp) == length(IMU.(char(sensor_names(jj))).time)      
            IMU.(char(sensor_names(jj))).temp = interp1(IMU.(char(sensor_names(jj))).time,IMU.(char(sensor_names(jj))).temp,IMU.(char(sensor_names(indmax))).time,'linear','extrap');
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'temp');
        end
          
        [match, ~] = ismember(IMU.(char(sensor_names(indmax))).time,IMU.(char(sensor_names(jj))).time);
        index_createALL = index(~match);
        
        new_q = zeros(length(IMU.(char(sensor_names(indmax))).time),4);
        new_q(match,:) = IMU.(char(sensor_names(jj))).q;

        for qq = 1:length(index_createALL)
            index_create = index_createALL(qq);     
            t1 = IMU.(char(sensor_names(indmax))).time(index_create - 1);
            t2 = IMU.(char(sensor_names(indmax))).time(index_create + 1);
            tc = IMU.(char(sensor_names(indmax))).time(index_create);

            t2 = t2 - t1;
            tc = tc - t1;
            tc = tc/t2;

            q1 = IMU.(char(sensor_names(jj))).q(index_create - 1,:);
            q2 = IMU.(char(sensor_names(jj))).q(index_create,:);
            qc = slerp(q1, q2, tc);

            new_q(index_create,:) = qc;

        end
        
        IMU.(char(sensor_names(jj))).q = new_q;      
    end
end

IMU.button.push_index = zeros(length(IMU.button.push_time),1);

for qq = 1:length(IMU.button.push_time)
    findmatch = IMU.(char(sensor_names(1))).time == IMU.button.push_time(qq);
    if sum(findmatch) == 1
        IMU.button.push_index(qq) = find(IMU.(char(sensor_names(1))).time == IMU.button.push_time(qq));
    else
        [~,IMU.button.push_index(qq)] = min(abs(double(IMU.(char(sensor_names(1))).time) - double(IMU.button.push_time(qq))));
    end
end

IMU.time = double((IMU.(char(sensor_names(1))).time - IMU.(char(sensor_names(1))).time(1)))/(10^6);
IMU.button.push_time = double((IMU.button.push_time - IMU.(char(sensor_names(1))).time(1)))/(10^6);

for ii = 1:num_sensors
    IMU.(char(sensor_names(ii))) = rmfield(IMU.(char(sensor_names(ii))),'time');
end

end