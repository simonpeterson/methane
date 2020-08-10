%{
Script to create plots for data analysis of methane flux
created by Simon Peterson for Nicholas Hasson
8/8/20
%}

%read in the results of the winter data
methane_data = readmatrix("results_winter_data_kwa-lab_07-28-20.xlsx",'NumHeaderLines',1,'OutputType','string','sheet',2);
%convert the id to integers and the data to doubles
flux_values = str2double(methane_data(:,3));
ID = str2double(methane_data(:,4));
%the measurement locations
locations = ["lake ice", 0; "littoral typha",1; "littoral sedge",2; "active margin",3; "upland forest", 4];
%create cell array to hold flux values for each location. 1x1 matrix is
%placeholding matrix
locations_data = {zeros(1),0;zeros(1),1;zeros(1),2;zeros(1),3;zeros(1),4};
%create the data for each location
single_flux_data = zeros(1);
for i = 0:4
    single_flux_data = double.empty;
    for p = 1:size(ID,1)
        if ID(p) == i
            single_flux_data = cat(2,single_flux_data,flux_values(p));
        end
    end
    %add to specified cell array
    locations_data{i+1,1} = single_flux_data;
end
%calculate the average values of the data
averages = zeros(5,2);
%set average to ID
averages(:,1) = locations(:,2);
for i = 0:4
    avg = 0;
    count = 0;
    for p = 1:size(ID,1)
        if ID(p) == i
            avg = avg + flux_values(p);
            count = count + 1;
        end
    end
    %shifted because matlab starts counting at one
    averages(i+1,2) = avg/count;
end
%% make plot of the averages
scatter(averages(:,1),averages(:,2),'filled');
set(gca,'xtick',[0:4],'xticklabel',locations(:,1));
title("Flux vs. Measurement Location: Average");
ylabel('CH_{4} umol m ^{-2} s ^{-1}','fontsize',20);
xlabel('Surface Type','fontsize',20);

%% Plot each of the datasets, match to a gaussian
%plot the data for 
histogram(locations_data{1,1},20);
