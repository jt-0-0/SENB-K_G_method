% ----------------------------------------------------------------------- %
% ------------------- Code to determine K_IC and G_IC ------------------- %
% ------------ In accordance to ISO 527-1 and ASTM D5045-14 ------------- %
% ----------------------------- Version 1.0 ----------------------------- %
% ------------------ Written by Sammy He and Joe Terry ------------------ %
% ---- Department of Mechanical Engineering, Imperial College London ---- %
% ------------------- sch12@ic.ac.uk, jst114@ic.ac.uk ------------------- %
% --- Copyright © 2022, Imperial College, London, All rights reserved --- %
% ---------------------- Last Updated: 09-03-2023 ----------------------- %
% ----------------------------------------------------------------------- %

% ------------------------------- SUMMARY ------------------------------- %

% This code will find the values for K_IC, and G_IC (via energy and LEFM 
% methods) for tests using multiple brittle plastic specimens, similar  
% to that in ISO 527-1 and ASTM D5045 - 14

% The code allows for a linear correction to the start of the load vs.
% extension data for compliance such that it passes through 0,0 and
% removes any initial load ramping errors.

% The user can select the region to calculate the linear gradient for
% specimen compliance and to then calculate P5.

% The code determines, Pmax, P5 and determines which to use. It then 
% calculates energy under the curve for sample and compliance, and 
% determines K_IC and G_IC (via energy and LEFM methods) for each specimen. 

% It will apply the validity criteria, specify valid/invalid specimens and
% it will present the mean average values and standard deviation of K_IC 
% and G_IC (via energy and LEFM methods) for valid specimens.

% It is difficult to control initial crack length when razor tapping 
% specimens, therefore the criteria specifying the length and variation of  
% this length, can be slightly loosened if desired by the user.

% ----------------------------------------------------------------------- %

%------------------------------ INSTRUCTIONS -----------------------------%

% It is reccomended that you make a new copy of this matlab script for each 
% type of specimen, so that it easy to re-run, should you so wish.

% Specimen data files should be comma seperated value files (.csv)
% Load should be in N and displacement in mm, otherwise there will be
% discprepancies with the units in subsequent calculations
% Specimen dimensions should be an excel spreadsheet (.xlsx)
% Dimension data should be inserted into a document following the 
% format of 'Specimen Dimensions K G Template'.
% Dimensions should be in mm.

% If you wish to add extra specimens to the dimensions excel sheet, you can
% copy and paste but do check the equation formulas are correct, especially 
% those in row F which were attached to fixed cells e.g. $C$5, $C$9 etc.

% To reduce probability of code errors, the compliance data should be 
% put into the same folder as the specimen test data. 
% File name for specimens should follow a system where the number at the
% end updates as the code iterates through the sample. e.g.
% 'Specimen_RawData_'. 
% This will become Specimen_RawData_1, Specimen_RawData_2 etc. 
% If preferred, the compliance (indentation) sample can follow this naming 
% scheme, as the last specimen.

% Often things can go wrong with testing (e.g. samples break before
% testing). For this code, the numbers of files must be consecutive.
% It is suggested to keep track of which samples correspond to the new 
% naming system in a readme file, if required for future reference.

% In this matlab script, in the USER INPUT section:
% - Write in a name for specimen type and add it's values for yield strength 
%   and Young's modulus from uniaxial tensile testing, these are required for
%   validity checks and G_IC_LEFM.
% - Insert file locations into the variables; directory_loc, indentdata_loc
%   and dimension data.
% - Input the start row of your data and the columns of load and extension 
%   from your specimen data files
% - Input the threshold of load drop compared to maximum at which you
%   would like to truncate test data, best to set this low e.g. 0.5 and increase if
%   required for your particular dataset.
% - Input the maxiumum load of compliance data to use during analysis
% - Make any adjustment to the validity criteria of the initial crack 
%   length (e.g. extend from a/w = 0.45-0.55 or increase allowed variation 
%   from average initial crack length from 10%.

% Then click run!

% A figure will appear with the indentation data. Select the linear region
% that is representative of the sample compliance.

% Afterwards, a figure will appear for each specimen, select two
% x-coordinates to define the range to determine the average linear 
% gradient, this will be used to find the Load 5 line.
% It is suggested to use this sensibly, select the linear section which is
% representative of the sample before the crack propagates. Some examples
% are shown in the help file.

% Once the code has run, it will produce a written summary of the result 
% for each specimen and mean values for valid specimens. It will also 
% produce comparison figures for all specimens, displaying their validity.

% If the user wishes to make the figure suitably sized for a report 
% (8 cm wide), they may open the figure and then paste the three below 
%  lines into the command window (uncomment first!)

% set(gcf,'Units', 'centimeters','Position',[10 10 8 6])
% set(gcf, [10 10 8 6]);
% set(gca, 'XMinorTick','on','YMinorTick','on','Layer', 'top','Units', 'centimeters','Position',[0.9 0.9 6.8 4.8])

% ----------------------------------------------------------------------- %

% Clear all in workspace, command window and close all other windows and
% figures.

clear
clc
clf
close all
warning('off') % this just removes a warning about missing items from legends

% ----------------------------------------------------------------------- %

%%-----------------------------USER INPUT--------------------------------%%

% Specimen type
spectype = 'Material_1';

% Tensile properties of the material from test data
E = 3; % Stiffness in GPa
sigma_y = 70; % Strength in MPa

% Add a link to the location with your raw inston results (usually ends in '_RawData')
% e.g. C:\Documents\Test_Results\Material_1\SENB\Material_1.is_flex_RawData
directory_loc = 'C:\Documents\Test_Results\Material_1\SENB\Material_1.is_flex_RawData';

% % Add a link to the file with your raw data for compliance test (ends in '_RawData_#.csv')
% e.g. C:\Documents\Test_Results\Material_1\SENB\Material_1.is_flex_RawData\Specimen_RawData_#.csv
indentdata_loc = 'C:\Documents\Test_Results\Material_1\SENB\Material_1.is_flex_RawData\Specimen_RawData_#.csv';

% Add a link to the location of the specimen dimension file (using provided template)
% e.g. C:\Documents\Test_Results\Material_1\SENB\Specimen_Dimensions
dimensiondata_loc = 'C:\Documents\Test_Results\Material_1\SENB\Specimen_Dimensions';

% Data start row (check in test data excel csv file which row your data starts for tests and indentation)
data_start_row = 6;

% Extension data column (check in test data excel csv file which column is strain for tests and indentation)
ext_col = 1;

% Load data column (check in test data excel csv file which column is load for tests and indentation)
load_col = 2;

% Set load drop threshold compared to max load (this is to delete extraneous data at the failure point, for a brittle material, 0.7-0.8 times maxload works well, but this can be altered if beneficial)
loaddropthresh = 0.8;

% Set max indent load (set to be 5-10 times larger than you max load from your testing, so there is chance to see the linear section)
max_indent_load = 500;

% Adjustment to validity criteria of samples

% Initial crack length min and max a/w (the standard suggests 0.45-0.55,
% however this can be difficult to control for brittle specimens which have
% been razor tapped hence can be potentially relaxed)
ICLmin = 0.45;
ICLmax = 0.55;

% Allowable variance of initial crack length (in %)(the standard uses 10%)
% this can be difficult to control for brittle specimens which have
% been razor tapped hence can be potentially relaxed
ICLalvar = 10;

%%-----------------------------END USER INPUT----------------------------%%

%%-----------------------------START OF CODE-----------------------------%%
% Define plot colours
col{1}  =   [0      0.447   0.741];     %'#0072BD'	standard blue
col{2}  =   [0.850  0.325   0.098];     %'#D95319'	orange
col{3}  =   [0.929  0.694   0.125];     %'#EDB120'	darker yellow
col{4}  =   [0.494  0.184   0.556];     %'#7E2F8E'	purple
col{5}  =   [0.466  0.674   0.188];     %'#77AC30'	green
col{6}  =   [0.301  0.745   0.933];     %'#4DBEEE'	lighter blue
col{7}  =   [0.635  0.078   0.184];     %'#A2142F'  burgandy
col{8}  =   [1      0       0];         %'#FF0000'  red
col{9}  =   [0      1       0];         %'#00FF00'	green
col{10} =   [0      0       1];         %'#0000FF'  blue
col{11} =   [0      1       1];         %'#00FFFF'	cyan
col{12} =   [1      0       1];         %'#FF00FF'	magenta
col{13} =   [1      1       0];         %'#FFFF00'	yellow
col{14} =   [0      0       0];         %'#000000'  black
col{15} =   [0      0.5     0];         %'#000000'  dark green
    
% Print specimen type 
disp(spectype);
disp(' ');

% ------------------- Indentation/Compliance Specimen ------------------- %

% Import indentation data
indentdata = xlsread(indentdata_loc); 

% Define indent load ext
indent_Ext = indentdata(data_start_row:end, ext_col);
indent_Load = indentdata(data_start_row:end, load_col);
indent_Ext(any(isnan(indent_Ext),2) , :) = []; % Remove all indentation data that are NaN 
indent_Load(any(isnan(indent_Load),2) , :) = []; % Remove all indentation data that are NaN 

% Delete indent data > max indent
[~,ind_del] = min(abs(max_indent_load-indent_Load));
indent_Ext((ind_del+1):end) = [];
indent_Load((ind_del+1):end) = [];

% Plot indentation data
figureindent = figure('Name','Compliance','Color','white','NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);% make figure full screen for accurate selection of gradient
axesindent = axes('Parent',figureindent);
plot(indent_Ext, indent_Load, 'color', col{4})
hold on
set(axesindent,'FontName','Calibri','FontSize',8)
set(gca, 'Layer', 'top')
box(axesindent,'on');
grid(axesindent,'on');
xlabel('Extension (mm)','FontWeight','bold','FontSize',10,...
    'FontName','Calibri');% x-axis label
ylabel('Load (N)','FontWeight','bold','FontSize',10,'FontName','Calibri') % y-axis label  

% Finding gradient of indentation

% Choose the x-axis range to search for the average gradient using ginput
[indentM_x1,~] = ginput(1);
[indentM_x2,~] = ginput(1);

% Find corresponding cell closest to the g-input value
[~,indentM_x1] = min(abs(indent_Ext-indentM_x1));
[~,indentM_x2] = min(abs(indent_Ext-indentM_x2));

% Switch gradient points if selected in incorrect order
if indentM_x1>indentM_x2
    indent_dummy = indentM_x1;
    indentM_x1 = indentM_x2;
    indentM_x2 = indent_dummy;
end

% Create new array over this defined interval
windowedxindent = indent_Ext(indentM_x1:indentM_x2);
windowedyindent = indent_Load(indentM_x1:indentM_x2);

% Fit a line to data within the interval 
[n, ~] = polyfit(windowedxindent, windowedyindent, 1);
grad_indent = n(1,1);
y_intercept_indent = n(1,2);

% Plot the linear gradient fit
fplot(@(x)grad_indent*x+y_intercept_indent, [0 max(indent_Ext)],'--','color',col{5});

% Plot the linear gradient fit with corrected displacement
compliancemaxlength=max(indent_Ext)-(- y_intercept_indent)/grad_indent;
indent_Ext_corr= 0:0.00001:compliancemaxlength;
indent_Load_corr=grad_indent*indent_Ext_corr;

% Plot corrected indentation data
plot(indent_Ext_corr,indent_Load_corr,'color',col{10});
set(gcf,'units','pixels','Position',[540 500 480 320]) % return figure to small size
legend('Location', 'northwest','NumColumnsMode','manual','NumColumns',1)
legend('Indent Data','Linear Fit', 'Corrected Indent Data')
savefig('Indentation');

% ----------------------------------------------------------------------- %

% -------------- Main loop to determine specimen toughness -------------- %

% 'dir' will count how many .csv files are in this folder
directory = dir([directory_loc, '\*.csv']);

% The number of specimens is one less than the size of directory (last is compliance test)
No_Spec=length(directory)-1;

% Pre-allocate for computing efficiency
senbdata = cell(1,No_Spec);    
W = zeros(1,No_Spec);
B = zeros(1,No_Spec);
a = zeros(1,No_Spec);
x = zeros(1,No_Spec);
f = zeros(1,No_Spec);
A = zeros(1,No_Spec);
dAdx = zeros(1,No_Spec);
phi = zeros(1,No_Spec);
Pmax = zeros(1,No_Spec);
Pmaxlength = zeros(1,No_Spec);
figureLE = zeros(1,No_Spec);
axesLE = zeros(1,No_Spec);
figurecomparison = zeros(1,No_Spec);
axescomparison = zeros(1,No_Spec);
M_x1 = zeros(1,No_Spec);
M_x2 = zeros(1,No_Spec);
M_dummy = zeros(1,No_Spec);
grad = zeros(1,No_Spec);
y_intercept_fit = zeros(1,No_Spec);
x_intercept_fit = zeros(1,No_Spec);
y_intercept_corr = zeros(1,No_Spec);
x_intercept_corr = zeros(1,No_Spec);
grad5 = zeros(1,No_Spec);
Load5 = zeros(1,No_Spec);
Pq = zeros(1,No_Spec);
P5length = zeros(1,No_Spec);
Uq = zeros(1,No_Spec);
Ui = zeros(1,No_Spec);
U = zeros(1,No_Spec);
Cq = zeros(1,No_Spec);
Ci = zeros(1,No_Spec);
Cc = zeros(1,No_Spec);
indentlength = zeros(1,No_Spec);
K_q = zeros(1,No_Spec);
G_q_energy = zeros(1,No_Spec);
G_q_LEFM = zeros(1,No_Spec);
v = zeros(1,No_Spec);
P = zeros(1,No_Spec);
K_IC = zeros(1,No_Spec);
G_IC_energy = zeros(1,No_Spec);
G_IC_LEFM = zeros(1,No_Spec);
E_stiff = zeros(1,No_Spec);
E_fract = zeros(1,No_Spec);
ACCratio = zeros(1,No_Spec);
r_bar = zeros(1,No_Spec);

% File name (do not write the number nor extension) 
% The default name is 'Specimen_RawData_#' if using raw Instron results
% the number will be added as the code iterates through the specimens.
filename = strcat(directory_loc,'\Specimen_RawData_'); 

% Import dimensions of  specimens from a spreadsheet in standard format. 
dimensiondata = xlsread(dimensiondata_loc);

for i = 1:No_Spec % Run for all specimens
% for i=5 % To run for one specimen, specify i as a single value, i.e. i=1 for specimen 1

    % Initially set sample as valid (this will change if fails
    % validity criteria later in flow)
    v(i)=1;
    
    % Read dimensions for the respective specimen
    W(i) = dimensiondata((4*i),3); % width (mm)
    B(i) = dimensiondata((4*i),4); % thickness (mm)
    a(i) = dimensiondata((4*i),6); % crack length (mm)
    x(i) = a(i)/W(i);

    % Calculate geometry based factors
    f(i) = (6*x(i)^0.5)*(1.99-x(i)*(1-x(i))*(2.15-3.93*x(i)+2.7*x(i)^2))/((1+2*x(i))*(1-x(i))^(3/2)); 
   
    A(i) = (16*x(i)^2/(1-x(i))^2)*...
    (8.9-33.717*x(i)+79.616*x(i)^2-112.952*x(i)^3+84.815*x(i)^4-25.672*x(i)^5);

    dAdx(i) = (16*x(i)^2/(1-x(i))^2)*...
    (-33.717+159.232*x(i)-338.856*x(i)^2+339.26*x(i)^3-128.36*x(i)^4)+16*...
    (8.9-33.717*x(i)+79.616*x(i)^2-112.952*x(i)^3+84.815*x(i)^4-25.672*x(i)^5)*...
    ((2*x(i)*(1-x(i))+2*x(i)^2)/(1-x(i))^3);

    phi(i) = (A(i)+18.64)/dAdx(i);
    
    % Import load extension data and define the load, extension and peak load.
    % each iteration, this will update \Specimen_RawData_1,Specimen_RawData_2 etc.
    num = num2str(i);
    file = append(filename,num,'.csv');  
    
    senbdata{i}= xlsread(file);      % Stores all raw data in a cell
    data = cell2mat(senbdata(i));
    disp(['Specimen ',num2str(i)])
    
    Ext = data(data_start_row:end, ext_col);
    Ext(any(isnan(Ext),2) , :) = []; % Remove all indentation data that are NaN 
    Load = data(data_start_row:end, load_col);
    Load(any(isnan(Load),2) , :) = []; % Remove all indentation data that are NaN 
    [Pmax(i),Pmaxlength(i)] = max(Load);

    % Delete extraneous data after specimen failure
    for j = Pmaxlength(i) : length(Load)
        if (Load(j)<loaddropthresh*Pmax(i))
            Load(j:end) = [];
            Ext(j:end) = [];
            break
        end
    end

    % Add a final array point at 0N, same displacement (to stop failure of
    % P5 line)
    d = size(Load)+1;
    Load(d) = 0;
    Ext(d) = Ext(d-1);
    
    % Plot load extension from test
    spectitle = ['Specimen ',num2str(i)];
    figureLE(i) = figure('Name',spectitle,'NumberTitle','off','Color','white','units','normalized','outerposition',[0 0 1 1]);
    axesLE(i) = axes('Parent',figureLE(i));
    plot(Ext, Load, 'color', col{8})
    hold on
    set(axesLE(i),'FontName','Calibri','FontSize',8)
    set(gca, 'Layer', 'top')
    box(axesLE(i),'on');
    grid(axesLE(i),'on');
    xlabel('Extension (mm)','FontWeight','bold','FontSize',10,...
        'FontName','Calibri');% x-axis label
    ylabel('Load (N)','FontWeight','bold','FontSize',10,'FontName','Calibri') % y-axis label  
    
    % Choose the x-axis range to search for the average gradient using ginput
    [M_x1(i),~] = ginput(1);
    [M_x2(i),~] = ginput(1);
    
    set(gcf,'units','pixels','Position',[540 500 480 320]) % return figure to small size

    % Find corresponding cell closest to the g-input value
    [~,M_x1(i)] = min(abs(Ext-M_x1(i)));
    [~,M_x2(i)] = min(abs(Ext-M_x2(i)));
    
    % Switch gradient points if selected in incorrect order (larger then
    % smaller)
    if  M_x1(i)>M_x2(i)
        M_dummy(i) = M_x1(i);
        M_x1(i) = M_x2(i);
        M_x2(i) = M_dummy(i);
    end

    % Create new array over this defined interval
    windowedx = Ext(M_x1(i):M_x2(i));
    windowedy = Load(M_x1(i):M_x2(i));

    % Fit a line to data within the interval 
    [n, ~] = polyfit(windowedx, windowedy, 1);
    grad(i) = n(1,1);
    y_intercept_fit(i) = n(1,2);
    x_intercept_fit(i) = -y_intercept_fit(i)/grad(i);

    % Plot the gradient fit on the graph
    fplot(@(x)grad(i)*x+y_intercept_fit(i), [0 max(Ext)],'--','color',col{9});
      
    % Find 0.95 gradient line and define Cq
    Cq(i) = 1/grad(i);
    grad5(i) = grad(i)*0.95;

    % 0.95 gradient line
    Ext_Load5=cat(1,Ext,(Ext(2:end)+ max(Ext))); % create array twice as long as extension
    Load5=grad5(i)*(Ext_Load5-x_intercept_fit(i)); % produce load5 line data
    
    % Find intersection between 5% and load ext plot, when is Load5 > Load?
    Loaddiff = Load5(1:max(size(Load))) - Load;
    is = find(Loaddiff>0);
    
    % Find position of largest continuous crossing by considering the
    % largest region that the Load5 line is less than the load ext plot.
    isdiff = diff(is);

    if  max(isdiff)>1 % If the line crosses the load ext plot more than once.
        [~,isloc] = max(isdiff);
        P5length(i) = is(isloc+1);
    else  % If the line only passes the load ext plot at Pq
        P5length(i) = is(1);
    end

    % Plot the 5% offset line
    plot(Ext(1:P5length(i)),Load5(1:P5length(i)),'-.','color',col{1});
    
    % Find the load at the point of crossing (Pq)
    Pq(i) = grad5(i)* (Ext(P5length(i))-x_intercept_fit(i)); 

    % Create variables for truncated load and extension
    Load_trunc=Load;
    Ext_trunc=Ext;
    
    if  P5length(i) < Pmaxlength(i) % If Pq happens before Pmax, use Pq            
        P(i) = Pq(i);

        % Truncate load ext data beyond P5length nd plot markers denoting
        % Pmax and P = Pq
        plot(Ext(Pmaxlength(i)),Load(Pmaxlength(i)),'o','MarkerSize',4,'color','black') % plot P max small
        Load_trunc((P5length(i)+1):end) = [];
        Ext_trunc((P5length(i)+1):end) = [];
        plot(Ext(P5length(i)),Load5(P5length(i)),'*','MarkerSize',10, 'color','black')
        legend('Location', 'northwest','NumColumnsMode','manual','NumColumns',2)
        legend('Test Data','S','0.95S','P_{max}','P = P_q')

    else % If Pmax happens before Pq, use Pmax            
        P(i) = Pmax(i);

        % Truncate load ext data beyond P5length and plot markers denoting
        % Pq and P = Pmax
        plot(Ext(P5length(i)),Load5(P5length(i)),'o','MarkerSize',4, 'color','black') % plot Pq small
        Load_trunc((Pmaxlength(i)+1):end) = [];
        Ext_trunc((Pmaxlength(i)+1):end) = [];  
        plot(Ext(Pmaxlength(i)),Load(Pmaxlength(i)),'*','MarkerSize',10, 'color','black'); 
        legend('Location', 'northwest','NumColumnsMode','manual','NumColumns',2)
        legend('Test Data','S','0.95S','P_q','P = P_{max}')
        
    end
              
    % Add axis limits for plot    
    xlim([0 (max(Ext)*1.05)]);
    ylim([0 (max(Load)*1.05)]);
   
    % Save Specimen # figure
    savefig(spectitle);
    
    % Truncate indent load and extenstion after P5length or Pmaxlength
    [~, indentlength(i)]  = min(abs(indent_Load_corr-P(i)));
    indent_Load_trunc=indent_Load_corr;
    indent_Ext_trunc=indent_Ext_corr;
    indent_Load_trunc((indentlength(i)+1):end) = [];
    indent_Ext_trunc((indentlength(i)+1):end) = [];
    
    % Plot indent + specimen load ext
    indentspecfig = ['Indent and specimen ',num2str(i)];
    figurecomparison(i)= figure('Name',indentspecfig,'color','white','NumberTitle','off');
    axescomparison(i)= axes('Parent',figurecomparison(i)); 
    plot(indent_Ext_trunc,indent_Load_trunc,'color',col{10});   
    hold on
    plot(Ext_trunc,Load_trunc,'color',col{8});
    set(axescomparison(i),'FontName','Calibri','FontSize',8)
    set(gcf,'Position',[540 500 480 320])
    set(gca, 'Layer', 'top')
    box(axescomparison(i),'on');
    grid(axescomparison(i),'on');
    xlabel('Extension (mm)','FontWeight','bold','FontSize',10,...
    'FontName','Calibri');% x-axis label
    ylabel('Load (N)','FontWeight','bold','FontSize',10,'FontName','Calibri') % y-axis label   
         
    % Plot markers for P, and Pmax or Pq
    if P5length(i) < Pmaxlength(i) % If Pq happens before Pmax, use Pq    
        plot(Ext(Pmaxlength(i)),Load(Pmaxlength(i)),'o','MarkerSize',4, 'Color','black') % plot P max small
        plot(Ext(P5length(i)),Load5(P5length(i)),'*','MarkerSize',10, 'Color','black')
        legend('Location', 'southeast','NumColumnsMode','manual','NumColumns',1)
        legend('Compliance','Test Data','P_{max}','P = P_q')
        lgd.NumColumns = 2;

    else % If Pmax happens before Pq, use Pmax 
        plot(Ext(P5length(i)),Load5(P5length(i)),'o','MarkerSize',4, 'Color','black') % plot Pq small
        plot(Ext(Pmaxlength(i)),Load(Pmaxlength(i)),'*','MarkerSize',10, 'Color','black');
        legend('Location', 'southeast','NumColumnsMode','manual','NumColumns',1)
        legend('Compliance','Test Data','P_q','P = P_{max}')
    end
            
    % Add axis limits for plot
    xlim([0 (max(Ext)*1.05)]);
    ylim([0 (Pmax(i)*1.05)]);     
   
    savefig(indentspecfig);
    
    % Calculate area under the respective curves for G_IC_Energy
    Uq(i) = trapz(Ext_trunc,Load_trunc)*10^-3; % units Joules
    Ui(i) = trapz(indent_Ext_trunc,indent_Load_trunc)*10^-3; % units Joules
    U(i) = Uq(i) - Ui(i); % units Joules
    
    % Calculate the corrected compliance
    Ci(i) = 1/grad_indent;
    Cc(i) = Cq(i) - Ci(i);

    % Final toughness calculation
    K_q(i) = ((P(i))/((B(i)*10^-3)*((W(i)*10^-3)^0.5)))*f(i)*10^-6; % units MPa m^1/2
    G_q_LEFM(i) = (1-0.35^2)*(K_q(i)*10^6)^2/(E*10^9); % units J m^-2
    G_q_energy(i) = (U(i)/((B(i)*10^-3)*(W(i)*10^-3)*phi(i))); % units J m^-2
    
    % Final toughness calculation
    disp(append('a =',' ',num2str(a(i)),' mm'))
    disp(append('K_IC =',' ',num2str(K_q(i)),' MPa m^{1/2}'))  
    disp(append('G_IC LEFM =',' ',num2str(G_q_LEFM(i)),' J m^{-2}'))  
    disp(append('G_IC Energy =',' ',num2str(G_q_energy(i)),' J m^{-2}'))  
    
    % Accuracy cross check calculations
    E_stiff(i) = ((2*f(i)^2)*phi(i))/((B(i)*10^-3)*(Cc(i)*10^3)); 
    E_fract(i) =  (K_q(i)*10^3)^2/G_q_LEFM(i); %K_IC^2/G_IC
    ACCratio(i) = E_stiff(i)/E_fract(i); % accuracy ratio
    
    % Calculate characteristic length r_bar
    r_bar(i) = ((K_q(i)*10^6)^2)/(sigma_y*10^6)^2;
    
    % Perform validity checks/give user warnings
    if ((x(i)<0.45||x(i)>0.55))
        disp(append('NOTE: x (a/W) is not between 0.45 and 0.55, x = ',num2str(x(i))))  
    end
    
    if dimensiondata((4*i),12)==1
    disp('NOTE: initial crack length varies by more than 20% across specimen width')
    
       elseif dimensiondata((4*i),11) == 1
       disp('NOTE: initial crack length varies by more than 10% across the specimen width') 
        
    end
    
    if E_fract(i)> E_stiff(i)
        disp('NOTE: E_fract is larger then E_stiff (usually vice-versa)')
    end
    
    if ACCratio(i) > 1.15
        disp('NOTE: E_fract is larger then E_stiff by more than 15%, examine results for possible errors')
    end
   
    if 2.5*( K_q(i)/sigma_y)^2 > B(i)/1000  % B in metres
    disp('Test is INVALID: size criteria have not been met, 2.5*( K_q/sigma_y)^2 > B/1000')
    v(i) = 0; 
    end
        
    if B(i)<2.5*r_bar(i)
        disp('Test is INVALID: thickness is less than 2.5 times the characteristic length')
        v(i) = 0; 
    end
        
    if a(i)<2.5*r_bar(i)
        disp('Test is INVALID: initial crack length is less than 2.5 times the characteristic length')
        v(i) = 0; 
    end
        
    if(W(i)-a(i))<2.5*r_bar(i)
        disp('Test is INVALID: ligament width is less than 2.5 times the characteristic length')
        v(i) = 0;  
    end
        
    if x(i)<ICLmin||x(i)>ICLmax % the values may be altered to 0.45 + 0.55 as per the standard
        % Display 'Test is INVALID: x (a/W) is not between user specified min and max values,',num2str(ICLmin), ' and ', num2str(ICLmax)), 'x = ',num2str(x(i)))         
        str1 = append('Test is INVALID: x (a/W) is not between user specified min and max values, ', num2str(ICLmin));
        str1 = append(str1, ' and ');
        str1 = append(str1, num2str(ICLmax));
        str1 = append(str1,', x = ');
        str1 = append(str1,num2str(x(i)));
        disp(str1)
        v(i)=0;
    end
        
    if dimensiondata((4*i-3),10) >= ICLalvar
        % Display 'Test is INVALID: variance of initial crack length, a, is above user specified max allowable value,', num2str(ICLalvar)),'%')
        str2 = append('Test is INVALID: variance of initial crack length, a, is above user specified max allowable value, ', num2str(ICLalvar));
        str2 = append(str2, ' %');
        disp(str2)
        v(i)=0;
    end
           
    if Pmaxlength(i) > P5length(i) && Pmax(i)/Pq(i) > 1.1
        disp('Test is INVALID: Pmax/Pq > 1.1')
        v(i) = 0;
    end
    
    if W/B > 4
        disp('Test is INVALID: width/thickness > 4')
        v(i) = 0;
    end
        
    if v(i)==1
        disp('Test is VALID')
        K_IC(i) = K_q(i);
        G_IC_LEFM(i) = G_q_LEFM(i);
        G_IC_energy(i) = G_q_energy(i);
        
    end
    
    disp(' ');
end

% ----------------------------------------------------------------------- %

% -------------- Display average values and overview plots -------------- %

% Calculate and display mean and standard deviation only for valid test
% results.
K_IC_avg = sum(K_IC)/sum(v);
G_IC_energy_avg = sum(G_IC_energy)/sum(v);
G_IC_LEFM_avg = sum(G_IC_LEFM)/sum(v);

% Remove invalidity values for std
K_IC_v = K_IC;
K_IC_v(K_IC_v==0) = [];
G_IC_energy_v = G_IC_energy;
G_IC_energy_v(G_IC_energy_v==0) = [];
G_IC_LEFM_v = G_IC_LEFM;
G_IC_LEFM_v(G_IC_LEFM_v==0) = [];

% Calculate standard deviation
K_IC_std = std(K_IC_v);
G_IC_energy_std = std(G_IC_energy_v);
G_IC_LEFM_std = std(G_IC_LEFM_v);

% Produce cell with results and export file
Spec_No = 1:No_Spec;

% SENBResults  = ["Spec. No.","Thickness (mm)","Width (mm)","Crack Length (mm)","K_IC (MPa m^{1/2})","G_IC Energy (J m^{-2})","G_IC LEFM (J m^{-2})";Spec_No.',B.',W.',a.',K_IC.',G_IC_Energy.', G_IC_LEFM.';"Average",0,0,0,K_IC_avg,G_IC_Energy_avg,G_IC_LEFM_avg;"Standard Deviation",0,0,0,K_IC_std,G_IC_Energy_std,G_IC_LEFM_std;]
SENBResults = table(Spec_No.',B.',W.',a.',K_IC.',G_IC_energy.', G_IC_LEFM.',v.');
header = {'Specimen','Thickness','Width','CrackLength','K_IC','G_ICEnergy','G_ICLEFM','Valid'};
SENBResults.Properties.VariableNames = header;

% Display the values for K_IC and G_IC via Energy and LEFM methods
disp(append('Summary of SENB results for ', spectype));
disp(append('Initial crack length, a/w between ',num2str(ICLmin),' and ',num2str(ICLmax)));
disp(append('Initial crack length, varience less than ',num2str(ICLalvar), '%'));
disp(append(num2str(No_Spec),' specimens tested, (',num2str(sum(v)),' valid and ', num2str(No_Spec-sum(v)),' invalid)' ));
disp(append('K_IC:        ',num2str(K_IC_avg),' +- ', num2str(K_IC_std),' MPa m^{1/2}'));
disp(append('G_IC Energy: ',num2str(G_IC_energy_avg),' +- ', num2str(G_IC_energy_std),' J m^{-2}'));
disp(append('G_IC LEFM:   ',num2str(G_IC_LEFM_avg),' +- ', num2str(G_IC_LEFM_std),' J m^{-2}'));
    
% Close all intermediate figures
close all

% Plot K_IC, G_IC LEFM and G_IC Energy for each specimen, including
% validity, mean and standard deviation

% K_IC plot
K_IC_overview = figure('Name','K_IC Overview','Color','white','NumberTitle','off');
K_IC_axes = axes('Parent',K_IC_overview);

% Plot mean and standard deviation
plot([0, max(Spec_No)*1.1],[K_IC_avg, K_IC_avg],'-','color', col{10});
hold on
plot([0, max(Spec_No)*1.1],[(K_IC_avg + K_IC_std), (K_IC_avg + K_IC_std)],'-.','color', col{10});
plot([0, max(Spec_No)*1.1],[(K_IC_avg - K_IC_std), (K_IC_avg - K_IC_std)],'-.','color', col{10},'HandleVisibility','off');
 
% Set up visuals (gridlines, font etc)
set(K_IC_axes,'FontName','Calibri','FontSize',8)
set(gcf,'Position',[540 500 480 320])
set(gca, 'Layer', 'top')
box(K_IC_axes,'on');
grid(K_IC_axes,'on');
xlabel('Specimen #','FontWeight','bold','FontSize',10,...
'FontName','Calibri');% x-axis label
ylabel('K_{IC} (MPa m^{1/2})','FontWeight','bold','FontSize',10,'FontName','Calibri') % y-axis label   

% Plot all specimens assuming invalid in red
plot(Spec_No, K_q,'*','color','red')

% Plot valid values on top in black
K_IC_plot = [K_IC;Spec_No];
K_IC_plot(:,any(K_IC_plot == 0))=[];
plot(K_IC_plot(2,:), K_IC_plot(1,:) ,'*','color',col{14})

% Set legend
legend('Location', 'south','NumColumnsMode','manual','NumColumns',2)
legend('Mean','SD','Invalid','Valid')

% Add axis limits for plot
xlim([0 (max(Spec_No)*1.05)]);
ylim([0 (max(K_q)*1.05)]);    

savefig('K_IC Overview.fig');

% G_IC_LEFM plot
G_IC_LEFM_overview = figure('Name','G_IC_LEFM Overview','Color','white','NumberTitle','off');
G_IC_LEFM_axes = axes('Parent',G_IC_LEFM_overview);

% Plot mean and standard deviation
plot([0, max(Spec_No)*1.1],[G_IC_LEFM_avg, G_IC_LEFM_avg],'-','color', col{10});
hold on
plot([0, max(Spec_No)*1.1],[(G_IC_LEFM_avg + G_IC_LEFM_std), (G_IC_LEFM_avg + G_IC_LEFM_std)],'-.','color', col{10});
plot([0, max(Spec_No)*1.1],[(G_IC_LEFM_avg - G_IC_LEFM_std), (G_IC_LEFM_avg - G_IC_LEFM_std)],'-.','color', col{10},'HandleVisibility','off');

% Set up visuals (gridlines, font etc)
set(G_IC_LEFM_axes,'FontName','Calibri','FontSize',8)
set(gcf,'Position',[540 500 480 320])
set(gca, 'Layer', 'top')
box(G_IC_LEFM_axes,'on');
grid(G_IC_LEFM_axes,'on');
xlabel('Specimen #','FontWeight','bold','FontSize',10,...
'FontName','Calibri');% x-axis label
ylabel('G_{IC, LEFM} (J m^{-2})','FontWeight','bold','FontSize',10,'FontName','Calibri') % y-axis label   

% Plot all specimens assuming invalid in red
plot(Spec_No, G_q_LEFM,'*','color','red')

% Plot valid values on top in black
G_IC_LEFM_plot = [G_IC_LEFM;Spec_No];
G_IC_LEFM_plot(:,any(G_IC_LEFM_plot == 0))=[];
plot(G_IC_LEFM_plot(2,:), G_IC_LEFM_plot(1,:) ,'*','color',col{14})

% Set Legend
legend('Location', 'south','NumColumnsMode','manual','NumColumns',2)
legend('Mean','SD','Invalid','Valid')

% Add axis limits for plot
xlim([0.8 (max(Spec_No)+0.2)]);

if max(G_q_energy)>max(G_q_LEFM)
    ylim([0.8 (max(G_q_energy)*1.05)]);  
else
    ylim([0.8 (max(G_q_LEFM)*1.05)]);
end

savefig('G_IC_LEFM Overview.fig');

% G_IC_energy plot
G_IC_energy_overview = figure('Name','G_IC_energy Overview','Color','white','NumberTitle','off');
G_IC_energy_axes = axes('Parent',G_IC_energy_overview);

% Plot mean and standard deviation
plot([0, max(Spec_No)*1.1],[G_IC_energy_avg, G_IC_energy_avg],'-','color', col{10});
hold on
plot([0, max(Spec_No)*1.1],[(G_IC_energy_avg + G_IC_energy_std), (G_IC_energy_avg + G_IC_energy_std)],'-.','color', col{10});
plot([0, max(Spec_No)*1.1],[(G_IC_energy_avg - G_IC_energy_std), (G_IC_energy_avg - G_IC_energy_std)],'-.','color', col{10},'HandleVisibility','off');

% Set up visuals (gridlines, font etc)
set(G_IC_energy_axes,'FontName','Calibri','FontSize',8)
set(gcf,'Position',[540 500 480 320])
set(gca, 'Layer', 'top')
box(G_IC_energy_axes,'on');
grid(G_IC_energy_axes,'on');
xlabel('Specimen #','FontWeight','bold','FontSize',10,...
'FontName','Calibri');% x-axis label
ylabel('G_{IC, energy} (J m^{-2})','FontWeight','bold','FontSize',10,'FontName','Calibri') % y-axis label   

% Plot all specimens assuming invalid in red
plot(Spec_No, G_q_energy,'*','color','red')

% Plot valid values on top in black
G_IC_energy_plot = [G_IC_energy;Spec_No];
G_IC_energy_plot(:,any(G_IC_energy_plot == 0))=[];
plot(G_IC_energy_plot(2,:), G_IC_energy_plot(1,:) ,'*','color',col{14})

% Set legend
legend('Location', 'south','NumColumnsMode','manual','NumColumns',2)
legend('Mean','SD','Invalid','Valid')

% Add axis limits for plot
xlim([0.8 (max(Spec_No)+0.2)]);

if max(G_q_energy)>max(G_q_LEFM)
    ylim([0.8 (max(G_q_energy)*1.05)]);  
else
    ylim([0.8 (max(G_q_LEFM)*1.05)]);
end

savefig('G_IC_energy Overview.fig');

% ----------------------------------------------------------------------- %
    
% Close figures and save workspace
close all
save(spectype);

%%----------------------------- END OF CODE -----------------------------%%

% ----------------------------------------------------------------------- %
% ------------------- Code to determine K_IC and G_IC ------------------- %
% ------------ In accordance to ISO 527-1 and ASTM D5045-14 ------------- %
% ----------------------------- Version 1.0 ----------------------------- %
% ------------------ Written by Sammy He and Joe Terry ------------------ %
% ---- Department of Mechanical Engineering, Imperial College London ---- %
% ------------------- sch12@ic.ac.uk, jst114@ic.ac.uk ------------------- %
% --- Copyright © 2022, Imperial College, London, All rights reserved --- %
% ---------------------- Last Updated: 09-03-2023 ----------------------- %
% ----------------------------------------------------------------------- %
