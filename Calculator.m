clear; clc; close all;
%% 0. Initialize Parameters
L = 1200; % Length of bridge
n = 1200; % Discretize into 1 mm seg.
dL = L/n; % Length of each segment
P = 400; % Total weight of train [N]
x = linspace(0, L, n+1); % x-axis

%% 1. SFD, BMD under train loading
%% 1.1 Define Train Loading
x_train = [52 228 392 568 732 908]; % Train Load Locations (the 6 wheels)
l_train = 960; % Train Length
P_factors = [1 1 1 1 1 1]; % Load Factors (the 6 wheels)
P_train = P_factors .* P/sum(P_factors) * -1; % load of each wheel

%% 1.2 Solve for SFD and BMD with the train at different locations
% set up the SFD and BMD arrays of the bridge for every discretized location of the train
% the rows are the locations of the train, with the first row being the train completely off of the bridge on the left and the last row being the train completely off of the bridge on the right
% the columns are the forces/moments at each cut, with the first column being the left end of the bridge and the last column being the right end of the bridge
SFDi = zeros(n + 1 + l_train / dL, n + 1); % SFDs of the bridge
BMDi = zeros(n + 1 + l_train / dL, n + 1); % BMDs of the bridge

% shift the train entirely off of the bridge to the left
x_train = x_train - l_train;

% for every discretized location of the train, find the SFD and BMD
for i = 1:n+1+l_train
    R1 = 0; % Reaction 1
    R2 = 0; % Reaction 2
    
    % find the loads that are on the bridge
    x_train_on_bridge = x_train(x_train >= 0 & x_train <= L);
    P_train_on_bridge = P_train(x_train >= 0 & x_train <= L);

    % if there are loads on the bridge, solve for the reactions
    if ~isempty(x_train_on_bridge)
        % solve for the reactions
        [R1, R2] = reaction(P_train_on_bridge, x_train_on_bridge);
    end

    % solve for the shear force and bending moment at each location
    for j = 1:n+1
        % j is the cut number
        % find the location of the cut
        cut = x(j);

        % find the loads that are within the cut
        x_train_in_cut = [x_train_on_bridge(x_train_on_bridge <= cut), 0, 0];
        P_train_in_cut = [P_train_on_bridge(x_train_on_bridge <= cut), R1, 0];

        % if the cut is at the end of the bridge, add the reaction at the other end
        if cut == x(end)
            P_train_in_cut(end) = R2;
            x_train_in_cut(end) = L;
        end

        % if there are loads within the cut, solve for the internal forces
        if ~isempty(x_train_in_cut)
            % solve for the internal forces
            [shear, moment] = internal_forces(P_train_in_cut, x_train_in_cut, cut);
        else
            % if there are no loads within the cut, the internal forces are 0
            shear = 0;
            moment = 0;
        end

        % add the internal forces to the SFD and BMD
        SFDi(i,j) = shear;
        BMDi(i,j) = moment;

    end

    x_train = x_train + dL; % move the train to the right by 1 increment
end

%% 1.3 Plot SFD and BMD
SFD = max(abs(SFDi)); % SFD envelope
BMD = max(BMDi); % BMD envelope

% plot the SFD and BMD envelopes on separate figures
% SFD is red
% BMD is blue
figure
hold on; grid on; grid minor;
plot(x, SFD, 'r', 'LineWidth', 2)
%plot(x, -SFD, 'r')
plot(x, zeros(1, n+1), 'k', 'LineWidth', 4)

title("Shear Force Envelope")
xlabel("Location on Bridge (mm)")
ylabel("Shear Force (N)")
legend('Shear Force')

figure
hold on; grid on; grid minor;
% invert the y axis so that positive bending moments are plotted downwards
set(gca, 'YDir','reverse')
plot(x, BMD, 'b', 'LineWidth', 2)
%plot(x, -BMD, 'b')
plot(x, zeros(1, n+1), 'k', 'LineWidth', 4)

title("Bending Moment Envelope")
xlabel("Location on Bridge (mm)")
ylabel("Bending Moment (N*mm)")
legend('Bending Moment')

%% 1.4 Display the maximum shear force and bending moment
max_shear = max(SFD);
disp("Maximum shear force: " + sigfig(max_shear) + " N")

max_moment = max(BMD);
disp("Maximum bending moment: " + sigfig(max_moment) + " N*mm")

disp(" ")

% find when the moment is at a maximum
max_moment_location = x(BMD == max_moment);

% find which row of the BMDi matrix the maximum moment is in column max_moment_location(1)

%% 2. Define Cross Section
%% 2.1 Define Individual Cross Sections
% x_sections is stored as a dictionary with the keys being the x location along the bridge and the values being the cross section at that location
% each cross section is stored as an array of subsections with the following values:
% [x, y, dx, dy, lc, id; ...]
% x is the x location of the reference point of the subsection
% y is the y location of the reference point of the subsection
% dx is the width of the subsection
% dy is the height of the subsection
% lc is the load case of the subsection (from 1 to 3)
% lc 0 is used for subsections that do not need to be considered
% id is for subsections that are stacked on top of each other
% subsections stacked on top of each other have the same id
% the reference point is the bottom left corner of the subsection
% the sections are constructed from the bottom up with 0, 0 being the bottom left corner of the bridge
% x_change = [0]; % x locations of cross section changes
% x_sections = {[10+1.27, 0, 80-2*1.27, 1.27, 0, 0; % bottom flange
%                10, 0, 1.27, 75+1.27, 3, 0; % left web
%                90-1.27, 0, 1.27, 75+1.27, 0, 0; % right web 
%                10+1.27, 75-1.27, 5, 1.27, 0, 0; % left glue connection
%                90-5-1.27, 75-1.27, 5, 1.27, 0, 0; % right glue connection
%                0, 75, 10, 1.27, 2, 0; % left top flange
%                10 + 1.27, 75, 100 - 20 - 2*1.27, 1.27, 1, 0; % center top flange
%                90, 75, 10, 1.27, 2, 1]}; % right top flange

% location of the diaphragms along the bridge
% diaphragms = [0, 400, 800, 1200];

%% 2.2 Define Glue
% The locations of the glue are stored as a dictionary with the keys being the x location along the bridge and the values being the locations of the glue
% each glue location is stored as an array of subsections with the following values:
% [dir, x, y, length/width, calculate; ...]
% dir is the direction of the glue with 0 being horizontal and 1 being vertical (vertical not implemented)
% x is the x location of the reference point of the glue
% y is the y location of the reference point of the glue
% length/width is the length or width of the glue
% calculate is either 0 or 1 with 1 being the glue that needs to be considered for shear
% the reference point is the bottom left corner of the glue
% % the glues are constructed from the bottom up with 0, 0 being the bottom left corner of the bridge 
% glue_locations = {[0, 10, 75, 1.27 + 5, 1;
%                    0, 10 + 80 - 5 - 1.27, 75, 1.27 + 5, 1]};

% x_section_params = dictionary(x_change, x_sections); % dictionary of the cross sections
% glue_params = dictionary(x_change, glue_locations); % dictionary of the glue locations

%% 2.3 Design 0
x_change = [0]; % x locations of cross section changes
x_sections = {[10+1.27, 0, 80-2*1.27, 1.27, 0, 0; % bottom flange
               10, 0, 1.27, 75+1.27, 3, 0; % left web
               90-1.27, 0, 1.27, 75+1.27, 0, 0; % right web 
               10+1.27, 75-1.27, 5, 1.27, 0, 0; % left glue connection
               90-5-1.27, 75-1.27, 5, 1.27, 0, 0; % right glue connection
               0, 75, 10, 1.27, 2, 0; % left top flange
               10 + 1.27, 75, 100 - 20 - 2*1.27, 1.27, 1, 0; % center top flange
               90, 75, 10, 1.27, 2, 1]}; % right top flange
diaphragms = [0, 400, 800, 1200];
glue_locations = {[0, 10, 75, 1.27 + 5, 1;
                   0, 10 + 80 - 5 - 1.27, 75, 1.27 + 5, 1]};
% %% 2.4 Final Design
% x_change = [0,600];
% x_sections = {[
%     %webs from left to right
%     21.865, 0, 1.27, 120, 3, 0;
%     96.865, 0, 1.27, 120, 3, 0;
%     %left upper flanges
%     20, 116.19, 1.865, 1.27, 2, 1;
%     0, 117.46, 21.865, 1.27, 2, 1;
%     0, 118.73, 21.865, 1.27, 2, 1;
%     %center flanges
%     23.135, 116.19, 73.73, 1.27, 1, 2;
%     23.135, 117.46, 73.73, 1.27, 1, 2;
%     23.135, 118.73, 73.73, 1.27, 1, 2;
%     %upper right flanges
%     98.135, 116.19, 1.865, 1.27, 2, 3;
%     98.135, 117.46, 21.865, 1.27, 2, 3;
%     98.135, 118.73, 21.865, 1.27, 2, 3;
%     %upper glue tabs
%     23.135, 114.92, 10, 1.27, 0, 0;
%     86.865, 114.92, 10, 1.27, 0, 0;
%     ],[
%     %webs from left to right
%     21.865, 0, 1.27, 120, 3, 0;
%     96.865, 0, 1.27, 120, 3, 0;
%     %left upper flanges
%     20, 116.19, 1.865, 1.27, 2, 1;
%     0, 117.46, 21.865, 1.27, 2, 1;
%     0, 118.73, 21.865, 1.27, 2, 1;
%     %center flanges
%     23.135, 116.19, 73.73, 1.27, 1, 2;
%     23.135, 117.46, 73.73, 1.27, 1, 2;
%     23.135, 118.73, 73.73, 1.27, 1, 2;
%     %upper right flanges
%     98.135, 116.19, 1.865, 1.27, 2, 3;
%     98.135, 117.46, 21.865, 1.27, 2, 3;
%     98.135, 118.73, 21.865, 1.27, 2, 3;
%     %upper glue tabs
%     23.135, 114.92, 10, 1.27, 0, 0;
%     86.865, 114.92, 10, 1.27, 0, 0;
%     ]};

% glue_locations = {[
%     0, 21.865, 116.19, 11.27, 1;
%     0, 86.855, 116.19, 11.27, 1;
%     0, 20, 117.46, 80, 0;
%     0, 20, 118.73, 80, 0;
%     ],[
%     0, 21.865, 116.19, 11.27, 1;
%     0, 86.855, 116.19, 11.27, 1;
%     0, 20, 117.46, 80, 0;
%     0, 20, 118.73, 80, 0;
%     ]};
% diaphragms = [0, 50, 96, 188, 326, 510, 740, 924, 1062, 1154, 1200, 1250];
%% 2.5 Mapping the Values
x_section_params = dictionary(x_change, x_sections); % dictionary of the cross sections
glue_params = dictionary(x_change, glue_locations); % dictionary of the glue locations
%% 3. Calculate Sectional Properties
%% 3.1 Setup arrays to store the sectional properties
A =     zeros(length(x_change), 1); % areas
ybar =  zeros(length(x_change), 1); % location of centroidal axis from the bottom of the cross section
ybot =  zeros(length(x_change), 1); % location of bottom of cross section from the bottom of the cross section
ytop =  zeros(length(x_change), 1); % location of top of cross section from the bottom of the cross section
I =     zeros(length(x_change), 1); % second moment of area
Qcent = zeros(length(x_change), 1); % Q at centroidal axes
Qglue = zeros(length(x_change), 1); % Q at glue location
gL =    zeros(length(x_change), 1); % length of glue
bg =    zeros(length(x_change), 1); % width of glue
b =     zeros(length(x_change), 1); % width of cross section that crosses the centroidal axis

%% 3.2 Calculate the sectional properties for each cross section
for i = 1:length(x_change)
    % find the cross section
    x_section = x_section_params(x_change(i));
    x_section = x_section{1, 1};

    % setup a table for the cross section
    x_section = array2table(x_section, 'VariableNames', {'x', 'y', 'dx', 'dy', 'lc', 'id'});

    % add a column for the area of each subsection
    x_section.area = x_section.dx.*x_section.dy;

    % add a column for the centroidal axis of each subsection
    x_section.ybar = x_section.y + x_section.dy/2;

    % add a column for the second moment of area of each subsection
    x_section.I = x_section.dx.*x_section.dy.^3/12;

    % add a column for the bottom of each subsection
    x_section.ybot = x_section.y;

    % add a column for the top of each subsection
    x_section.ytop = x_section.y + x_section.dy;

    % find the total area
    A(i) = sum(x_section.area);

    % find the centroidal axis
    ybar(i) = sum(x_section.area.*x_section.ybar)/A(i);

    % find the bottom of the cross section
    ybot(i) = min(x_section.y);

    % find the top of the cross section
    ytop(i) = max(x_section.y + x_section.dy);

    % find the second moment of area using the parallel axis theorem
    % I = I_cent + A*d^2
    I(i) = sum(x_section.I + x_section.area.*(x_section.ybar - ybar(i)).^2);

    % find Q at the centroidal axis
    % get the subsections with ybar > ybot
    x_section_cent_bot = x_section(ybar(i) > x_section.ybot, :);

    % if the top of the subsection is above the centroidal axis, cut the subsection to include only parts below the centroidal axis
    x_section_cent_bot.dy(x_section_cent_bot.ytop > ybar(i)) = ybar(i) - x_section_cent_bot.ybot(x_section_cent_bot.ytop > ybar(i));

    % for the subsections that were cut, add up their widths to get b
    b(i) = sum(x_section_cent_bot.dx(x_section_cent_bot.ytop > ybar(i)));

    % recalculate the area of the subsections
    x_section_cent_bot.area = x_section_cent_bot.dx.*x_section_cent_bot.dy;

    % recalculate the centroidal axis of the subsections
    x_section_cent_bot.ybar = x_section_cent_bot.y + x_section_cent_bot.dy/2;

    % calculate the centroidal axis of the subsections combined
    ybar_cent_bot = sum(x_section_cent_bot.area.*x_section_cent_bot.ybar)/sum(x_section_cent_bot.area);

    % find Q at the centroidal axis
    Qcent(i) = sum(x_section_cent_bot.area.*(ybar(i) - ybar_cent_bot));

    % find Q at the glue location
    % find the glue locations of the cross section that need to be calculated
    glue = glue_params(x_change(i));
    glue = glue{1, 1};
    % sum the widths of the glue
    gL(i) = sum(glue(:,4));

    % get the glue locations that need to be calculated
    glue = glue(glue(:,5) == 1, :);
    
    % if there are glue locations, calculate Q at the glue locations

    % check if the glue is horizontal or vertical
    % assuming that all glue for a given cross section is either all horizontal or all vertical
    if ~isempty(glue)
        bg = 0; % width of the glue
        if glue(1, 1) == 0
            % only calculate the glue if calculate is 1
                % if the glue is horizontal, find the width of the glue
                % find the subsections with y less than the glue y
                x_section_glue = x_section(x_section.y < glue(1,3), :);

                % trim the subsections that are above the glue by setting it so that their y + dy is equal to the glue y
                x_section_glue.dy(x_section_glue.y + x_section_glue.dy > glue(1,3)) = glue(1,3) - x_section_glue.y(x_section_glue.y + x_section_glue.dy > glue(1,3));

                % recalculate the area of the subsections
                x_section_glue.area = x_section_glue.dx.*x_section_glue.dy;

                % recalculate the centroidal axis of the subsections
                x_section_glue.ybar = x_section_glue.y + x_section_glue.dy/2;

                % find the width of the glue
                bg = sum(glue(:,4));

                % find Q at the glue location
                Qglue(i) = sum(x_section_glue.area.*(ybar(i) - x_section_glue.ybar));

                bg(i) = bg;
        else
            % not implemented
        end
        
    end

end

%% 3.3 Fill the table with the sectional properties
bridge_properties = array2table(zeros(n+1, 8), "VariableNames", ["x", "A", "ybar", "ybot", "ytop", "I", "Qcent", "Qglue"]);

bridge_properties.x = x.'; % x locations

for i = 1:length(x_change)
    % find the x locations that are within the cross section
    x_in_section = bridge_properties.x >= x_change(i);

    % add the properties to the table
    bridge_properties.A(x_in_section) = A(i);
    bridge_properties.ybar(x_in_section) = ybar(i);
    bridge_properties.ybot(x_in_section) = ybot(i);
    bridge_properties.ytop(x_in_section) = ytop(i);
    bridge_properties.I(x_in_section) = I(i);
    bridge_properties.Qcent(x_in_section) = Qcent(i);
    bridge_properties.Qglue(x_in_section) = Qglue(i);
    bridge_properties.b(x_in_section) = b(i);
    bridge_properties.gL(x_in_section) = gL(i);
    bridge_properties.bg(x_in_section) = bg(i);

    disp("Cross Section at x = " + x_change(i) + " mm")
    disp("Area: " + sigfig(A(i)) + " mm^2")
    disp("Centroidal Axis: " + sigfig(ybar(i)) + " mm")
    disp("Bottom of Cross Section: " + sigfig(ybot(i)) + " mm")
    disp("Top of Cross Section: " + sigfig(ytop(i)) + " mm")
    disp("Second Moment of Area: " + sigfig(I(i)) + " mm^4")
    disp("Q at Centroidal Axis: " + sigfig(Qcent(i)) + " mm^3")
    disp("Q at Glue Location: " + sigfig(abs(Qglue(i))) + " mm^3")
    disp("Width of Cross Section that Crosses Centroidal Axis: " + sigfig(b(i)) + " mm")
    disp("Length of Glue: " + sigfig(gL(i)) + " mm")
    disp(" ")
end

%% 3.4 Calculate the amount of material needed
matboard = sum(bridge_properties.A);
disp("Matboard needed: " + sigfig(matboard) + " mm^3")

%% 3.5 Calculate the amount of glue needed
glue_needed = sum(bridge_properties.gL);
disp("Glue needed: " + sigfig(glue_needed) + " mm^2")
disp(" ")

%% 3.6 Plotting the Cross Sections
for i = 1:length(x_change)
    % find the cross section
    x_section = x_section_params(x_change(i));
    x_section = x_section{1, 1};

    % find the glue locations
    glue = glue_params(x_change(i));
    glue = glue{1, 1};

    % plot the cross section
    figure
    hold on; grid on; grid minor;
    for j = 1:size(x_section,1)
        % determine the color of the subsection based on the load case
        if x_section(j,5) == 0
            color = 'k';
        elseif x_section(j,5) == 1
            color = 'b';
        elseif x_section(j,5) == 2
            color = 'r';
        elseif x_section(j,5) == 3
            color = 'g';
        end
        rectangle('Position', [x_section(j,1), x_section(j,2), x_section(j,3), x_section(j,4)], 'FaceColor', color)
    end

    % plot the glue
    for j = 1:size(glue,1)
        if glue(j,1) == 0
            rectangle('Position', [glue(j,2), glue(j,3) - 0.05, glue(j,4), 0.1], 'FaceColor', 'y')
        else
            rectangle('Position', [glue(j,2) - 0.05, glue(j,3), 0.1, glue(j,4)], 'FaceColor', 'y')
        end
    end

    % plot the centroidal axis
    plot([0, max(x_section(:,1) + x_section(:,3))], [bridge_properties.ybar(bridge_properties.x == x_change(i)), bridge_properties.ybar(bridge_properties.x == x_change(i))], 'k--', 'LineWidth', 2, "DisplayName", "Centroidal Axis")

    % plot the x axis up to the width of the cross section
    % find the width of the cross section
    width = max(x_section(:,1) + x_section(:,3));

    % plot the x axis
    plot([0, width], [0, 0], 'k', 'LineWidth', 2)

    % plot the y axis
    % find the height of the cross section
    height = max(x_section(:,2) + x_section(:,4));

    % plot the y axis
    plot([0, 0], [0, height], 'k', 'LineWidth', 2)

    % set the axis limits
    xlim([0, width])
    ylim([0, height])

    % label the axes
    xlabel('x (mm)')
    ylabel('y (mm)')

    % title the figure
    title(['Cross Section at x = ', sigfig(x_change(i)), ' mm'])
end

%% 4. Calculate Applied Stress
%% 4.1 Calculating Applied Stresses
% stress at the top across the entire bridge
bridge_properties.sigma_top = BMD.'.*abs(bridge_properties.ytop - bridge_properties.ybar)./bridge_properties.I;

% stress at the bottom across the entire bridge
bridge_properties.sigma_bot = BMD.'.*abs(bridge_properties.ybot - bridge_properties.ybar)./bridge_properties.I;

% shear stress across the entire bridge at the centroidal axis
bridge_properties.tau_xy = SFD.'.*bridge_properties.Qcent./bridge_properties.I./bridge_properties.b;

% shear stress across the entire bridge at the glue locations
bridge_properties.tau_g = SFD.'.*bridge_properties.Qglue./bridge_properties.I./bridge_properties.bg;

%% 4.2 Displaying the Maximum Applied Stresses
% maximum stress at the top
max_sigma_top = max(bridge_properties.sigma_top);
disp("Maximum stress at the top: " + sigfig(max_sigma_top) + " MPa")

% maximum stress at the bottom
max_sigma_bot = max(bridge_properties.sigma_bot);
disp("Maximum stress at the bottom: " + sigfig(max_sigma_bot) + " MPa")

% maximum shear stress at the centroidal axis
max_tau_xy = max(abs(bridge_properties.tau_xy));
disp("Maximum shear stress at the centroidal axis: " + sigfig(max_tau_xy) + " MPa")

% maximum shear stress at the glue locations
max_tau_g = max(abs(bridge_properties.tau_g));
disp("Maximum shear stress at the glue locations: " + sigfig(max_tau_g) + " MPa")

disp(" ")

%% 5. Material and Thin Plate Buckling Capacities
%% 5.1 Setup a table to store the capacities of the material
material_properties = array2table(zeros(n+1, 11), "VariableNames", ["x", "E", "mu", "sigma_tens", "sigma_comp", "tau_max", "tau_gmax", "sigma_buck_flange_webs", "sigma_buck_flange_tips", "sigma_buck_webs", "tau_buck_webs"]);

material_properties.E = 4000 * ones(n+1, 1);
material_properties.mu = 0.2 * ones(n+1, 1);

material_properties.sigma_tens = 30 * ones(n+1, 1);
material_properties.sigma_comp = 6 * ones(n+1, 1);
material_properties.tau_max = 4 * ones(n+1, 1);

material_properties.tau_gmax = 2 * ones(n+1, 1);

%% 5.2 Look through each cross section and calculate the buckling capacities
for i = 1:length(x_change)
    x_section = x_section_params(x_change(i));

    % get all subsections that are assigned load case 1
    x_sections_lc1 = x_section{1, 1}(x_section{1, 1}(:,5) == 1, :);

    % combine subsections that are stacked on top of each other
    % use the id column to combine subsections that are stacked on top of each other
    x_sections_lc1_combined = [];
    
    % get the unique ids
    ids = unique(x_sections_lc1(:,6));
    for j = 1:length(ids)
        % get the subsections with the current id
        x_sections_lc1_id = x_sections_lc1(x_sections_lc1(:,6) == ids(j), :);
        % combine the subsections
        x_sections_lc1_combined = [x_sections_lc1_combined; sum(x_sections_lc1_id(:,3)), x_sections_lc1_id(1,2), x_sections_lc1_id(1,3), sum(x_sections_lc1_id(:,4)), x_sections_lc1_id(1,5)];
    end

    x_sections_lc1 = x_sections_lc1_combined;

    sigma_buck1 = inf;
    for j = 1:size(x_sections_lc1, 1)
        % calculate buckling capacity for each subsection
        % if it is the lowest buckling capacity, set buck1 to that value and highlight the subsection in black
        % if it is not the lowest buckling capacity, highlight the subsection in blue
        % sigma = (4 * pi^2 * E) / (12 * (1 - mu^2)) * (t / b)^2
        % E and mu are constants
        % t is the height of the horizontal flange
        % b is the width of the horizontal flange
        sigma_buck1_i = (4 * pi^2 * material_properties.E(1)) / (12 * (1 - material_properties.mu(1)^2)) * (x_sections_lc1(j,4) / x_sections_lc1(j,3))^2;
        if sigma_buck1_i < sigma_buck1
            sigma_buck1 = sigma_buck1_i;
        end

    end

    % get all subsections that are assigned load case 2
    x_sections_lc2 = x_section{1, 1}(x_section{1, 1}(:,5) == 2, :);

    % combine subsections that are stacked on top of each other
    % use the id column to combine subsections that are stacked on top of each other
    x_sections_lc2_combined = [];

    % get the unique ids
    ids = unique(x_sections_lc2(:,6));
    for j = 1:length(ids)
        % get the subsections with the current id
        x_sections_lc2_id = x_sections_lc2(x_sections_lc2(:,6) == ids(j), :);
        % combine the subsections
        x_sections_lc2_combined = [x_sections_lc2_combined; sum(x_sections_lc2_id(:,3)), x_sections_lc2_id(1,2), x_sections_lc2_id(1,3), sum(x_sections_lc2_id(:,4)), x_sections_lc2_id(1,5)];
    end

    x_sections_lc2 = x_sections_lc2_combined;

    sigma_buck2 = inf;
    for j = 1:size(x_sections_lc2, 1)
        % calculate buckling capacity for each subsection
        % if it is the lowest buckling capacity, set buck2 to that value and highlight the subsection in black
        % if it is not the lowest buckling capacity, highlight the subsection in blue
        % sigma = (0.425 * pi^2 * E) / (12 * (1 - mu^2)) * (t / b)^2
        % E and mu are constants
        % t is the height of the horizontal flange
        % b is the width of the horizontal flange
        sigma_buck2_i = (0.425 * pi^2 * material_properties.E(1)) / (12 * (1 - material_properties.mu(1)^2)) * (x_sections_lc2(j,4) / x_sections_lc2(j,3))^2;
        if sigma_buck2_i < sigma_buck2
            sigma_buck2 = sigma_buck2_i;
        end
    end

    % get all subsections that are assigned load case 3
    x_sections_lc3 = x_section{1, 1}(x_section{1, 1}(:,5) == 3, :);

    sigma_buck3 = inf;
    for j = 1:size(x_sections_lc3, 1)
        % calculate buckling capacity for each subsection
        % if it is the lowest buckling capacity, set buck3 to that value and highlight the subsection in black
        % if it is not the lowest buckling capacity, highlight the subsection in blue
        % sigma = (6 * pi^2 * E) / (12 * (1 - mu^2)) * (t / b)^2
        % b is the distance between the centroidal axis and the top of the subsection
        % E and mu are constants
        % find the ybar of the bridge at this location
        % look for the ybar of the bridge at x = x_change(i)
        ybar_i = bridge_properties.ybar(bridge_properties.x == x_change(i));

        sigma_buck3_i = (6 * pi^2 * material_properties.E(1)) / (12 * (1 - material_properties.mu(1)^2)) * (x_sections_lc3(j,3) / (x_sections_lc3(j,2) + x_sections_lc3(j, 4) - ybar_i))^2;
        if sigma_buck3_i < sigma_buck3
            sigma_buck3 = sigma_buck3_i;
        end
    end

    % fill the table with the buckling capacities until the next cross section without overfilling
    if i < length(x_change)
        material_properties.sigma_buck_flange_webs(x_change(i) <= bridge_properties.x & bridge_properties.x < x_change(i+1)) = sigma_buck1;
        material_properties.sigma_buck_flange_tips(x_change(i) <= bridge_properties.x & bridge_properties.x < x_change(i+1)) = sigma_buck2;
        material_properties.sigma_buck_webs(x_change(i) <= bridge_properties.x & bridge_properties.x < x_change(i+1)) = sigma_buck3;
    else
        material_properties.sigma_buck_flange_webs(x_change(i) <= bridge_properties.x) = sigma_buck1;
        material_properties.sigma_buck_flange_tips(x_change(i) <= bridge_properties.x) = sigma_buck2;
        material_properties.sigma_buck_webs(x_change(i) <= bridge_properties.x) = sigma_buck3;
    end
end

%% 5.3 Calculate the shear buckling capacity of the webs
% tau = (5 * pi^2 * E) / (12 * (1 - mu^2)) * ( (t/h) ^ 2 + (t/a) ^ 2)
% E and mu are constants
% t is the thickness of the web
% h is the height of the web
% a is the spacing between diaphragms
% calculate the shear buckling capacity of the webs between each diaphragm
for i = 1:length(diaphragms)-1
    start_x = diaphragms(i);
    end_x = diaphragms(i+1);
    a = end_x - start_x;
    h = bridge_properties.ytop(bridge_properties.x == start_x) - bridge_properties.ybot(bridge_properties.x == start_x);
    t = 1.27;
    tau_buck_webs = (5 * pi^2 * material_properties.E(1)) / (12 * (1 - material_properties.mu(1)^2)) * ( (t/h) ^ 2 + (t/a) ^ 2);
    material_properties.tau_buck_webs(bridge_properties.x >= start_x & bridge_properties.x <= end_x) = tau_buck_webs;
end

%% 5.4 Display the Minimum Critical Stresses
crit_stresses = array2table(zeros(1, 8), 'VariableNames', ["min_sigma_tens", "min_sigma_comp", "min_tau_max", "min_tau_gmax", "min_sigma_buck_flange_webs", "min_sigma_buck_flange_tips", "min_sigma_buck_webs", "min_tau_buck_webs"]);
crit_stresses.min_sigma_tens = min(material_properties.sigma_tens);
crit_stresses.min_sigma_comp = min(material_properties.sigma_comp);
crit_stresses.min_tau_max = min(material_properties.tau_max);
crit_stresses.min_tau_gmax = min(material_properties.tau_gmax);
crit_stresses.min_sigma_buck_flange_webs = min(material_properties.sigma_buck_flange_webs);
crit_stresses.min_sigma_buck_flange_tips = min(material_properties.sigma_buck_flange_tips);
crit_stresses.min_sigma_buck_webs = min(material_properties.sigma_buck_webs);
crit_stresses.min_tau_buck_webs = min(material_properties.tau_buck_webs);
disp("All Stress Values are in MPa")
disp("Summary of Critical Stresses (not rounded)")
disp(crit_stresses)
% display the minimum critical stresses with the correct number of significant figures
disp("Critical Tensile Stress: " + sigfig(crit_stresses.min_sigma_tens))
disp("Critical Compressive Stress: " + sigfig(crit_stresses.min_sigma_comp))
disp("Critical Shear Stress: " + sigfig(crit_stresses.min_tau_max))
disp("Critical Glue Shear Stress: " + sigfig(crit_stresses.min_tau_gmax))
disp("Critical Buckling Stress for Flange Between Webs: " + sigfig(crit_stresses.min_sigma_buck_flange_webs))
disp("Critical Buckling Stress for Flange Tips: " + sigfig(crit_stresses.min_sigma_buck_flange_tips))
disp("Critical Buckling Stress for Webs: " + sigfig(crit_stresses.min_sigma_buck_webs))
disp("Critical Shear Buckling Stress for Webs: " + sigfig(crit_stresses.min_tau_buck_webs))
disp(" ")

%% 6. FOS
%% 6.1 Calculate FOS
bridge_properties.FOS_tens = abs(material_properties.sigma_tens./bridge_properties.sigma_bot);
bridge_properties.FOS_comp = abs(material_properties.sigma_comp./bridge_properties.sigma_top);
bridge_properties.FOS_shear = abs(material_properties.tau_max./bridge_properties.tau_xy);
bridge_properties.FOS_glue = abs(material_properties.tau_gmax./bridge_properties.tau_g);
bridge_properties.FOS_flange_webs = abs(material_properties.sigma_buck_flange_webs./bridge_properties.sigma_top);
bridge_properties.FOS_flange_tips = abs(material_properties.sigma_buck_flange_tips./bridge_properties.sigma_top);
bridge_properties.FOS_webs = abs(material_properties.sigma_buck_webs./bridge_properties.sigma_top);
bridge_properties.FOS_webs_shear = abs(material_properties.tau_buck_webs./bridge_properties.tau_xy);
% replace all NaNs with inf
bridge_properties.FOS_tens(isnan(bridge_properties.FOS_tens)) = inf;
bridge_properties.FOS_comp(isnan(bridge_properties.FOS_comp)) = inf;
bridge_properties.FOS_shear(isnan(bridge_properties.FOS_shear)) = inf;
bridge_properties.FOS_glue(isnan(bridge_properties.FOS_glue)) = inf;
bridge_properties.FOS_flange_webs(isnan(bridge_properties.FOS_flange_webs)) = inf;
bridge_properties.FOS_flange_tips(isnan(bridge_properties.FOS_flange_tips)) = inf;
bridge_properties.FOS_webs(isnan(bridge_properties.FOS_webs)) = inf;
bridge_properties.FOS_webs_shear(isnan(bridge_properties.FOS_webs_shear)) = inf;

%% 6.2 Find the minimum FOS
min_FOSs = array2table(zeros(1, 9), 'VariableNames', ["min_FOS_tens", "min_FOS_comp", "min_FOS_shear", "min_FOS_glue", "min_FOS_flange_webs", "min_FOS_flange_tips", "min_FOS_webs", "min_FOS_webs_shear", "overall_min_FOS"]);
min_FOSs.min_FOS_tens = min(bridge_properties.FOS_tens);
min_FOSs.min_FOS_comp = min(bridge_properties.FOS_comp);
min_FOSs.min_FOS_shear = min(bridge_properties.FOS_shear);
min_FOSs.min_FOS_glue = min(bridge_properties.FOS_glue);
min_FOSs.min_FOS_flange_webs = min(bridge_properties.FOS_flange_webs);
min_FOSs.min_FOS_flange_tips = min(bridge_properties.FOS_flange_tips);
min_FOSs.min_FOS_webs = min(bridge_properties.FOS_webs);
min_FOSs.min_FOS_webs_shear = min(bridge_properties.FOS_webs_shear);
min_FOSs.overall_min_FOS = min([min_FOSs.min_FOS_tens, min_FOSs.min_FOS_comp, min_FOSs.min_FOS_shear, min_FOSs.min_FOS_glue, min_FOSs.min_FOS_flange_webs, min_FOSs.min_FOS_flange_tips, min_FOSs.min_FOS_webs, min_FOSs.min_FOS_webs_shear]);
disp("Summary of Minimum FOSs (not rounded)")
disp(min_FOSs)
% display the minimum FOSs with the correct number of significant figures
disp("Minimum FOS Against Tensile Failure of Walls: " + sigfig(min_FOSs.min_FOS_tens))
disp("Minimum FOS Against Compressive Failure of Walls: " + sigfig(min_FOSs.min_FOS_comp))
disp("Minimum FOS Against Shear Failure of Walls: " + sigfig(min_FOSs.min_FOS_shear))
disp("Minimum FOS Against Glue Shear Failure: " + sigfig(min_FOSs.min_FOS_glue))
disp("Minimum FOS Against Flange Between Webs Buckling Failure: " + sigfig(min_FOSs.min_FOS_flange_webs))
disp("Minimum FOS Against Flange Tips Buckling Failure: " + sigfig(min_FOSs.min_FOS_flange_tips))
disp("Minimum FOS Against Webs Buckling Failure: " + sigfig(min_FOSs.min_FOS_webs))
disp("Minimum FOS Against Webs Shear Buckling Failure: " + sigfig(min_FOSs.min_FOS_webs_shear))
disp("Overall Minimum FOS: " + sigfig(min_FOSs.overall_min_FOS))
disp(" ")

%% 6.3 Plotting the FOSs on separate graphs
figure
subplot(2, 3, 1)
hold on; grid on; grid minor;
plot(x, bridge_properties.FOS_shear, 'r:', 'LineWidth', 2)
title("Shear Failure of Walls")
xlabel("Location on Bridge (mm)")
ylabel("Factor of Safety")
plot(x(bridge_properties.FOS_shear == min_FOSs.min_FOS_shear), min_FOSs.min_FOS_shear, 'rx')
% label the first instance of the minimum FOS
for i = 1:length(x)
    if bridge_properties.FOS_shear(i) == min_FOSs.min_FOS_shear
        text(x(i), min_FOSs.min_FOS_shear, sigfig(min_FOSs.min_FOS_shear))
        break
    end
end


legend("FOS", "Minimum FOS")
% change the limits to not include excessively high FOSs
ylim([0, max([min_FOSs.min_FOS_shear])*5])

subplot(2, 3, 2)
hold on; grid on; grid minor;
plot(x, bridge_properties.FOS_glue, 'g:', 'LineWidth', 2)
title("Glue Shear Failure")
xlabel("Location on Bridge (mm)")
ylabel("Factor of Safety")
plot(x(bridge_properties.FOS_glue == min_FOSs.min_FOS_glue), min_FOSs.min_FOS_glue, 'rx')
for i = 1:length(x)
    if bridge_properties.FOS_glue(i) == min_FOSs.min_FOS_glue
        text(x(i), min_FOSs.min_FOS_glue, sigfig(min_FOSs.min_FOS_glue))
        break
    end
end
legend("FOS", "Minimum FOS")
ylim([0, max([min_FOSs.min_FOS_glue])*5])

subplot(2, 3, 3)
hold on; grid on; grid minor;
plot(x, bridge_properties.FOS_webs_shear, 'b:', 'LineWidth', 2)
title("Webs Buckling Failure")
xlabel("Location on Bridge (mm)")
ylabel("Factor of Safety")
plot(x(bridge_properties.FOS_webs_shear == min_FOSs.min_FOS_webs_shear), min_FOSs.min_FOS_webs_shear, 'rx')
for i = 1:length(x)
    if bridge_properties.FOS_webs_shear(i) == min_FOSs.min_FOS_webs_shear
        text(x(i), min_FOSs.min_FOS_webs_shear, sigfig(min_FOSs.min_FOS_webs_shear))
        break
    end
end
legend("FOS", "Minimum FOS")
ylim([0, max([min_FOSs.min_FOS_webs_shear])*5])

subplot(2, 3, 4)
hold on; grid on; grid minor;
plot(x, bridge_properties.FOS_tens, 'r', 'LineWidth', 2)
plot(x, bridge_properties.FOS_comp, 'g', 'LineWidth', 2)
title("Flexural Failure of Walls")
xlabel("Location on Bridge (mm)")
ylabel("Factor of Safety")
plot(x(bridge_properties.FOS_tens == min_FOSs.min_FOS_tens), min_FOSs.min_FOS_tens, 'rx')
for i = 1:length(x)
    if bridge_properties.FOS_tens(i) == min_FOSs.min_FOS_tens
        text(x(i), min_FOSs.min_FOS_tens, sigfig(min_FOSs.min_FOS_tens))
        break
    end
end
plot(x(bridge_properties.FOS_comp == min_FOSs.min_FOS_comp), min_FOSs.min_FOS_comp, 'rx')
for i = 1:length(x)
    if bridge_properties.FOS_comp(i) == min_FOSs.min_FOS_comp
        text(x(i), min_FOSs.min_FOS_comp, sigfig(min_FOSs.min_FOS_comp))
        break
    end
end
legend("FOS Against Flexture Tensile Failure of Walls", "FOS Against Flexural Compressive Failure of Walls", "Minimum FOS", "")
ylim([0, max([min_FOSs.min_FOS_tens, min_FOSs.min_FOS_comp])*5])

subplot(2, 3, 5)
hold on; grid on; grid minor;
plot(x, bridge_properties.FOS_flange_webs, 'r-.', 'LineWidth', 2)
plot(x, bridge_properties.FOS_flange_tips, 'g-.', 'LineWidth', 2)
title("Flange Buckling Failure")
xlabel("Location on Bridge (mm)")
ylabel("Factor of Safety")
plot(x(bridge_properties.FOS_flange_webs == min_FOSs.min_FOS_flange_webs), min_FOSs.min_FOS_flange_webs, 'rx')
for i = 1:length(x)
    if bridge_properties.FOS_flange_webs(i) == min_FOSs.min_FOS_flange_webs
        text(x(i), min_FOSs.min_FOS_flange_webs, sigfig(min_FOSs.min_FOS_flange_webs))
        break
    end
end
plot(x(bridge_properties.FOS_flange_tips == min_FOSs.min_FOS_flange_tips), min_FOSs.min_FOS_flange_tips, 'rx')
for i = 1:length(x)
    if bridge_properties.FOS_flange_tips(i) == min_FOSs.min_FOS_flange_tips
        text(x(i), min_FOSs.min_FOS_flange_tips, sigfig(min_FOSs.min_FOS_flange_tips))
        break
    end
end
legend("FOS Against Flange Between Webs Buckling Failure", "FOS Against Flange Tips Buckling Failure", "Minimum FOS", "")
ylim([0, max([min_FOSs.min_FOS_flange_webs, min_FOSs.min_FOS_flange_tips])*5])

subplot(2, 3, 6)
hold on; grid on; grid minor;
plot(x, bridge_properties.FOS_webs, 'b-.', 'LineWidth', 2)
title("Webs Buckling Failure")
xlabel("Location on Bridge (mm)")
ylabel("Factor of Safety")
plot(x(bridge_properties.FOS_webs == min_FOSs.min_FOS_webs), min_FOSs.min_FOS_webs, 'rx')
for i = 1:length(x)
    if bridge_properties.FOS_webs(i) == min_FOSs.min_FOS_webs
        text(x(i), min_FOSs.min_FOS_webs, sigfig(min_FOSs.min_FOS_webs))
        break
    end
end
legend("FOS", "Minimum FOS")
ylim([0, max([min_FOSs.min_FOS_webs])*5])

sgtitle("Factors of Safety")

%% 6.4 Plotting the FOSs on the same graph
fosfig = figure;
hold on; grid on; grid minor;
fos1 = plot(x, bridge_properties.FOS_tens.', 'r-', 'LineWidth', 2);
% display the minimum FOS as an x
plot(x(bridge_properties.FOS_tens == min_FOSs.min_FOS_tens), min_FOSs.min_FOS_tens, 'rx')
% label the minimum FOS
for i = 1:length(x)
    if bridge_properties.FOS_tens(i) == min_FOSs.min_FOS_tens
        text(x(i), min_FOSs.min_FOS_tens, sigfig(min_FOSs.min_FOS_tens))
        break
    end
end

fos2 = plot(x, bridge_properties.FOS_comp.', 'g-', 'LineWidth', 2);
plot(x(bridge_properties.FOS_comp == min_FOSs.min_FOS_comp), min_FOSs.min_FOS_comp, 'rx')
for i = 1:length(x)
    if bridge_properties.FOS_comp(i) == min_FOSs.min_FOS_comp
        text(x(i), min_FOSs.min_FOS_comp, sigfig(min_FOSs.min_FOS_comp))
        break
    end
end

fos3 = plot(x, abs(bridge_properties.FOS_shear.'), 'r:', 'LineWidth', 2);
plot(x(bridge_properties.FOS_shear == min_FOSs.min_FOS_shear), min_FOSs.min_FOS_shear, 'rx')
for i = 1:length(x)
    if bridge_properties.FOS_shear(i) == min_FOSs.min_FOS_shear
        text(x(i), min_FOSs.min_FOS_shear, sigfig(min_FOSs.min_FOS_shear))
        break
    end
end

fos4 = plot(x, abs(bridge_properties.FOS_glue.'), 'g:', 'LineWidth', 2);
plot(x(bridge_properties.FOS_glue == min_FOSs.min_FOS_glue), min_FOSs.min_FOS_glue, 'rx')
for i = 1:length(x)
    if bridge_properties.FOS_glue(i) == min_FOSs.min_FOS_glue
        text(x(i), min_FOSs.min_FOS_glue, sigfig(min_FOSs.min_FOS_glue))
        break
    end
end

fos5 = plot(x, bridge_properties.FOS_flange_webs.', 'r-.', 'LineWidth', 2);
plot(x(bridge_properties.FOS_flange_webs == min_FOSs.min_FOS_flange_webs), min_FOSs.min_FOS_flange_webs, 'rx')
for i = 1:length(x)
    if bridge_properties.FOS_flange_webs(i) == min_FOSs.min_FOS_flange_webs
        text(x(i), min_FOSs.min_FOS_flange_webs, sigfig(min_FOSs.min_FOS_flange_webs))
        break
    end
end

fos6 = plot(x, bridge_properties.FOS_flange_tips.', 'g-.', 'LineWidth', 2);
plot(x(bridge_properties.FOS_flange_tips == min_FOSs.min_FOS_flange_tips), min_FOSs.min_FOS_flange_tips, 'rx')
for i = 1:length(x)
    if bridge_properties.FOS_flange_tips(i) == min_FOSs.min_FOS_flange_tips
        text(x(i), min_FOSs.min_FOS_flange_tips, sigfig(min_FOSs.min_FOS_flange_tips))
        break
    end
end

fos7 = plot(x, bridge_properties.FOS_webs.', 'b-.', 'LineWidth', 2);
plot(x(bridge_properties.FOS_webs == min_FOSs.min_FOS_webs), min_FOSs.min_FOS_webs, 'rx')
for i = 1:length(x)
    if bridge_properties.FOS_webs(i) == min_FOSs.min_FOS_webs
        text(x(i), min_FOSs.min_FOS_webs, sigfig(min_FOSs.min_FOS_webs))
        break
    end
end

fos8 = plot(x, abs(bridge_properties.FOS_webs_shear.'), 'b:', 'LineWidth', 2);
plot(x(bridge_properties.FOS_webs_shear == min_FOSs.min_FOS_webs_shear), min_FOSs.min_FOS_webs_shear, 'rx')
for i = 1:length(x)
    if bridge_properties.FOS_webs_shear(i) == min_FOSs.min_FOS_webs_shear
        text(x(i), min_FOSs.min_FOS_webs_shear, sigfig(min_FOSs.min_FOS_webs_shear))
        break
    end
end

% plot the minimum FOS line
fos9 = plot(x, ones(n+1, 1)*min_FOSs.overall_min_FOS, 'r--', 'LineWidth', 4);
text(x(1), min_FOSs.overall_min_FOS, sigfig(min_FOSs.overall_min_FOS))

% plot the FOS = 1 line
fos10 = plot(x, ones(n+1, 1), 'k--', 'LineWidth', 2);

title("Factors of Safety")
xlabel("Location on Bridge (mm)")
ylabel("Factor of Safety")
legend([fos1, fos2, fos3, fos4, fos5, fos6, fos7, fos8, fos9, fos10], ["Flexture Tensile Failure of Walls", "Flexural Compressive Failure of Walls", "Shear Failure of Walls", "Glue Shear Failure", "Flange Between Webs Buckling Failure", "Flange Tips Buckling Failure", "Webs Buckling Failure", "Webs Shear Buckling Failure", "Minimum FOS", "FOS = 1"])

ylim([0, min_FOSs.overall_min_FOS*10])

set(findall(fosfig,'-property','FontSize'),'FontSize',16)

%% 7. PFail
%% 7.1 Calculating using FOS
% Failure Mode 1: Flexural Tensile Failure of Walls
bridge_properties.P_flex_tens = abs(bridge_properties.FOS_tens.*P);
% Failure Mode 2: Flexural Compressive Failure of Walls
bridge_properties.P_flex_comp = abs(bridge_properties.FOS_comp.*P);
% Failure Mode 3: Shear Failure of Walls
bridge_properties.P_shear = abs(bridge_properties.FOS_shear.*P);
% Failure Mode 4: Glue Shear Failure
bridge_properties.P_glue = abs(bridge_properties.FOS_glue.*P);
% Failure Mode 5: Flange Between Webs Buckling Failure
bridge_properties.P_flange_webs = abs(bridge_properties.FOS_flange_webs.*P);
% Failure Mode 6: Flange Tips Buckling Failure
bridge_properties.P_flange_tips = abs(bridge_properties.FOS_flange_tips.*P);
% Failure Mode 7: Webs Buckling Failure
bridge_properties.P_webs = abs(bridge_properties.FOS_webs.*P);
% Failure Mode 8: Webs Shear Buckling Failure
bridge_properties.P_webs_shear = abs(bridge_properties.FOS_webs_shear.*P);

%% 7.2 Finding minimum failure load for each failure mode
minP = array2table(zeros(1, 9), 'VariableNames', ["minP_flex_tens", "minP_flex_comp", "minP_shear", "minP_glue", "minP_flange_webs", "minP_flange_tips", "minP_webs", "minP_webs_shear", "overall_minP"]);
minP.minP_flex_tens = min(bridge_properties.P_flex_tens);
minP.minP_flex_comp = min(bridge_properties.P_flex_comp);
minP.minP_shear = min(bridge_properties.P_shear);
minP.minP_glue = min(bridge_properties.P_glue);
minP.minP_flange_webs = min(bridge_properties.P_flange_webs);
minP.minP_flange_tips = min(bridge_properties.P_flange_tips);
minP.minP_webs = min(bridge_properties.P_webs);
minP.minP_webs_shear = min(bridge_properties.P_webs_shear);
minP.overall_minP = min([minP.minP_flex_tens, minP.minP_flex_comp, minP.minP_shear, minP.minP_glue, minP.minP_flange_webs, minP.minP_flange_tips, minP.minP_webs, minP.minP_webs_shear]);
disp("All values are in Newtons")
disp("Summary of Minimum Failure Loads (not rounded)")
disp(minP)
% display the minimum failure loads with the correct number of significant figures
disp("Minimum Failure Load Against Tensile Failure of Walls: " + sigfig(minP.minP_flex_tens))
disp("Minimum Failure Load Against Compressive Failure of Walls: " + sigfig(minP.minP_flex_comp))
disp("Minimum Failure Load Against Shear Failure of Walls: " + sigfig(minP.minP_shear))
disp("Minimum Failure Load Against Glue Shear Failure: " + sigfig(minP.minP_glue))
disp("Minimum Failure Load Against Flange Between Webs Buckling Failure: " + sigfig(minP.minP_flange_webs))
disp("Minimum Failure Load Against Flange Tips Buckling Failure: " + sigfig(minP.minP_flange_tips))
disp("Minimum Failure Load Against Webs Buckling Failure: " + sigfig(minP.minP_webs))
disp("Minimum Failure Load Against Webs Shear Buckling Failure: " + sigfig(minP.minP_webs_shear))
disp("Overall Minimum Failure Load: " + sigfig(minP.overall_minP))
disp(" ")

%% 7.3 Plotting the failure loads on separate graphs
figure
subplot(2, 3, 1)
hold on; grid on; grid minor;
plot(x, bridge_properties.P_shear, 'r:', 'LineWidth', 2)
title("Shear Failure of Walls")
xlabel("Location on Bridge (mm)")
ylabel("Failure Load (N)")
plot(x(bridge_properties.P_shear == minP.minP_shear), minP.minP_shear, 'rx')
for i = 1:length(x)
    if bridge_properties.P_shear(i) == minP.minP_shear
        text(x(i), minP.minP_shear, sigfig(minP.minP_shear))
        break
    end
end
legend("Failure Load", "Minimum Failure Load")
ylim([0, max([minP.minP_shear])*5])

subplot(2, 3, 2)
hold on; grid on; grid minor;
plot(x, bridge_properties.P_glue, 'g:', 'LineWidth', 2)
title("Glue Shear Failure")
xlabel("Location on Bridge (mm)")
ylabel("Failure Load (N)")
plot(x(bridge_properties.P_glue == minP.minP_glue), minP.minP_glue, 'rx')
for i = 1:length(x)
    if bridge_properties.P_glue(i) == minP.minP_glue
        text(x(i), minP.minP_glue, sigfig(minP.minP_glue))
        break
    end
end
legend("Failure Load", "Minimum Failure Load")
ylim([0, max([minP.minP_glue])*5])

subplot(2, 3, 3)
hold on; grid on; grid minor;
plot(x, bridge_properties.P_webs_shear, 'b:', 'LineWidth', 2)
title("Webs Buckling Failure")
xlabel("Location on Bridge (mm)")
ylabel("Failure Load (N)")
plot(x(bridge_properties.P_webs_shear == minP.minP_webs_shear), minP.minP_webs_shear, 'rx')
for i = 1:length(x)
    if bridge_properties.P_webs_shear(i) == minP.minP_webs_shear
        text(x(i), minP.minP_webs_shear, sigfig(minP.minP_webs_shear))
        break
    end
end
legend("Failure Load", "Minimum Failure Load")
ylim([0, max([minP.minP_webs_shear])*5])

subplot(2, 3, 4)
hold on; grid on; grid minor;
plot(x, bridge_properties.P_flex_tens, 'r', 'LineWidth', 2)
plot(x, bridge_properties.P_flex_comp, 'g', 'LineWidth', 2)
title("Flexural Failure of Walls")
xlabel("Location on Bridge (mm)")
ylabel("Failure Load (N)")
plot(x(bridge_properties.P_flex_tens == minP.minP_flex_tens), minP.minP_flex_tens, 'rx')
for i = 1:length(x)
    if bridge_properties.P_flex_tens(i) == minP.minP_flex_tens
        text(x(i), minP.minP_flex_tens, sigfig(minP.minP_flex_tens))
        break
    end
end
plot(x(bridge_properties.P_flex_comp == minP.minP_flex_comp), minP.minP_flex_comp, 'rx')
for i = 1:length(x)
    if bridge_properties.P_flex_comp(i) == minP.minP_flex_comp
        text(x(i), minP.minP_flex_comp, sigfig(minP.minP_flex_comp))
        break
    end
end
legend("Flexture Tensile Failure Load of Walls", "Flexural Compressive Failure Load of Walls", "Minimum Failure Load", "")
ylim([0, max([minP.minP_flex_tens, minP.minP_flex_comp])*5])

subplot(2, 3, 5)
hold on; grid on; grid minor;
plot(x, bridge_properties.P_flange_webs, 'r-.', 'LineWidth', 2)
plot(x, bridge_properties.P_flange_tips, 'g-.', 'LineWidth', 2)
title("Flange Buckling Failure")
xlabel("Location on Bridge (mm)")
ylabel("Failure Load (N)")
plot(x(bridge_properties.P_flange_webs == minP.minP_flange_webs), minP.minP_flange_webs, 'rx')
for i = 1:length(x)
    if bridge_properties.P_flange_webs(i) == minP.minP_flange_webs
        text(x(i), minP.minP_flange_webs, sigfig(minP.minP_flange_webs))
        break
    end
end
plot(x(bridge_properties.P_flange_tips == minP.minP_flange_tips), minP.minP_flange_tips, 'rx')
for i = 1:length(x)
    if bridge_properties.P_flange_tips(i) == minP.minP_flange_tips
        text(x(i), minP.minP_flange_tips, sigfig(minP.minP_flange_tips))
        break
    end
end
legend("Flange Between Webs Buckling Failure", "Flange Tips Buckling Failure", "Minimum Failure Load", "")
ylim([0, max([minP.minP_flange_webs, minP.minP_flange_tips])*5])

subplot(2, 3, 6)
hold on; grid on; grid minor;
plot(x, bridge_properties.P_webs, 'b-.', 'LineWidth', 2)
title("Webs Buckling Failure")
xlabel("Location on Bridge (mm)")
ylabel("Failure Load (N)")
plot(x(bridge_properties.P_webs == minP.minP_webs), minP.minP_webs, 'rx')
for i = 1:length(x)
    if bridge_properties.P_webs(i) == minP.minP_webs
        text(x(i), minP.minP_webs, sigfig(minP.minP_webs))
        break
    end
end
legend("Failure Load", "Minimum Failure Load")
ylim([0, max([minP.minP_webs])*5])

sgtitle("Failure Loads")
%% 7.4 Plotting Failure Loads on the same graph
pfailfig = figure;
hold on; grid on; grid minor;
Pf1 = plot(x, bridge_properties.P_flex_tens.', "r-", 'LineWidth', 2);
% mark the minimum failure load
plot(x(bridge_properties.P_flex_tens == minP.minP_flex_tens), minP.minP_flex_tens, 'rx')
% label the minimum failure load
for i = 1:length(x)
    if bridge_properties.P_flex_tens(i) == minP.minP_flex_tens
        text(x(i), minP.minP_flex_tens, sigfig(minP.minP_flex_tens))
        break
    end
end

Pf2 = plot(x, bridge_properties.P_flex_comp.', "g-", 'LineWidth', 2);
plot(x(bridge_properties.P_flex_comp == minP.minP_flex_comp), minP.minP_flex_comp, 'rx')
for i = 1:length(x)
    if bridge_properties.P_flex_comp(i) == minP.minP_flex_comp
        text(x(i), minP.minP_flex_comp, sigfig(minP.minP_flex_comp))
        break
    end
end

Pf3 = plot(x, abs(bridge_properties.P_shear.'), "r:", 'LineWidth', 2);
plot(x(bridge_properties.P_shear == minP.minP_shear), minP.minP_shear, 'rx')
for i = 1:length(x)
    if bridge_properties.P_shear(i) == minP.minP_shear
        text(x(i), minP.minP_shear, sigfig(minP.minP_shear))
        break
    end
end

Pf4 = plot(x, abs(bridge_properties.P_glue.'), "g:", 'LineWidth', 2);
plot(x(bridge_properties.P_glue == minP.minP_glue), minP.minP_glue, 'rx')
for i = 1:length(x)
    if bridge_properties.P_glue(i) == minP.minP_glue
        text(x(i), minP.minP_glue, sigfig(minP.minP_glue))
        break
    end
end

Pf5 = plot(x, bridge_properties.P_flange_webs.', "r-.", 'LineWidth', 2);
plot(x(bridge_properties.P_flange_webs == minP.minP_flange_webs), minP.minP_flange_webs, 'rx')
for i = 1:length(x)
    if bridge_properties.P_flange_webs(i) == minP.minP_flange_webs
        text(x(i), minP.minP_flange_webs, sigfig(minP.minP_flange_webs))
        break
    end
end

Pf6 = plot(x, bridge_properties.P_flange_tips.', "g-.", 'LineWidth', 2);
plot(x(bridge_properties.P_flange_tips == minP.minP_flange_tips), minP.minP_flange_tips, 'rx')
for i = 1:length(x)
    if bridge_properties.P_flange_tips(i) == minP.minP_flange_tips
        text(x(i), minP.minP_flange_tips, sigfig(minP.minP_flange_tips))
        break
    end
end

Pf7 = plot(x, bridge_properties.P_webs.', "b-.", 'LineWidth', 2);
plot(x(bridge_properties.P_webs == minP.minP_webs), minP.minP_webs, 'rx')
for i = 1:length(x)
    if bridge_properties.P_webs(i) == minP.minP_webs
        text(x(i), minP.minP_webs, sigfig(minP.minP_webs))
        break
    end
end

Pf8 = plot(x, bridge_properties.P_webs_shear.', "b:", 'LineWidth', 2);
plot(x(bridge_properties.P_webs_shear == minP.minP_webs_shear), minP.minP_webs_shear, 'rx')
for i = 1:length(x)
    if bridge_properties.P_webs_shear(i) == minP.minP_webs_shear
        text(x(i), minP.minP_webs_shear, sigfig(minP.minP_webs_shear))
        break
    end
end

% plot the overall minimum failure load
Pf9 = plot(x, ones(1, n+1)*minP.overall_minP, 'r--', 'LineWidth', 4);
% label the overall minimum failure load
text(x(1), minP.overall_minP, sigfig(minP.overall_minP))

plot(x, zeros(1, n+1), "k", 'LineWidth', 4)

title("Failure Loads")
xlabel("Location on Bridge (mm)")
ylabel("Failure Load (N)")
legend([Pf1, Pf2, Pf3, Pf4, Pf5, Pf6, Pf7, Pf8, Pf9], ["Flexural Tensile Failure of Walls", "Flexural Compressive Failure of Walls", "Shear Failure of Walls", "Glue Shear Failure", "Flange Between Webs Buckling Failure", "Flange Tips Buckling Failure", "Webs Buckling Failure", "Webs Shear Buckling Failure", "Minimum Failure Load"], 'Location', 'best')

ylim([0, minP.overall_minP*10])

set(findall(pfailfig,'-property','FontSize'),'FontSize',16)

%% 8. Bending Moment and Shear Force Capacities
%% 8.1 Calculate the bending moment and shear force capacities using the minimum FoS for each failure mode
bridge_properties.M_tens = abs(min_FOSs.min_FOS_tens.*BMD.');
bridge_properties.M_comp = abs(min_FOSs.min_FOS_comp.*BMD.');
bridge_properties.V_shear = abs(min_FOSs.min_FOS_shear.*SFD.');
bridge_properties.V_glue = abs(min_FOSs.min_FOS_glue.*SFD.');
bridge_properties.M_flange_webs = abs(min_FOSs.min_FOS_flange_webs.*BMD.');
bridge_properties.M_flange_tips = abs(min_FOSs.min_FOS_flange_tips.*BMD.');
bridge_properties.M_webs = abs(min_FOSs.min_FOS_webs.*BMD.');
bridge_properties.V_webs_shear = abs(min_FOSs.min_FOS_webs_shear.*SFD.');

%% 8.2 Plotting the capacities on separate graphs
figure
subplot(2, 3, 1)
hold on; grid on; grid minor;
plot(x, bridge_properties.V_shear, 'r:', 'LineWidth', 2)
title("Shear Capacity of Walls")
xlabel("Location on Bridge (mm)")
ylabel("Shear Capacity (N)")
plot(x, SFD, "k", 'LineWidth', 4)
plot(x, -bridge_properties.V_shear, 'r:', 'LineWidth', 2)
plot(x, zeros(1, n+1), "k", "LineWidth", 4)
legend("Shear Force Capacity", "Shear Force Diagram", "", "")

subplot(2, 3, 2)
hold on; grid on; grid minor;
plot(x, bridge_properties.V_glue, 'g:', 'LineWidth', 2)
title("Glue Shear Capacity")
xlabel("Location on Bridge (mm)")
ylabel("Shear Capacity (N)")
plot(x, SFD, "k", 'LineWidth', 4)
plot(x, -bridge_properties.V_glue, 'g:', 'LineWidth', 2)
plot(x, zeros(1, n+1), "k", "LineWidth", 4)
legend("Shear Force Capacity", "Shear Force Diagram", "", "")

subplot(2, 3, 3)
hold on; grid on; grid minor;
plot(x, bridge_properties.V_webs_shear, 'b:', 'LineWidth', 2)
title("Webs Shear Buckling Capacity")
xlabel("Location on Bridge (mm)")
ylabel("Shear Capacity (N)")
plot(x, SFD, "k", 'LineWidth', 4)
plot(x, -bridge_properties.V_webs_shear, 'b:', 'LineWidth', 2)
plot(x, zeros(1, n+1), "k", "LineWidth", 4)
legend("Shear Force Capacity", "Shear Force Diagram", "", "")

subplot(2, 3, 4)
hold on; grid on; grid minor;
plot(x, bridge_properties.M_tens, 'r', 'LineWidth', 2)
plot(x, bridge_properties.M_comp, 'g', 'LineWidth', 2)
title("Flexural Capacity of Walls")
xlabel("Location on Bridge (mm)")
ylabel("Bending Moment Capacity (N*mm)")
plot(x, BMD, "k", 'LineWidth', 4)
plot(x, zeros(1, n+1), "k", "LineWidth", 4)
legend("Flexural Tensile Capacity of Walls", "Flexural Compressive Capacity of Walls", "Bending Moment Diagram", "")

subplot(2, 3, 5)
hold on; grid on; grid minor;
plot(x, bridge_properties.M_flange_webs, 'r-.', 'LineWidth', 2)
plot(x, bridge_properties.M_flange_tips, 'g-.', 'LineWidth', 2)
title("Flange Buckling Capacity")
xlabel("Location on Bridge (mm)")
ylabel("Bending Moment Capacity (N*mm)")
plot(x, BMD, "k", 'LineWidth', 4)
plot(x, zeros(1, n+1), "k", "LineWidth", 4)
legend("Flange Between Webs Buckling Capacity", "Flange Tips Buckling Capacity", "Bending Moment Diagram", "")

subplot(2, 3, 6)
hold on; grid on; grid minor;
plot(x, bridge_properties.M_webs, 'b-.', 'LineWidth', 2)
title("Webs Buckling Capacity")
xlabel("Location on Bridge (mm)")
ylabel("Bending Moment Capacity (N*mm)")
plot(x, BMD, "k", 'LineWidth', 4)
plot(x, zeros(1, n+1), "k", "LineWidth", 4)
legend("Webs Buckling Capacity", "Bending Moment Diagram", "")

sgtitle("Bending Moment and Shear Force Capacities")

%% 8.3 Plotting the Shear Force Capacilities together
shearcapacityfig = figure;
hold on; grid on; grid minor;
Vc1 = plot(x, bridge_properties.V_shear.', "r-", 'LineWidth', 2);

Vc2 = plot(x, abs(bridge_properties.V_glue.'), "g-", 'LineWidth', 2);

Vc3 = plot(x, bridge_properties.V_webs_shear.', "b-", 'LineWidth', 2);

% plot the SFD as a reference
Vc4 = plot(x, SFD, "k", 'LineWidth', 4);

plot(x, zeros(1, n+1), "k", "LineWidth", 4)

title("Shear Force Capacities")
xlabel("Location on Bridge (mm)")
ylabel("Shear Force Capacity (N)")
legend([Vc1, Vc2, Vc3, Vc4], ["Shear Failure of Walls", "Glue Shear Failure", "Webs Shear Buckling Failure", "Shear Force Envelope"], 'Location', 'best')

set(findall(shearcapacityfig,'-property','FontSize'),'FontSize',16)

%% 8.4 Plotting the Bending Moment Capacilities together
momentcapacityfig = figure;
hold on; grid on; grid minor;
Mc1 = plot(x, bridge_properties.M_tens.', "r-", 'LineWidth', 2);

Mc2 = plot(x, bridge_properties.M_comp.', "g-", 'LineWidth', 2);

Mc3 = plot(x, bridge_properties.M_flange_webs.', "r-.", 'LineWidth', 2);

Mc4 = plot(x, bridge_properties.M_flange_tips.', "g-.", 'LineWidth', 2);

Mc5 = plot(x, bridge_properties.M_webs.', "b-.", 'LineWidth', 2);

% plot the BMD as a reference
Mc6 = plot(x, BMD, "k", 'LineWidth', 4);

plot(x, zeros(1, n+1), "k", "LineWidth", 4)

title("Bending Moment Capacities")
xlabel("Location on Bridge (mm)")
ylabel("Bending Moment Capacity (N*mm)")
legend([Mc1, Mc2, Mc3, Mc4, Mc5, Mc6], ["Flexural Tensile Failure of Walls", "Flexural Compressive Failure of Walls", "Flange Between Webs Buckling Failure", "Flange Tips Buckling Failure", "Webs Buckling Failure", "Bending Moment Diagram"], 'Location', 'best')

set(findall(momentcapacityfig,'-property','FontSize'),'FontSize',16)

%% 9. Functions
%% 9.1 Reaction Solver
function [R1, R2] = reaction(P, x)
% REACTION solve for reaction forces at the ends of the bridge
% the supports are at x = 0 and x = 1200
% P are the loads on the bridge
% P = [P1, P2, ...]
% x are the x locations of the loads
% x = [x1, x2, ...]
% R1, R2 are the reaction forces at the ends of the bridge
% net force in y direction = 0 -> R1 + R2 + sum(P) = 0 -> R1 + R2 = -sum(P)
% net moment about R1 = 0 -> R2*1200 + sum(P*x) = 0 -> R2 = -sum(P*x)/1200
    % solve for R2
    % multiply each load by its x location
    Mx = -(P.*x);
    % sum the moments
    M = sum(Mx);
    % solve for R2
    R2 = M/1200;

    % solve for R1
    R1 = -sum(P) - R2;
end

%% 9.2 Internal Forces Solver
function [shear, moment] = internal_forces(P, x, cut)
    % INTERNAL_FORCES solve for the internal forces at each location by cutting the bridge at the cut location
    % P are the loads that are within the cut
    % P = [P1, P2, ...]
    % x are the x locations of the loads that are within the cut
    % x = [x1, x2, ...]
    % cut is the location of the cut
    % shear is the shear force at the cut
    % moment is the bending moment at the cut

    % solve for the shear force at the cut
    % net force in y direction = 0 -> shear + sum(P) = 0
    shear = sum(P);

    % solve for the bending moment at the cut
    % net moment about the cut = 0 -> moment + sum(P*x) = 0
    % multiply each load by its distance from the cut
    Mx = P.*(x - cut);
    % sum the moments
    moment = -sum(Mx);
end

%% 9.3 Significant Figures Formatter
function output = sigfig(input)
    % SIGFIG format the input to 3 significant figures if the first digit is not 1, and 4 significant figures if the first digit is 1
    % input is the number to be formatted as a number
    % output is the formatted number as a string

    % convert the input to scientific notation
    input = sprintf("%e", input);

    % convert the input to a character array
    input = char(input);

    % check if the first digit is 1
    if input(1) == "1"
        input = string(input);
        % format the input to 4 significant figures
        output = char(sprintf("%.4g", input));
        
    else
        input = string(input);
        % format the input to 3 significant figures
        output = char(sprintf("%.3g", input));
    end

    % shift the decimal point until the exponent is a multiple of 3
    % find the index of the decimal point
    decimal_index = strfind(output, ".");
    % find the index of the e
    e_index = strfind(output, "e");
    % find the exponent
    exponent = str2double(output(e_index+1:end));
    if exponent > 3
        while exponent > 3 && mod(exponent, 3) ~= 0
            % shift the decimal point to the right
            output(decimal_index) = output(decimal_index+1);
            output(decimal_index+1) = ".";
            % update the exponent
            exponent = exponent - 1;

            % update the decimal index
            decimal_index = decimal_index + 1;

            % if the decimal index is at the e index, remove the decimal
            if (decimal_index == e_index)
                output = output(1:e_index-1) + output(e_index+1:end);
                % update the e index
                e_index = e_index - 1;
            end
        end
        % remove the e and the exponent
        output = output(1:e_index-1);

        % format the exponent to 2 characters
        exponent = sprintf("%02d", exponent);

        % add the exponent back
        output = output + "e+" + exponent;
    end
end