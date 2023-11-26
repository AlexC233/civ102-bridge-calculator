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
P_factors = [1.35 1.35 1 1 1 1]; % Load Factors (the 6 wheels)
P_train = P_factors .* P/sum(P_factors) * -1; % load of each wheel

%% 1.2 Solve for SFD and BMD with the train at different locations
% set up the SFD and BMD arrays of the bridge for every discretized location of the train
% the rows are the locations of the train, with the first row being the train completely off of the bridge on the left and the last row being the train completely off of the bridge on the right
% the columns are the forces/moments at each cut, with the first column being the left end of the bridge and the last column being the right end of the bridge
SFDi = zeros(n + 1 + l_train, n + 1); % SFDs of the bridge
BMDi = zeros(n + 1 + l_train, n + 1); % BMDs of the bridge

% shift the train entirely off of the bridge to the left
x_train = x_train - L;

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

%% 1.2 Plot SFD and BMD
SFD = max(abs(SFDi)); % SFD envelope
BMD = max(BMDi); % BMD envelope

% plot the SFD and BMD envelopes on separate figures
% SFD is red
% BMD is blue
figure
hold on
plot(x, SFD, 'r')
%plot(x, -SFD, 'r')
plot(x, zeros(1, n+1), 'k')

title("Shear Force Envelope")
xlabel("Location on Bridge (mm)")
ylabel("Shear Force (N)")
legend('Shear Force')

figure
hold on
% invert the y axis so that positive bending moments are plotted downwards
set(gca, 'YDir','reverse')
plot(x, BMD, 'b')
%plot(x, -BMD, 'b')
plot(x, zeros(1, n+1), 'k')

title("Bending Moment Envelope")
xlabel("Location on Bridge (mm)")
ylabel("Bending Moment (N*mm)")
legend('Bending Moment')

%% 2. Define Cross Section
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

%% Glue
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

%% Design 0
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

x_section_params = dictionary(x_change, x_sections); % dictionary of the cross sections
glue_params = dictionary(x_change, glue_locations); % dictionary of the glue locations

%% 2.2 Plotting the Cross Sections
for i = 1:length(x_change)
    % find the cross section
    x_section = x_section_params(x_change(i));
    x_section = x_section{1, 1};

    % find the glue locations
    glue = glue_params(x_change(i));
    glue = glue{1, 1};

    % plot the cross section
    figure
    hold on
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
    title(['Cross Section at x = ', num2str(x_change(i)), ' mm'])
end

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

    % DEBUG add the ybar to the figure
    plot([0, 100], [ybar(i), ybar(i)], 'k', 'LineWidth', 2)

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
    glue = glue{1, 1}(glue{1, 1}(:,5) == 1, :);

    % sum the widths of the glue
    gL(i) = sum(glue(:,4));
    
    % if there are glue locations, calculate Q at the glue locations

    % check if the glue is horizontal or vertical
    % assuming that all glue for a given cross section is either all horizontal or all vertical
    if ~isempty(glue)
        bg = 0; % width of the glue
        if glue(1, 1) == 0
            % only calculate the glue if calculate is 1
                % if the glue is horizontal, find the width of the glue
                % find the subsections that are within the glue
                x_section_glue = x_section(x_section.x >= glue(1,2) & x_section.x <= glue(1,2) + glue(1,4), :);

                % find the width of the glue
                bg = sum(glue(:,4));

                % find the centroidal axis of the glue
                ybar_glue = sum(x_section_glue.area.*x_section_glue.ybar)/sum(x_section_glue.area);

                % find Q at the glue location
                Qglue(i) = sum(x_section_glue.area.*(ybar(i) - ybar_glue));
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
end

%% 3.4 Calculate the amount of material needed
matboard = sum(bridge_properties.A);
disp("Matboard needed: " + matboard + " mm^3")

%% 3.5 Calculate the amount of glue needed
glue_needed = sum(bridge_properties.gL);
disp("Glue needed: " + glue_needed + " mm^2")

%% 4. Calculate Applied Stress
% stress at the top across the entire bridge
bridge_properties.sigma_top = BMD.'.*abs(bridge_properties.ytop - bridge_properties.ybar)./bridge_properties.I;

% stress at the bottom across the entire bridge
bridge_properties.sigma_bot = BMD.'.*abs(bridge_properties.ybot - bridge_properties.ybar)./bridge_properties.I;

% shear stress across the entire bridge at the centroidal axis
bridge_properties.tau_xy = SFD.'.*bridge_properties.Qcent./bridge_properties.I./bridge_properties.b;

% shear stress across the entire bridge at the glue locations
bridge_properties.tau_g = SFD.'.*bridge_properties.Qglue./bridge_properties.I./bridge_properties.b;

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
disp(min_FOSs)

%% 6.3 Plotting the FOSs
figure
hold on
fos1 = plot(x, bridge_properties.FOS_tens.', 'r-');
% display the minimum FOS as an x
plot(x(bridge_properties.FOS_tens == min_FOSs.min_FOS_tens), min_FOSs.min_FOS_tens, 'rx')
% label the minimum FOS
text(x(bridge_properties.FOS_tens == min_FOSs.min_FOS_tens), min_FOSs.min_FOS_tens, num2str(min_FOSs.min_FOS_tens))

fos2 = plot(x, bridge_properties.FOS_comp.', 'g-');
plot(x(bridge_properties.FOS_comp == min_FOSs.min_FOS_comp), min_FOSs.min_FOS_comp, 'rx')
text(x(bridge_properties.FOS_comp == min_FOSs.min_FOS_comp), min_FOSs.min_FOS_comp, num2str(min_FOSs.min_FOS_comp))

fos3 = plot(x, abs(bridge_properties.FOS_shear.'), 'r:');
plot(x(bridge_properties.FOS_shear == min_FOSs.min_FOS_shear), min_FOSs.min_FOS_shear, 'rx')
text(x(bridge_properties.FOS_shear == min_FOSs.min_FOS_shear), min_FOSs.min_FOS_shear, num2str(min_FOSs.min_FOS_shear))

fos4 = plot(x, abs(bridge_properties.FOS_glue.'), 'g:');
plot(x(bridge_properties.FOS_glue == min_FOSs.min_FOS_glue), min_FOSs.min_FOS_glue, 'rx')
text(x(bridge_properties.FOS_glue == min_FOSs.min_FOS_glue), min_FOSs.min_FOS_glue, num2str(min_FOSs.min_FOS_glue))

fos5 = plot(x, bridge_properties.FOS_flange_webs.', 'r-.');
plot(x(bridge_properties.FOS_flange_webs == min_FOSs.min_FOS_flange_webs), min_FOSs.min_FOS_flange_webs, 'rx')
text(x(bridge_properties.FOS_flange_webs == min_FOSs.min_FOS_flange_webs), min_FOSs.min_FOS_flange_webs, num2str(min_FOSs.min_FOS_flange_webs))

fos6 = plot(x, bridge_properties.FOS_flange_tips.', 'g-.');
plot(x(bridge_properties.FOS_flange_tips == min_FOSs.min_FOS_flange_tips), min_FOSs.min_FOS_flange_tips, 'rx')
text(x(bridge_properties.FOS_flange_tips == min_FOSs.min_FOS_flange_tips), min_FOSs.min_FOS_flange_tips, num2str(min_FOSs.min_FOS_flange_tips))

fos7 = plot(x, bridge_properties.FOS_webs.', 'b-.');
plot(x(bridge_properties.FOS_webs == min_FOSs.min_FOS_webs), min_FOSs.min_FOS_webs, 'rx')
text(x(bridge_properties.FOS_webs == min_FOSs.min_FOS_webs), min_FOSs.min_FOS_webs, num2str(min_FOSs.min_FOS_webs))

fos8 = plot(x, abs(bridge_properties.FOS_webs_shear.'), 'b:');
plot(x(bridge_properties.FOS_webs_shear == min_FOSs.min_FOS_webs_shear), min_FOSs.min_FOS_webs_shear, 'rx')
text(x(bridge_properties.FOS_webs_shear == min_FOSs.min_FOS_webs_shear), min_FOSs.min_FOS_webs_shear, num2str(min_FOSs.min_FOS_webs_shear))

% plot the overall minimum FOS
fos9 = plot(x, min_FOSs.overall_min_FOS*ones(n+1, 1), 'r--');
text(x(1), min_FOSs.overall_min_FOS, num2str(min_FOSs.overall_min_FOS))

% plot the FOS = 1 line
fos10 = plot(x, ones(n+1, 1), 'k--');

title("Factors of Safety")
xlabel("Location on Bridge (mm)")
ylabel("Factor of Safety")
legend([fos1, fos2, fos3, fos4, fos5, fos6, fos7, fos8, fos9, fos10], ["Tensile Failure of Walls", "Compressive Failure of Walls", "Shear Failure of Walls", "Glue Shear Failure", "Flange Between Webs Buckling Failure", "Flange Tips Buckling Failure", "Webs Buckling Failure", "Webs Shear Buckling Failure", "Minimum FOS", "FOS = 1"], 'Location', 'best')

% resize to only show FOS up to 10
ylim([0, 10])

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
disp(minP)

%% 7.3 Plotting Failure Loads
figure
hold on
Pf1 = plot(x, bridge_properties.P_flex_tens.', "r-");
% mark the minimum failure load
plot(x(bridge_properties.P_flex_tens == minP.minP_flex_tens), minP.minP_flex_tens, 'rx')
% label the minimum failure load
text(x(bridge_properties.P_flex_tens == minP.minP_flex_tens), minP.minP_flex_tens, num2str(minP.minP_flex_tens))

Pf2 = plot(x, bridge_properties.P_flex_comp.', "g-");
plot(x(bridge_properties.P_flex_comp == minP.minP_flex_comp), minP.minP_flex_comp, 'rx')
text(x(bridge_properties.P_flex_comp == minP.minP_flex_comp), minP.minP_flex_comp, num2str(minP.minP_flex_comp))    

Pf3 = plot(x, abs(bridge_properties.P_shear.'), "r:");
plot(x(bridge_properties.P_shear == minP.minP_shear), minP.minP_shear, 'rx')
text(x(bridge_properties.P_shear == minP.minP_shear), minP.minP_shear, num2str(minP.minP_shear))

Pf4 = plot(x, abs(bridge_properties.P_glue.'), "g:");
plot(x(bridge_properties.P_glue == minP.minP_glue), minP.minP_glue, 'rx')
text(x(bridge_properties.P_glue == minP.minP_glue), minP.minP_glue, num2str(minP.minP_glue))

Pf5 = plot(x, bridge_properties.P_flange_webs.', "r-.");
plot(x(bridge_properties.P_flange_webs == minP.minP_flange_webs), minP.minP_flange_webs, 'rx')
text(x(bridge_properties.P_flange_webs == minP.minP_flange_webs), minP.minP_flange_webs, num2str(minP.minP_flange_webs))

Pf6 = plot(x, bridge_properties.P_flange_tips.', "g-.");
plot(x(bridge_properties.P_flange_tips == minP.minP_flange_tips), minP.minP_flange_tips, 'rx')
text(x(bridge_properties.P_flange_tips == minP.minP_flange_tips), minP.minP_flange_tips, num2str(minP.minP_flange_tips))

Pf7 = plot(x, bridge_properties.P_webs.', "b-.");
plot(x(bridge_properties.P_webs == minP.minP_webs), minP.minP_webs, 'rx')
text(x(bridge_properties.P_webs == minP.minP_webs), minP.minP_webs, num2str(minP.minP_webs))

Pf8 = plot(x, bridge_properties.P_webs_shear.', "b:");
plot(x(bridge_properties.P_webs_shear == minP.minP_webs_shear), minP.minP_webs_shear, 'rx')
text(x(bridge_properties.P_webs_shear == minP.minP_webs_shear), minP.minP_webs_shear, num2str(minP.minP_webs_shear))

% plot the overall minimum failure load
Pf9 = plot(x, ones(1, n+1)*minP.overall_minP, 'k--');
% label the overall minimum failure load
text(x(1), minP.overall_minP, num2str(minP.overall_minP))

plot(x, zeros(1, n+1), "k")

title("Failure Loads")
xlabel("Location on Bridge (mm)")
ylabel("Failure Load (N)")
legend([Pf1, Pf2, Pf3, Pf4, Pf5, Pf6, Pf7, Pf8, Pf9], ["Flexural Tensile Failure of Walls", "Flexural Compressive Failure of Walls", "Shear Failure of Walls", "Glue Shear Failure", "Flange Between Webs Buckling Failure", "Flange Tips Buckling Failure", "Webs Buckling Failure", "Webs Shear Buckling Failure", "Minimum Failure Load"], 'Location', 'best')

% only display failure loads up to 2000 N
ylim([0, 2000])

% %% 8. Bending Moment and Shear Force Capacities
% %% 8.1 Calculate the bending moment and shear force capacities using the FOS and the BMD and SFD
% bridge_properties.M_tens = abs(bridge_properties.FOS_tens.*BMD.');
% bridge_properties.M_comp = abs(bridge_properties.FOS_comp.*BMD.');
% bridge_properties.V_shear = abs(bridge_properties.FOS_shear.*SFD.');
% bridge_properties.V_glue = abs(bridge_properties.FOS_glue.*SFD.');
% bridge_properties.M_flange_webs = abs(bridge_properties.FOS_flange_webs.*BMD.');
% bridge_properties.M_flange_tips = abs(bridge_properties.FOS_flange_tips.*BMD.');
% bridge_properties.M_webs = abs(bridge_properties.FOS_webs.*BMD.');
% bridge_properties.V_webs_shear = abs(bridge_properties.FOS_webs_shear.*SFD.');

% %% 8.2 Find the minimum bending moment and shear force capacities
% minM = array2table(zeros(1, 9), 'VariableNames', ["minM_tens", "minM_comp", "minV_shear", "minV_glue", "minM_flange_webs", "minM_flange_tips", "minM_webs", "minV_webs_shear", "overall_minM"]);
% minM.minM_tens = min(bridge_properties.M_tens);
% minM.minM_comp = min(bridge_properties.M_comp);
% minM.minV_shear = min(bridge_properties.V_shear);
% minM.minV_glue = min(bridge_properties.V_glue);
% minM.minM_flange_webs = min(bridge_properties.M_flange_webs);
% minM.minM_flange_tips = min(bridge_properties.M_flange_tips);
% minM.minM_webs = min(bridge_properties.M_webs);
% minM.minV_webs_shear = min(bridge_properties.V_webs_shear);
% minM.overall_minM = min([minM.minM_tens, minM.minM_comp, minM.minV_shear, minM.minV_glue, minM.minM_flange_webs, minM.minM_flange_tips, minM.minM_webs, minM.minV_webs_shear]);
% disp(minM)

% %% 8.3 Plotting the Shear Force Capacilities
% figure
% hold on
% Vc1 = plot(x, bridge_properties.V_shear.', "r-");
% % mark the minimum shear force capacity
% plot(x(bridge_properties.V_shear == minM.minV_shear), minM.minV_shear, 'rx')
% % label the minimum shear force capacity
% % text(x(bridge_properties.V_shear == minM.minV_shear), minM.minV_shear, num2str(minM.minV_shear))

% Vc2 = plot(x, abs(bridge_properties.V_glue.'), "g-");
% plot(x(bridge_properties.V_glue == minM.minV_glue), minM.minV_glue, 'rx')
% % text(x(bridge_properties.V_glue == minM.minV_glue), minM.minV_glue, num2str(minM.minV_glue))

% Vc3 = plot(x, bridge_properties.V_webs_shear.', "b-");
% plot(x(bridge_properties.V_webs_shear == minM.minV_webs_shear), minM.minV_webs_shear, 'rx')
% % text(x(bridge_properties.V_webs_shear == minM.minV_webs_shear), minM.minV_webs_shear, num2str(minM.minV_webs_shear))

% % plot the overall minimum shear force capacity
% Vc4 = plot(x, ones(1, n+1)*minM.overall_minM, 'k--');
% % label the overall minimum shear force capacity
% % text(x(1), minM.overall_minM, num2str(minM.overall_minM))

% % plot the SFD as a reference
% Vc5 = plot(x, SFD, "k");

% title("Shear Force Capacities")
% xlabel("Location on Bridge (mm)")
% ylabel("Shear Force Capacity (N)")
% legend([Vc1, Vc2, Vc3, Vc4, Vc5], ["Shear Failure of Walls", "Glue Shear Failure", "Webs Shear Buckling Failure", "Minimum Shear Force Capacity", "Shear Force Diagram"], 'Location', 'best')

% % only display shear force capacities up to 2000 N
% ylim([0, 2000])

% %% 8.4 Plotting the Bending Moment Capacilities
% figure
% hold on
% Mc1 = plot(x, bridge_properties.M_tens.', "r-");
% % % mark the minimum bending moment capacity
% % plot(x(bridge_properties.M_tens == minM.minM_tens), minM.minM_tens, 'rx')
% % % label the minimum bending moment capacity
% % text(x(bridge_properties.M_tens == minM.minM_tens), minM.minM_tens, num2str(minM.minM_tens))

% Mc2 = plot(x, bridge_properties.M_comp.', "g-");
% % plot(x(bridge_properties.M_comp == minM.minM_comp), minM.minM_comp, 'rx')
% % text(x(bridge_properties.M_comp == minM.minM_comp), minM.minM_comp, num2str(minM.minM_comp))

% Mc3 = plot(x, bridge_properties.M_flange_webs.', "r-.");
% % plot(x(bridge_properties.M_flange_webs == minM.minM_flange_webs), minM.minM_flange_webs, 'rx')
% % text(x(bridge_properties.M_flange_webs == minM.minM_flange_webs), minM.minM_flange_webs, num2str(minM.minM_flange_webs))

% Mc4 = plot(x, bridge_properties.M_flange_tips.', "g-.");
% % plot(x(bridge_properties.M_flange_tips == minM.minM_flange_tips), minM.minM_flange_tips, 'rx')
% % text(x(bridge_properties.M_flange_tips == minM.minM_flange_tips), minM.minM_flange_tips, num2str(minM.minM_flange_tips))

% Mc5 = plot(x, bridge_properties.M_webs.', "b-.");
% % plot(x(bridge_properties.M_webs == minM.minM_webs), minM.minM_webs, 'rx')
% % text(x(bridge_properties.M_webs == minM.minM_webs), minM.minM_webs, num2str(minM.minM_webs))

% % plot the overall minimum bending moment capacity
% Mc6 = plot(x, ones(1, n+1)*minM.overall_minM, 'k--');
% % label the overall minimum bending moment capacity
% % text(x(1), minM.overall_minM, num2str(minM.overall_minM))

% % plot the BMD as a reference
% Mc7 = plot(x, BMD, "k");

% title("Bending Moment Capacities")
% xlabel("Location on Bridge (mm)")
% ylabel("Bending Moment Capacity (N*mm)")
% legend([Mc1, Mc2, Mc3, Mc4, Mc5, Mc6, Mc7], ["Flexural Tensile Failure of Walls", "Flexural Compressive Failure of Walls", "Flange Between Webs Buckling Failure", "Flange Tips Buckling Failure", "Webs Buckling Failure", "Minimum Bending Moment Capacity", "Bending Moment Diagram"], 'Location', 'best')

% % only display bending moment capacities up to 10e5 N*mm
% ylim([0, 10e5])

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