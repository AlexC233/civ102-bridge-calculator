clear; clc; close all;
%% 0. Initialize Parameters
L = 1200; % Length of bridge
n = 1200; % Discretize into 1 mm seg.
dL = L/n; % Length of each segment
P = 446.666666; % Total weight of train [N]
x = linspace(0, L, n+1); % x-axis

%% 1. SFD, BMD under train loading
x_train = [52 228 392 568 732 908]; % Train Load Locations (the 6 wheels)
P_factors = [1.35 1.35 1 1 1 1]; % Load Factors (the 6 wheels)
P_train = P_factors .* P/sum(P_factors) * -1; % load of each wheel
n_train = 3; % num of train locations
SFDi = zeros(n+1); % 1 SFD for each train loc.
BMDi = zeros(n+1); % 1 BMD for each train loc.
%% 1.1 Solve for SFD and BMD with the train at different locations

% for every discretized location of the train, find the SFD and BMD
for i = 1:n+1
    R1 = 0; % Reaction 1
    R2 = 0; % Reaction 2
    
    % find the loads that are on the bridge
    x_train_on_bridge = x_train(x_train <= x(end));
    P_train_on_bridge = P_train(x_train <= x(end));

    % if there are loads on the bridge, solve for the reactions
    if ~isempty(x_train_on_bridge)
        % solve for the reactions
        [R1, R2] = reaction(P_train_on_bridge, x_train_on_bridge);
    else
        % if there are no loads on the bridge, the reactions are 0
        R1 = 0;
        R2 = 0;
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

SFD = max(abs(SFDi)); % SFD envelope
BMD = max(BMDi); % BMD envelope

% plot the SFD and BMD of the train centered on the bridge
figure
plot(x, SFDi(121,:), 'r')
xlabel('Distance along bridge (mm)')
ylabel('Shear Force (N)')
title('SFD of Train Centered on Bridge')
figure
plot(x, BMDi(121,:), 'r')
xlabel('Distance along bridge (mm)')
ylabel('Bending Moment (N*mm)')
title('BMD of Train Centered on Bridge')

%% 2. Define Cross Section
% x_sections is stored as a dictionary with the keys being the x location along the bridge and the values being the cross section at that location
% each cross section is stored as an array of subsections with the following values:
% [x, y, dx, dy; ...]
% x is the x location of the reference point of the subsection
% y is the y location of the reference point of the subsection
% dx is the width of the subsection
% dy is the height of the subsection
% the reference point is the bottom left corner of the subsection
% the sections are constructed from the bottom up with 0, 0 being the bottom left corner of the bridge

x_change = [0]; % x locations of cross section changes
x_sections = {[10, 0, 80, 1.27; 
               10, 1.27, 1.27, 75-1.27*2; 
               90-1.27, 1.27, 1.27, 75-2*1.27; 
               10, 75-1.27, 5 + 1.27, 1.27;
               90-5-1.27, 75-1.27, 5+1.27, 1.27; 
               0, 75, 100, 1.27]}; % cross sections of the bridge

% DEBUG draw the first cross section
figure
hold on; grid on; grid minor;
axis equal
axis([0, 100, 0, 100])
for i = 1:size(x_sections{1},1)
    rectangle('Position', x_sections{1}(i,:), 'FaceColor', 'b')
end


params = dictionary(x_change, x_sections); % dictionary of the cross sections

%% 3. Calculate Sectional Properties
% set up arrays to store the sectional properties
A =     zeros(length(x_change), 1); % areas
ybar =  zeros(length(x_change), 1); % location of centroidal axis from the bottom of the cross section
ybot =  zeros(length(x_change), 1); % location of bottom of cross section from the bottom of the cross section
ytop =  zeros(length(x_change), 1); % location of top of cross section from the bottom of the cross section
I =     zeros(length(x_change), 1); % moment of inertia
Qcent = zeros(length(x_change), 1); % Q at centroidal axes
Qglue = zeros(length(x_change), 1); % Q at glue location
b =     zeros(length(x_change), 1); % width of cross section that crosses the centroidal axis

% for each cross section, calculate the sectional properties
for i = 1:length(x_change)
    % find the cross section
    x_section = params(x_change(i));
    x_section = x_section{1, 1};

    % setup a table for the cross section
    x_section = array2table(x_section, 'VariableNames', {'x', 'y', 'dx', 'dy'});

    % add a column for the area of each subsection
    x_section.area = x_section.dx.*x_section.dy;

    % add a column for the centroidal axis of each subsection
    x_section.ybar = x_section.y + x_section.dy/2;

    % add a column for the moment of inertia of each subsection
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

    % find the moment of inertia
    % using parallel axis theorem
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

    % DEBUG plot these subsections
    for j = 1:size(x_section_cent_bot,1)
        rectangle('Position', [x_section_cent_bot.x(j), x_section_cent_bot.ybot(j), x_section_cent_bot.dx(j), x_section_cent_bot.dy(j)], 'FaceColor', 'r')
    end

    % calculate the centroidal axis of the subsections combined
    ybar_cent_bot = sum(x_section_cent_bot.area.*x_section_cent_bot.ybar)/sum(x_section_cent_bot.area);

    % find Q at the centroidal axis
    Qcent(i) = sum(x_section_cent_bot.area.*(ybar(i) - ybar_cent_bot));

    % find Q at the glue location
    % Qglue(i) = sum(x_section(:,3).*(x_section(:,2) - ybot(i)).*x_section(:,4).^2);

end

% setup a table to store these properties for every discretized location
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
    %bridge_properties.Qglue(x_in_section) = Qglue(i);
    bridge_properties.b(x_in_section) = b(i);
end


%% 4. Calculate Applied Stress
% stress at the top across the entire bridge
bridge_properties.sigma_top = BMD.*abs(bridge_properties.ytop - bridge_properties.ybar)./bridge_properties.I;

% stress at the bottom across the entire bridge
bridge_properties.sigma_bot = BMD.*abs(bridge_properties.ybot - bridge_properties.ybar)./bridge_properties.I;

% shear stress across the entire bridge at the centroidal axis
bridge_properties.tau_xy = SFD.*bridge_properties.Qcent./bridge_properties.I./bridge_properties.b;

% T_glue =

%% 5. Material and Thin Plate Buckling Capacities
% setup a table to store the capacities of the material
material_properties = array2table(zeros(1, 10), "VariableNames", ["E", "mu", "sigma_tens", "sigma_comp", "tau_max", "tau_gmax", "sigma_buck_flange_webs", "sigma_buck_flange_tips", "sigma_buck_webs", "tau_buck_webs"]);

material_properties.E = 4000;
material_properties.mu = 0.2;

material_properties.sigma_tens = 30;
material_properties.sigma_comp = 6;
material_properties.tau_max = 4;

material_properties.tau_gmax = 2;

% S_buck1 =
% S_buck2 =
% S_buck3 =
% T_buck =
%% 6. FOS
bridge_properties.FOS_tens = material_properties.sigma_tens./bridge_properties.sigma_top;
bridge_properties.FOS_comp = material_properties.sigma_comp./bridge_properties.sigma_bot;
bridge_properties.FOS_shear = material_properties.tau_max./bridge_properties.tau_xy;

% DELIVERABLE 1 for train centered on bridge
FOS_tens = material_properties.sigma_tens(1)/(max(BMDi(121,:))*abs(bridge_properties.ybot(1) - bridge_properties.ybar(1))/bridge_properties.I(1));
FOS_comp = material_properties.sigma_comp(1)/(max(BMDi(121,:))*abs(bridge_properties.ytop(1) - bridge_properties.ybar(1))/bridge_properties.I(1));

% bridge_properties.FOS_glue =
% bridge_properties.FOS_buck1 =
% bridge_properties.FOS_buck2 =
% bridge_properties.FOS_buck3 =
% bridge_properties.FOS_buckV =
% %% 7. Min FOS and the failure load Pfail
% minFOS =
% Pf =
% %% 8. Vfail and Mfail
% Mf_tens =
% Mf_comp =
% Vf_shear =
% Vf_glue =
% Mf_buck1 =
% Mf_buck2 =
% Mf_buck3 =
% Vf_buckV =
% %% 9. Output plots of Vfail and Mfail
% subplot(2,3,1)
% hold on; grid on; grid minor;
% plot(x, Vf_shear, 'r')
% plot(x, -Vf_shear.* SFD, 'r')
% plot(x, SFDi, 'k');
% plot([0, L], [0, 0], 'k', 'LineWidth', 2)
% legend('Matboard Shear Failure')
% xlabel('Distance along bridge (mm)')
% ylabel('Shear Force (N)')

%% Functions
% solve for reaction forces
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