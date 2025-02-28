


%% ________________________________Task 1__________________________________
%  Implement a script in MATLAB which generates the reference currents isd(Te,ωs), isq(Te,ωs) 
% for a specific operating point Te,ωs according to the principles taught in the lecture to operate a 
% generic permanent magnet synchronous machine safe and efficient in base speed operation, basic field weakening
%  and field weakening. Please use the following machine and system parameters as well as the respective naming conventions in your script.

% Set Reference Values by asking the user for input:
disp("%%%% Short-Circuit current is within the MA-Circle (kappa<1) %%%%% ") 
disp("%%%% Example: omega=10000rps and Te=70Nm %%%%% ")
disp("%%%% Press str+C to end the Input Loop %%%%% ")
prompt1 = 'Enter the speed omega_s in rps (range from 1-20943): ';
prompt2 = 'Enter the Torque in Nm (1-100): ';
omega_s =input(prompt1);  %rps
T_e    = input(prompt2);  %Electrical torque in Nm

is_dq_values = isdq_components(omega_s,T_e); %calling the function isdq_components(omega_s,T_e)

% ________________________________END_Task 1_______________________________



disp("Task 2/3 is loading...approximate: 1 Minute... ")
%% ________________________________Task 2 / Task 3_________________________
% -----------------------------> Task 2: i_smax = 500;---------------------
% -----------------------------> Task 3: i_smax = 400;---------------------
% ----------------------> Please change the values manually----------------
% _________________________________________________________________________ 
%% initialize
%___________________switch for Task 2/ Task 3______________________________
i_smax= 500;        %Maximum stator current (amplitude) [A]
%__________________________________________________________________________

Psi_f = 90e-3;        %Field flux linkage [Vs]
L_sd  = 200e-6;       %Stator inductance in d-axis [H]
L_sq  = 500e-6;       %Stator inductance in q-axis [H]
U_dc  = 350;          %Dc-link voltage [V]
u_smax= U_dc/sqrt(3); %Maximum stator voltage
P     = 4;            %Pole pairs

i_sc   = 450;                    %Short circuit current
kappa  = i_sc/i_smax;            %Short circuit current normalized on ismax
chi    = (L_sq-L_sd)/(2*L_sd);   %Saliency


%% Set Reference Values:


T_e    = 70;  %Electrical torque in Nm
speed_max = linspace(1, 50000, 1000);    % in rpm von 1-50000
omega_s = 2*pi.*speed_max./60.*P;        % in rps for formulas, consider P
omega_s_ellipse_set = 20943.95;          %sets the omega_s for the desired MF-Ellipse, changes the desired plotted ellipse

%% Create Environment


if kappa<1  % Differentiate between kappa<1 and kappa >1

 %% MA-Circle normalized to i_smax in parametric form

        theta = linspace(pi/2,pi,100);  
        
        a_ma = 0;
        b_ma = 0;
        r_ma = 1;
        
        i_sd_ma = a_ma + r_ma.*cos(theta);
        i_sq_ma = b_ma + r_ma.*sin(theta); 
        
        % interpolate MA-Circle for intersection points
        splines_ma = spline(i_sd_ma, i_sq_ma);
        x_interpolated_ma = linspace(min(i_sd_ma), max(i_sd_ma), 1000);
        y_interpolated_ma = ppval(splines_ma, x_interpolated_ma);
        
    %% MTPA Line
       
        i_s = linspace(0, 1, 15);                                                % Generate values for i_s
       
        i_sd_mtpa = (kappa/(8*chi)) - sqrt((i_s.^2/2 + (kappa/(8*chi))^2));      % Calculate i_sd and i_sq with formulas from the lecture
        i_sq_mtpa = sqrt(i_s.^2 - i_sd_mtpa.^2);
        
        % interpolate MTPA-Line for intersection points
        interp_func_mtpa = spline(i_sd_mtpa, i_sq_mtpa);
        evaluation_points_mtpa = linspace(-1, 0, 100);  
        interpolated_values_mtpa = ppval(interp_func_mtpa, evaluation_points_mtpa);
        
    %% MTPF Line

        speed_max1 = linspace(1, 1000000, 1000); % in rpm, increase the resolution to increase interpolation-quality
        omega_s_mtpf = 2*pi.*speed_max1./60.*P;  % in rps

        % Formulas from the lecture:
        i_o = u_smax./(omega_s_mtpf.*L_sd.*i_smax);
        i_sd_mtpf = -kappa + (1 + 2*chi).*kappa./(8*chi) - sqrt((i_o/2).^2 + ((1 + 2*chi)*kappa/(8*chi)).^2);
        i_sq_mtpf = 1./(2*chi+1).*sqrt(i_o.^2 - (kappa + i_sd_mtpf).^2);
        
       
        % interpolate MTPF-Line for intersection points
        x_interpolated_mtpf = linspace(-2, -i_sc/i_smax, 1000);
        splines_mtpf = spline(i_sd_mtpf, i_sq_mtpf);
        y_interpolated_mtpf = ppval(splines_mtpf, x_interpolated_mtpf);


        
    %% Constant Torque line
        
        i_sd_CT =  linspace(-i_smax, 0, 100); % Define the range of i_sd values for the 2nd quadrant
        
        i_sq_CT = -2/3 * T_e/Psi_f * (1 ./ (2 * chi / i_sc .* i_sd_CT - 1)).*(1/i_smax);  % Calculate i_sq by using formula from the lecture
        
        % interpolate CT-Line for intersection points
        interp_func = spline(i_sd_CT, i_sq_CT);

        % Evaluate the interpolated function of constant torque
        evaluation_points_constanttorque = linspace(-i_smax, 1, 100);  
        interpolated_values_constanttorque = ppval(interp_func, evaluation_points_constanttorque);

    %% MF-Ellipse (specific for operating point)
       	
        % Ellipse in parametric form
        omega_s_ellipse = omega_s_ellipse_set ;
        centerX = -i_sc;      % X-coordinate of the center
        centerY = 0;      % Y-coordinate of the center
        radiusX = ((u_smax/(omega_s_ellipse*L_sd)));                % Semi-major axis length by formula from the lecture
        radiusY = (u_smax/(omega_s_ellipse*L_sd)*1/(2*chi+1));      % Semi-minor axis length by formula from the lecture      
       
        theta = linspace(0, pi, 100);                               % Generate theta values from a range 0 to pi, 100 data-points
        
        % Parametric equation of the ellipse
        x = centerX + radiusX * cos(theta);
        y = centerY + radiusY * sin(theta);
      
    %% intersections
    % Using InterX Function ---> Not natively within MATLAB, but as ADD-on
   
        P_CT_mtpa = InterX([evaluation_points_constanttorque.*(1/i_smax);interpolated_values_constanttorque],[evaluation_points_mtpa;interpolated_values_mtpa]);  % Point A   CT-Line / MTPA-Line
        P_CT_MA = InterX([evaluation_points_constanttorque.*(1/i_smax);interpolated_values_constanttorque],[x_interpolated_ma;y_interpolated_ma]);                % Point B   CT-Line / MA-Circle
        P_mtpf_MA = InterX([x_interpolated_mtpf;y_interpolated_mtpf],[x_interpolated_ma;y_interpolated_ma]);                                                      % Point C MTPF-Line / MA-Circle
        

    %% Calculate Thresh-hold Values omega_s
        omega_spointa =  u_smax/(L_sd*sqrt((P_CT_mtpa(1,1)*i_smax+i_sc)^2+(P_CT_mtpa(2,1)*i_smax)^2*(2*chi+1)^2));       % Calculate omega_s(A)
        omega_spointb =  u_smax/(L_sd*sqrt((P_CT_MA(1,1)*i_smax+i_sc)^2+(P_CT_MA(2,1)*i_smax)^2*(2*chi+1)^2));           % Calculate omega_s(B)
        omega_spointc =  u_smax/(L_sd*sqrt((P_mtpf_MA(1,1)*i_smax+i_sc)^2+(P_mtpf_MA(2,1)*i_smax)^2*(2*chi+1)^2));       % Calculate omega_s(C)

   %% Output function
    %differentiate with if statements

% ------------Base - Speed Operation --> MTPA-Line------------

if omega_s_ellipse<=omega_spointa                      
      
        i_sq_output = P_CT_mtpa(2,1); %store the intersection point
        i_sd_output = P_CT_mtpa(1,1); %store the intersection point
end

%------------Base - Field Weakening --> Constant Torque Line------------

if omega_s_ellipse>omega_spointa && omega_s_ellipse<=omega_spointb  
        
        P_MF_CT = InterX([x/i_smax;y/i_smax],[evaluation_points_constanttorque.*(1/i_smax);interpolated_values_constanttorque]); % intersection between constant torque and MF line
        i_sq_output = P_MF_CT(2,1); %store the intersection point
        i_sd_output = P_MF_CT(1,1); %store the intersection point
end

%------------Field Weakening --> MA-Circle------------

if omega_s_ellipse>omega_spointb && omega_s_ellipse<=omega_spointc  %Field Weakening --> MA-Circle

       P_MF_MA = InterX([x/i_smax;y/i_smax],[i_sd_ma;i_sq_ma]);% intersection between constant torque and MF line
       i_sq_output = P_MF_MA(2,1); %store the intersection point
       i_sd_output = P_MF_MA(1,1); %store the intersection point
end

%------------Field Weakening --> MTPF Line------------

if omega_s_ellipse>omega_spointc 

       P_MF_CT = InterX([x/i_smax;y/i_smax],[x_interpolated_mtpf;y_interpolated_mtpf]);% intersection between constant torque and MF line
       i_sq_output = P_MF_CT(2,1); %store the intersection point
       i_sd_output = P_MF_CT(1,1); %store the intersection point

end
   



 %% Store Values of i_sd and i_sq into vector 

%Section 1
        upperlimit_mtpa_index = find(abs(interpolated_values_mtpa - P_CT_mtpa(2,1))<=0.0002);  % Getting the intersection Index Value to get a stop value for the MTPA Line, with tolerance of 0.0002, different approach?
        i_sq_currentpath =  interpolated_values_mtpa(1,upperlimit_mtpa_index:100);             % store the mtpa-line into the current-trajectory until the intersection point
        i_sd_currentpath =  evaluation_points_mtpa(1,upperlimit_mtpa_index:100);               % same for i_sd-component
        
        
%Section 2
        lower_limit = P_CT_mtpa(1,1)*i_smax;   % setting range, for how long operation on CT-line                                                               
        upper_limit = P_CT_MA(1,1).*i_smax;   
        evaluation_points_constanttorque_limit = linspace(lower_limit, upper_limit, 100);              % store the range    
        interpolated_values_currentpath = ppval(interp_func, evaluation_points_constanttorque_limit);  % store the values into current-trajectory section 2
  
%Section 3      
        lower_limit2 = P_CT_MA(1,1);           % setting range, for how long operation on MA-circle 
        upper_limit2 = P_mtpf_MA(1,1);
        section_MA_x = linspace(lower_limit2, upper_limit2, 100); % store the range 
        section_MA_y = ppval(splines_ma, section_MA_x);           % store the values into current-trajectory section 3
        
%Section 4
        index_intersection_mtpf_pointc = find(abs(y_interpolated_mtpf - P_mtpf_MA(2,1))<=0.0031);         % Getting the intersection Index Value to get a start value for the MTPF Line, with tolerance of 0.0031, different approach?
        index_intersection_mtpf_operatingpoint = find(abs(y_interpolated_mtpf - i_sq_output ) <= 0.008);  % Getting the intersection Index Value to get a stop value for the MTPF Line, with tolerance of 0.008, different approach?
        section_mtpf_x = x_interpolated_mtpf(index_intersection_mtpf_pointc(1,1):index_intersection_mtpf_operatingpoint(1,2)); % store the range 
        section_mtpf_y = y_interpolated_mtpf(index_intersection_mtpf_pointc(1,1):index_intersection_mtpf_operatingpoint(1,2)); % store the values into current-trajectory section 4

%% plot the result

close all;

sc_fc= i_smax; %Create a Scaling Factor, set to 1---> the plot is normalized on i_smax

figure;
plot(x_interpolated_ma.*sc_fc, y_interpolated_ma.*sc_fc,'red'); %MA-Circle
hold on
plot(evaluation_points_constanttorque.*(1/i_smax).*sc_fc ,interpolated_values_constanttorque.*sc_fc, Color='black'); %Constant Torque line
hold on
plot(evaluation_points_mtpa.*sc_fc,interpolated_values_mtpa.*sc_fc,'--', Color='red'); %MTPA
hold on
plot( x_interpolated_mtpf.*sc_fc , y_interpolated_mtpf.*sc_fc , '--',Color='green') % MTPF Line
hold on
plot( x/i_smax.*sc_fc,y/i_smax.*sc_fc,Color = 'green'); %MF-Ellipse 
hold on
plot( i_sd_currentpath.*sc_fc ,i_sq_currentpath.*sc_fc, 'LineWidth',2,Color = 'black'); % first current trajectory
hold on
plot(-i_sc/i_smax.*sc_fc,0 ,'o',Color='green'); %Short-Circuit Current
hold on
plot(P_CT_mtpa(1,1).*sc_fc, P_CT_mtpa(2,1).*sc_fc, 'o'); %Point A
hold on;
plot(P_CT_MA(1,1).*sc_fc,P_CT_MA(2,1).*sc_fc , 'o'); %Point B
hold on
plot(P_mtpf_MA(1).*sc_fc,P_mtpf_MA(2).*sc_fc , 'o'); %Point C
hold on
plot(i_sd_output.*sc_fc,i_sq_output.*sc_fc,'o'); %operating point
hold on
plot( evaluation_points_constanttorque_limit.*(1/i_smax).*sc_fc ,interpolated_values_currentpath.*sc_fc, 'LineWidth',2,Color = 'black'); %2nd current trajectory
hold on
plot( section_MA_x.*sc_fc  ,section_MA_y.*sc_fc, 'LineWidth',2,Color = 'black'); %third current trajectory
hold on
plot( section_mtpf_x.*sc_fc ,section_mtpf_y.*sc_fc, 'LineWidth',2,Color = 'black'); %4th current trajectory
axis equal;
xlabel('i_{sd}') % Add x-axis label
ylabel('i_{sq}') % Add y-axis label
legend("MA-Circle","Constant Torque line","MTPA-line","MTPF Line","MF-Ellipse","Current-Trajectory","Short-Circuit Current","Point A", "Point B","Point C","Operating Point", location ="best" );
grid on % Add gridlines
title("Task 1.2a)");


%% Task 1.2b)

T_e_plot = zeros(1,376);                  % Y-Axis, variable for constant Torque-values
omega_s_plot = linspace(0, 50000, 50000); % X-Axis, ranges from 0 to 50.000 rpm


omega_spointb_rpm =  omega_spointb *60/(2*pi*P); % until point A (omega values) store T_e=70Nm into T_e_plot
omega_spointc_rpm =  omega_spointc *60/(2*pi*P); % Umformung auf rpm


for i = 1:omega_spointb_rpm % Push constant T_e values until omega Point B is reached
T_e_plot(1,i) = T_e;
end

% i_sq and i_sd values differ based on operation mode
% create Torque-line for the dq-components generated in Field Weakening (MA-Circle)
omega_s_range = omega_spointb_rpm:1:omega_spointc_rpm;  

% Preallocate space for intersections
intersections = zeros(numel(omega_s_range), 2);
T_e_function_values = zeros(numel(omega_s_range), 1);

% Iterate over the specified range

for idx = 1:numel(omega_s_range)
    
    omega_s_ellipse = omega_s_range(idx)*2*pi/60*P;  %change from rpm to rad/s to obtain correct values from the formula

    % Generate ellipse for the intersection
    centerX = -i_sc;
    centerY = 0;
    radiusX = ((u_smax/(omega_s_ellipse*L_sd)));
    radiusY = (u_smax/(omega_s_ellipse*L_sd)*1/(2*chi+1));
    theta = linspace(0, pi, 100);
    x_ellipse = centerX + radiusX * cos(theta);
    y_ellipse = centerY + radiusY * sin(theta);
    
    % Generate MA-Circle for the intersection 
    theta = linspace(pi/2,pi,100);
    a_ma = 0;
    b_ma = 0;
    r_ma = 1;
    i_sd_ma = a_ma + r_ma.*cos(theta);
    i_sq_ma = b_ma + r_ma.*sin(theta);

    
    % Find intersection using InterX
    Point = InterX([x_ellipse/i_smax; y_ellipse/i_smax],[i_sd_ma; i_sq_ma]);
    
    if isempty(Point) %For Debugging
        disp('No intersection found for this value of omega_s') 
        continue;
    end
    
    % Store intersection (Taking the first intersection point, might not always be right)
    intersections(idx, :) = [Point(1,1), Point(2,1)];  % Point(1,1) is the x-coordinate and Point(2,1) is the y-coordinate of the first intersection

    % Calculate the according T_e value with the calculated dq-components
    % for each iteration
    T_e_function_values(idx) = 3/2*(Point(2,1)*i_smax*Psi_f+(L_sd-L_sq)*Point(2,1)*i_smax*Point(1,1)*i_smax);
end


% _________________________________________________________________________
%%%% The Following Code does the equivalent to the for-loop above, but in
%%%% the Field-Weakening Operation (MTPF-Line).
%__________________________________________________________________________

omega_s_ellipse_set_rpm = omega_s_ellipse_set *60/(2*pi*P);

% Initialize range of omega_s values
omega_s_range2 = omega_spointc_rpm:1:omega_s_ellipse_set_rpm;  

% Preallocate space for intersections
intersections2 = zeros(numel(omega_s_range2), 2);
T_e_function_values2 = zeros(numel(omega_s_range2), 1);

% Iterate over omega_s values
for idx = 1:numel(omega_s_range2)
    
    omega_s_ellipse2 = omega_s_range2(idx)*2*pi/60*P;

    % Generate ellipse
    centerX = -i_sc;
    centerY = 0;
    radiusX = ((u_smax/(omega_s_ellipse2*L_sd)));
    radiusY = (u_smax/(omega_s_ellipse2*L_sd)*1/(2*chi+1));
    theta = linspace(0, pi, 100);
    x_ellipse2 = centerX + radiusX * cos(theta);
    y_ellipse2 = centerY + radiusY * sin(theta);
    
    % MTPF Line

        speed_max1 = linspace(1, 1000000, 1000); % in rpm
        omega_s_mtpf = 2*pi.*speed_max1./60;
        i_o = u_smax./(omega_s_mtpf.*L_sd.*i_smax);
        i_sd_mtpf2 = -kappa + (1 + 2*chi).*kappa./(8*chi) - sqrt((i_o/2).^2 + ((1 + 2*chi)*kappa/(8*chi)).^2);
        i_sq_mtpf2 = 1./(2*chi+1).*sqrt(i_o.^2 - (kappa + i_sd_mtpf).^2);

        x_interpolated_mtpf2 = linspace(-2, -i_sc/i_smax, 1000);
        splines_mtpf = spline(i_sd_mtpf2, i_sq_mtpf2);
        y_interpolated_mtpf2 = ppval(splines_mtpf, x_interpolated_mtpf2);

        P2 = InterX([x_ellipse2/i_smax; y_ellipse2/i_smax],[x_interpolated_mtpf2;y_interpolated_mtpf2]);  % Find intersection using InterX

    
    if isempty(P2) % For Debugging
        disp('No intersection found for this value of omega_s')
        continue;
    end
    
    % Store intersection (Taking the first intersection point, might not always be right)
    intersections2(idx, :) = [P2(1,1), P2(2,1)];  % P(1,1) is the x-coordinate and P(2,1) is the y-coordinate of the first intersection

    % Calculate the according T_e value with the calculated dq-components
    % for each iteration
    T_e_function_values2(idx) = 3/2*(P2(2,1)*i_smax*Psi_f+(L_sd-L_sq)*P2(2,1)*i_smax*P2(1,1)*i_smax);
end

% ____________________________END OF FOLLOWING CODE________________________
%% plot 

% Create length of n, regarding the operation mode

x1=1:1:omega_spointb_rpm;
x2=omega_spointb_rpm:1:omega_spointc_rpm;
x3=omega_spointc_rpm:1:omega_s_ellipse_set*60/(2*pi*P); 

figure;
plot (x1,T_e_plot,color="blue"); % Constant T_e_plot
hold on
plot (x2, T_e_function_values,color="blue"); % Field-Weakening (MA-Circle)
hold on
plot (x3, T_e_function_values2,color="blue"); % Field-Weakening (MTPF-Line)
grid on
title("Task 2b)")
ylabel('Torque in Nm') % Add x-axis label
xlabel('Speed in rpm') % Add y-axis label

elseif kappa > 1
%__________________________________________________________________________
 %% Create Environment for kappa > 1
 %% MA-Circle normalized to i_smax
%__________________________________________________________________________
        theta = linspace(pi/2,pi,100);
        
        a_ma = 0;
        b_ma = 0;
        r_ma = 1;
        
        i_sd_ma = a_ma + r_ma.*cos(theta);
        i_sq_ma = b_ma + r_ma.*sin(theta); 
        
        % interpolate MA-Circle for intersection
        splines_ma = spline(i_sd_ma, i_sq_ma);
        x_interpolated_ma = linspace(min(i_sd_ma), max(i_sd_ma), 1000);
        y_interpolated_ma = ppval(splines_ma, x_interpolated_ma);
        
    %% MTPA Line
        
        % Generate values for i_s
        i_s = linspace(0, 1, 15);
        
        % Calculate i_sd and i_sq
        i_sd_mtpa = (kappa/(8*chi)) - sqrt((i_s.^2/2 + (kappa/(8*chi))^2));
        i_sq_mtpa = sqrt(i_s.^2 - i_sd_mtpa.^2);
        
        % interpolate MTPA-Line for intersection
        interp_func_mtpa = spline(i_sd_mtpa, i_sq_mtpa);
        evaluation_points_mtpa = linspace(-1, 0, 100); 
        interpolated_values_mtpa = ppval(interp_func_mtpa, evaluation_points_mtpa);
        
    %% MTPF Line
        speed_max1 = linspace(1, 1000000, 1000); % in rpm
        omega_s_mtpf = 2*pi.*speed_max1./60.*P; 

        i_o = u_smax./(omega_s_mtpf.*L_sd.*i_smax);
        i_sd_mtpf = -kappa + (1 + 2*chi).*kappa./(8*chi) - sqrt((i_o/2).^2 + ((1 + 2*chi)*kappa/(8*chi)).^2);
        i_sq_mtpf = 1./(2*chi+1).*sqrt(i_o.^2 - (kappa + i_sd_mtpf).^2);
        
        % interpolate MTPF-Line for intersection
        x_interpolated_mtpf = linspace(-2, -i_sc/i_smax, 1000);
        splines_mtpf = spline(i_sd_mtpf, i_sq_mtpf);
        y_interpolated_mtpf = ppval(splines_mtpf, x_interpolated_mtpf);


        
    %% Constant Torque line
        
        % Define the range of i_sd values for the left side
        i_sd_CT =  linspace(-i_smax, 0, 100);
        
        % Calculate i_sq using the given formula for the left side
        i_sq_CT = -2/3 * T_e/Psi_f * (1 ./ (2 * chi / i_sc .* i_sd_CT - 1)).*(1/i_smax);
        
        % Perform spline interpolation
        interp_func = spline(i_sd_CT, i_sq_CT);
        
        % interpolate CT-Line for intersection
        evaluation_points_constanttorque = linspace(-i_smax, 1, 100);  
        interpolated_values_constanttorque = ppval(interp_func, evaluation_points_constanttorque);

    %% MF-Ellipse (specific for operating point)
       	
        % Ellipse parameters
        omega_s_ellipse = omega_s_ellipse_set ;
        centerX = -i_sc;      % X-coordinate of the center
        centerY = 0;      % Y-coordinate of the center
        radiusX = ((u_smax/(omega_s_ellipse*L_sd)));                % Semi-major axis length
        radiusY = (u_smax/(omega_s_ellipse*L_sd)*1/(2*chi+1));      % Semi-minor axis length       
        
        % Generate theta values
        theta = linspace(0, pi, 100);
        
        % Parametric equations of the ellipse
        x = centerX + radiusX * cos(theta);
        y = centerY + radiusY * sin(theta);
      
    %% intersections (intersection point C not needed)
        P_CT_mtpa = InterX([evaluation_points_constanttorque.*(1/i_smax);interpolated_values_constanttorque],[evaluation_points_mtpa;interpolated_values_mtpa]);% Point A 
        P_CT_MA = InterX([evaluation_points_constanttorque.*(1/i_smax);interpolated_values_constanttorque],[x_interpolated_ma;y_interpolated_ma]);  % Point B

    %% Calculate Thresh-hold Values (omega (C) not needed)
        omega_spointa =  u_smax/(L_sd*sqrt((P_CT_mtpa(1,1)*i_smax+i_sc)^2+(P_CT_mtpa(2,1)*i_smax)^2*(2*chi+1)^2));       % Calculate omega_s(A)
        omega_spointb =  u_smax/(L_sd*sqrt((P_CT_MA(1,1)*i_smax+i_sc)^2+(P_CT_MA(2,1)*i_smax)^2*(2*chi+1)^2));           % Calculate omega_s(B)
       
   %% Output function - Does not contain Field Weakening (MTPF-Line)

%------------Base - Speed Operation --> MTPA-Line------------
if omega_s_ellipse<=omega_spointa 
      
        i_sq_output = P_CT_mtpa(2,1);
        i_sd_output = P_CT_mtpa(1,1);
end



%------------Base - Field Weakening --> Constant Torque Line------------
if omega_s_ellipse>omega_spointa && omega_s_ellipse<=omega_spointb  
        
        P_MF_CT = InterX([x/i_smax;y/i_smax],[evaluation_points_constanttorque.*(1/i_smax);interpolated_values_constanttorque]);% intersection between constant torque and MF line
        i_sq_output = P_MF_CT(2,1);
        i_sd_output = P_MF_CT(1,1);
end

 %------------Field Weakening --> MA-Circle------------
if omega_s_ellipse>omega_spointb  

       P_MF_MA = InterX([x/i_smax;y/i_smax],[i_sd_ma;i_sq_ma]);% intersection between constant torque and MF line

        if isempty(P_MF_MA)
            i_sq_output = 0;
            i_sd_output = -1;
        else
            i_sq_output = P_MF_MA(2,1);
            i_sd_output = P_MF_MA(1,1);
        end

end



 %% Store Values of i_sd and i_sq into vector (similar to above)

        upperlimit_mtpa_index = find(abs(interpolated_values_mtpa - P_CT_mtpa(2,1))<=0.01);    % Getting the intersection Index Value to get a stop variable for the MTPA Line
        i_sq_currentpath =  interpolated_values_mtpa(1,upperlimit_mtpa_index:100);
        i_sd_currentpath =  evaluation_points_mtpa(1,upperlimit_mtpa_index:100);
        
        
        lower_limit = P_CT_mtpa(1,1)*i_smax;
        upper_limit = P_CT_MA(1,1).*i_smax;   
        evaluation_points_constanttorque_limit = linspace(lower_limit, upper_limit, 100); 
        interpolated_values_currentpath = ppval(interp_func, evaluation_points_constanttorque_limit);
        
        lower_limit2 = P_CT_MA(1,1);
        upper_limit2 = i_sd_output;
        section_MA_x = linspace(lower_limit2, upper_limit2, 100);
        section_MA_y = ppval(splines_ma, section_MA_x);
        




%% plot the result

close all;

sc_fc= i_smax; %Create a Scaling Factor, set it to 1, to normalize the plot to i_smax
figure;

plot(x_interpolated_ma.*sc_fc, y_interpolated_ma.*sc_fc,'red'); %MA-Circle
hold on
plot(evaluation_points_constanttorque.*(1/i_smax).*sc_fc ,interpolated_values_constanttorque.*sc_fc, Color='black'); %Constant Torque line
hold on
plot(evaluation_points_mtpa.*sc_fc,interpolated_values_mtpa.*sc_fc,'--', Color='red'); %MTPA
hold on
plot( x_interpolated_mtpf.*sc_fc , y_interpolated_mtpf.*sc_fc , '--',Color='green') % MTPF Line
hold on
plot( x/i_smax.*sc_fc,y/i_smax.*sc_fc,Color = 'green'); %MF-Ellipse 
hold on
plot( i_sd_currentpath.*sc_fc ,i_sq_currentpath.*sc_fc, 'LineWidth',2,Color = 'black'); % first current trajectory
hold on
plot(-i_sc/i_smax.*sc_fc,0 ,'o',Color='green'); %Short-Circuit Current
hold on
plot(P_CT_mtpa(1,1).*sc_fc, P_CT_mtpa(2,1).*sc_fc, 'o'); %Point A
hold on;
plot(P_CT_MA(1,1).*sc_fc,P_CT_MA(2,1).*sc_fc , 'o'); %Point B
hold on
plot(i_sd_output.*sc_fc,i_sq_output.*sc_fc,'o'); %operating point
hold on
plot( evaluation_points_constanttorque_limit.*(1/i_smax).*sc_fc ,interpolated_values_currentpath.*sc_fc, 'LineWidth',2,Color = 'black'); %2nd current trajectory
hold on
plot( section_MA_x.*sc_fc  ,section_MA_y.*sc_fc, 'LineWidth',2,Color = 'black'); %third current trajectory
axis equal;
xlabel('i_{sd}') % Add x-axis label
ylabel('i_{sq}') % Add y-axis label
legend("MA-Circle","Constant Torque line","MTPA-line","MTPF Line","MF-Ellipse","Current-Trajectory","Short-Circuit Current","Point A", "Point B","Operating Point", location ="best" );
grid on % Add gridlines
title("Task 1.2a)");



%% Task 1.2b)______________________________________________________________ 
% Similar to case kappa<1, 
% but only one for-loop is needed, and the 
% case of a not intersecting MF-Ellipse with the MA-Circle is regarded
%__________________________________________________________________________

T_e_plot = zeros(1,376);                  % Y-Axis, preallocate space for T_e
omega_s_plot = linspace(0, 50000, 50000); % X-Axis, range from 0-50.000 rpm



omega_spointb_rpm =  omega_spointb *60/(2*pi*P);
omega_s_ellipse_set_rpm = omega_s_ellipse_set *60/(2*pi*P);


% until pointa (omega values) store 70Nm into T_e_plot
for i = 1:omega_spointb_rpm
T_e_plot(1,i) = T_e;
end


% Initialize range of omega_s values 
omega_s_range = omega_spointb_rpm:1:omega_s_ellipse_set_rpm;  

% Preallocate space for intersections
intersections = zeros(numel(omega_s_range), 2);
T_e_function_values = zeros(numel(omega_s_range), 1);


for idx = 1:numel(omega_s_range) % Iterate over omega_s values
    
    omega_s_ellipse = omega_s_range(idx)*2*pi/60*P;

    % Generate ellipse
    centerX = -i_sc;
    centerY = 0;
    radiusX = ((u_smax/(omega_s_ellipse*L_sd)));
    radiusY = (u_smax/(omega_s_ellipse*L_sd)*1/(2*chi+1));
    theta = linspace(0, pi, 100);
    x_ellipse = centerX + radiusX * cos(theta);
    y_ellipse = centerY + radiusY * sin(theta);
    
    % Generate MA-Circle
    theta = linspace(pi/2,pi,100);
    a_ma = 0;
    b_ma = 0;
    r_ma = 1;
    i_sd_ma = a_ma + r_ma.*cos(theta);
    i_sq_ma = b_ma + r_ma.*sin(theta);

    
    % Find intersection using InterX
    Point = InterX([x_ellipse/i_smax; y_ellipse/i_smax],[i_sd_ma; i_sq_ma]);

        if isempty(Point)    % if no intersection is found, the MF-Ellipse does not cross the MA-Circle anymore
            Point = [-1; 0]; % Set a Pre-defined Operation-Point within Limit in this case  
        else
            Point = InterX([x_ellipse/i_smax; y_ellipse/i_smax],[i_sd_ma; i_sq_ma]);
        end

    % Store intersection (Taking the first intersection point, might not always be right)
    intersections(idx, :) = [Point(1,1), Point(2,1)];  % Point(1,1) is the x-coordinate and Point(2,1) is the y-coordinate of the first intersection

    %Generate T_e-Values
    T_e_function_values(idx) = 3/2*(Point(2,1)*i_smax*Psi_f+(L_sd-L_sq)*Point(2,1)*i_smax*Point(1,1)*i_smax);
    disp("iteration done:"+idx); % For Debugging, how many iterations work
end


%% plot 

x1=1:1:omega_spointb_rpm;
x2=omega_spointb_rpm:1:omega_s_ellipse_set_rpm;

figure; % only two sections here, because no operation along MTPF-Line
plot (x1,T_e_plot,color="blue");
hold on
plot (x2, T_e_function_values,color="blue");
grid on
title("Task 2b)")
ylabel('Torque in Nm') % Add x-axis label
xlabel('Speed in rpm') % Add y-axis label



end 
%__________________________________________________________________________
%_______End of Kappa-Destinction and of the script for Task 2/3________%
%__________________________________________________________________________



%% __________________Function - section____________________________________ 



%_____________BEGIN dq-component generator-function________________________

function  [output] = isdq_components(omega_s, T_e)

%% initialize
clc;

Psi_f = 90e-3;        %Field flux linkage [Vs]
L_sd  = 200e-6;       %Stator inductance in d-axis [H]
L_sq  = 500e-6;       %Stator inductance in q-axis [H]
i_smax= 500;          %Maximum stator current (amplitude) [A]
U_dc  = 350;          %Dc-link voltage [V]
u_smax= U_dc/sqrt(3); %Maximum stator voltage
P     = 4;            %Pole pairs


i_sc   = 450;            %Short circuit current
kappa  = i_sc/i_smax;    %Short circuit current normalized on ismax
chi    = (L_sq-L_sd)/(2*L_sd);   %Saliency


%% Create Environment to Calculate the operation mode thresh-holds

        %% MA-Circle normalized to i_smax in parametric form

        theta = linspace(pi/2,pi,100);
        
        a_ma = 0;
        b_ma = 0;
        r_ma = 1;
        
        i_sd_ma = a_ma + r_ma.*cos(theta);
        i_sq_ma = b_ma + r_ma.*sin(theta); 
        
        % interpolate MA-Circle for intersection points
        splines_ma = spline(i_sd_ma, i_sq_ma);
        x_interpolated_ma = linspace(min(i_sd_ma), max(i_sd_ma), 1000);
        y_interpolated_ma = ppval(splines_ma, x_interpolated_ma);
        
        %% MTPA Line
        
        
        i_s = linspace(0, 1, 15);                                           % Generate values for i_s
       
        i_sd_mtpa = (kappa/(8*chi)) - sqrt((i_s.^2/2 + (kappa/(8*chi))^2)); % Calculate i_sd and i_sq with formulas from the lecture
        i_sq_mtpa = sqrt(i_s.^2 - i_sd_mtpa.^2);
        
        % interpolate MTPA-Line for intersection points
        interp_func_mtpa = spline(i_sd_mtpa, i_sq_mtpa);
        evaluation_points_mtpa = linspace(-1, 0, 100);  
        interpolated_values_mtpa = ppval(interp_func_mtpa, evaluation_points_mtpa);
        
        %% MTPF Line
        speed_max = linspace(1, 1000000, 1000); % Calculate i_sd and i_sq with formulas from the lecture
        omega_s_mtpf = 2*pi.*speed_max./60*P;   % in rps

        % Formulas from the lecture:
        i_o = u_smax./(omega_s_mtpf.*L_sd.*i_smax);
        i_sd_mtpf = -kappa + (1 + 2*chi).*kappa./(8*chi) - sqrt((i_o/2).^2 + ((1 + 2*chi)*kappa/(8*chi)).^2);
        i_sq_mtpf = 1./(2*chi+1).*sqrt(i_o.^2 - (kappa + i_sd_mtpf).^2);
        
        % interpolate MTPF-Line for intersection points
        x_interpolated_mtpf = linspace(-2, -i_sc/i_smax, 1000);
        splines_mtpf = spline(i_sd_mtpf, i_sq_mtpf);
        y_interpolated_mtpf = ppval(splines_mtpf, x_interpolated_mtpf);

        %% MF-Ellipse
    	
         % Ellipse in parametric form
        centerX = -i_sc;  % X-coordinate of the center
        centerY = 0;      % Y-coordinate of the center
        radiusX = ((u_smax/(omega_s*L_sd)));                % Semi-major axis length by formula from the lecture
        radiusY = (u_smax/(omega_s*L_sd)*1/(2*chi+1));      % Semi-minor axis length by formula from the lecture        
        
        % Generate theta values
        theta = linspace(0, pi, 100);
        
        % Parametric equation of the ellipse
        x = centerX + radiusX * cos(theta);
        y = centerY + radiusY * sin(theta);
        
        %% Constant Torque line
        
       % Define the range of i_sd values for the 2nd quadrant
        i_sd_CT =  linspace(-i_smax, 0, 100);
        
        % Calculate i_sq by using formula from the lecture
        i_sq_CT = -2/3 * T_e/Psi_f * (1 ./ (2 * chi / i_sc .* i_sd_CT - 1)).*(1/i_smax);
        
        % spline-interpolate CT-Line for intersection points
        interp_func = spline(i_sd_CT, i_sq_CT);
        
        % Evaluate the interpolated function of constant torque
        evaluation_points_constanttorque = linspace(-i_smax, 1, 100);  
        interpolated_values_constanttorque = ppval(interp_func, evaluation_points_constanttorque);
        
        %% intersections
         % Using InterX Function ---> Not natively within MATLAB, but as ADD-on
        P_CT_mtpa = InterX([evaluation_points_constanttorque.*(1/i_smax);interpolated_values_constanttorque],[evaluation_points_mtpa;interpolated_values_mtpa]);  % Point A   CT-Line / MTPA-Line
        P_mtpf_MA = InterX([x_interpolated_mtpf;y_interpolated_mtpf],[x_interpolated_ma;y_interpolated_ma]);                                                      % Point B   CT-Line / MA-Circle
        P_CT_MA = InterX([evaluation_points_constanttorque.*(1/i_smax);interpolated_values_constanttorque],[x_interpolated_ma;y_interpolated_ma]);                % Point C MTPF-Line / MA-Circle
        
        %% Calculate Thresh-hold Values
        omega_spointa =  u_smax/(L_sd*sqrt((P_CT_mtpa(1,1)*i_smax+i_sc)^2+(P_CT_mtpa(2,1)*i_smax)^2*(2*chi+1)^2));       % Calculate omega_s(A)
        omega_spointb =  u_smax/(L_sd*sqrt((P_CT_MA(1,1)*i_smax+i_sc)^2+(P_CT_MA(2,1)*i_smax)^2*(2*chi+1)^2));           % Calculate omega_s(B)
        omega_spointc =  u_smax/(L_sd*sqrt((P_mtpf_MA(1,1)*i_smax+i_sc)^2+(P_mtpf_MA(2,1)*i_smax)^2*(2*chi+1)^2));       % Calculate omega_s(C)

%% Output function
 %differentiate with if statements
% _________________________________________________________________________
% ------------Base - Speed Operation --> MTPA-Line-------------------------
if omega_s<=omega_spointa 
      
        i_sq_output = P_CT_mtpa(2,1); %store the intersection point
        i_sd_output = P_CT_mtpa(1,1); %store the intersection point
        disp("Base-Speed Operation")
end

%------------Base - Field Weakening --> Constant Torque Line---------------
if omega_s>omega_spointa && omega_s<=omega_spointb 
        
        P_MF_CT = InterX([x/i_smax;y/i_smax],[evaluation_points_constanttorque.*(1/i_smax);interpolated_values_constanttorque]); % intersection between constant torque and MF line
        i_sq_output = P_MF_CT(2,1); %store the intersection point
        i_sd_output = P_MF_CT(1,1); %store the intersection point
        disp("Base Field Weakening")
end

%------------Field Weakening --> MA-Circle---------------------------------
if omega_s>omega_spointb && omega_s<=omega_spointc  

       P_MF_MA = InterX([x/i_smax;y/i_smax],[i_sd_ma;i_sq_ma]); % intersection between MA and MF line
       i_sq_output = P_MF_MA(2,1); %store the intersection point
       i_sd_output = P_MF_MA(1,1); %store the intersection point
       disp("Field Weakening MA-Circle")
end

%------------Field Weakening --> MTPF Line---------------------------------
if omega_s>omega_spointc 

    P_MF_CT = InterX([x/i_smax;y/i_smax],[x_interpolated_mtpf;y_interpolated_mtpf]); % intersection between MTPF and MF line
    i_sq_output = P_MF_CT(2,1); %store the intersection point
    i_sd_output = P_MF_CT(1,1); %store the intersection point
    disp("Field Weakening MTPF Line")
end
 % ________________________________________________________________________

% Define the Output of the function:
output = [i_sd_output, i_sq_output];


% Display output Text
output_text1 = ("the corresping i_sd value is: ");
output_text2 = ("the corresping i_sq value is: ");
disp(output_text1+i_sd_output)
disp(output_text2+i_sq_output)



end 
%_____________END dq-component generator-function__________________________


%_____________InterX_Function (not self-written)___________________________
function intersection = InterX(L1,varargin)
%INTERX Intersection of curves
%   P = INTERX(L1,L2) returns the intersection points of two curves L1 
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply 
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   P = INTERX(L1) returns the self-intersection points of L1. To keep
%   the code simple, the points at which the curve is tangent to itself are
%   not included. P = INTERX(L1,L1) returns all the points of the curve 
%   together with any self-intersection points.
%   
%   Example:
%       t = linspace(0,2*pi);
%       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
%       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
%       P = InterX([x1;y1],[x2;y2]);
%       plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro')
%   Author : NS
%   Version: 3.0, 21 Sept. 2010
%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point. 
%   Each factor of the 'C' arrays is essentially a matrix containing 
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.
    %...Argument checks and assignment of L2
    error(nargchk(1,2,nargin));
    if nargin == 1,
        L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
    else
        L2 = varargin{1}; hF = @le;
    end
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';
    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i),P = zeros(2,0);return; end;
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    intersection = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
              
    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end
%_______________________END InterX_________________________________________
%% ______________________END Function-Section______________________________


