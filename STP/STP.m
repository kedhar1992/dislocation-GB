% Copyright 2013 Max-Planck-Institut für Eisenforschung GmbH
%% Example of function used to calculate mprime parameter
tabularasa;
% author: d.mercier@mpie.de

%%%%%%%%%%%%%% INPUT  - START
%%%%  112 TILT AXIS
% euler_1 = [90 35.26 225]; % Bunge notation (in degrees)  % CENTRAL GRAIN    225 (from stabix) instead of 135 (from Kedharnath)
% tilt_angle = [90.0000  101.5370   74.5013   66.9261   55.9523   50.7685   45.5847   32.8783   27.0357   11.5370    1.1694];  % in Euler notations
% tilt = [0 11.5 15.5 23.1  34.1  39.2   44.4  57  63   78  88 ];
% GBE = [0 977.06596 1204 1265.15862 1328 1401.74558 1286.81087 1244 1009.12984 1409.62395 1615.62124]; 

%%%%  110 TILT AXIS
% euler_1 = [35.2641 90 315]; % Bunge notation (in degrees)  % CENTRAL GRAIN    315  (from stabix) instead of 45 (from Kedharnath)
% tilt_angle = [35.2644   27.9384   23.8426  195.2144  188.7390  254.2068  164.7856   94.2568  324.7356  118.2103  128.6367];  % in Euler notations
% GBE = [0 859  1003   1302 1504  1423  1255  1356  1454   1203  1410];
% tilt = [0  7.3   11.4   20   27  39   50    58   70  83   86];

%%%%  111 TWIST AXIS
euler_1 = [180 54.7403 45]; % Bunge notation (in degrees)  % CENTRAL GRAIN
tilt_angle = [180.0000  189.4300  163.5736  158.2132  152.2042  141.7868  223.5736  129.4300  120.0000  110.5700  106.8264   98.2132   92.2042];  % in Euler notations
tilt = [0  9.4  16.4   21.8   27.8   38.2    43.6   50.6    60   69   73   81    87  ];
GBE = [0  880   1186  1276  1348   1333    1574   1598   1451   1599   1572    1333   1408];

% Index 11 = -110 slip plane 111 slip direction
% Index 24 = 11-2 slip plabe 111 slip direction
ind_slip_in = 24;                   % Change slip systems by using slip_systems.m
%%%%%%%%%%%%%% INPUT  - END

GBN = [1 1 1]';
mp_store = [0     0     0 0 0   0 0 0];
M1 = [0 0 0]';
M2 = [0 0 0]';

for i = 1:length(tilt_angle)
    euler_2 = [tilt_angle(i) euler_1(2) euler_1(3)]; % Bunge notation (in degrees)  % ROTATED GRAIN
    
    mp_all = [0     0 0 0   0 0 0];
    for j = 1:48
        %% Data of user
        material = 'Ta';
        phase = 'bcc';
        all_slips = slip_systems(phase);
        ind_slip_out = j;                    % Change slip systems
        incoming_slip = all_slips(:,:,ind_slip_in);   % first is plane, second is direction
        outgoing_slip = all_slips(:,:,ind_slip_out);
        
        %% Calculations
        % Eulers angle to rotation matrix
        rot_mat1 = eulers2g(euler_1);
        rot_mat2 = eulers2g(euler_2);
        
        % Lattice parameters
        ca_ratio = listLattParam(material, phase);
        
        % if strcmp(phase, 'bcc') == 1
        %     Miller-Bravais to Cartesian coordinates
        %     n1_cart = millerbravaisplane2cart(incoming_slip(1,:), ca_ratio(1))';
        %     d1_cart = millerbravaisdir2cart(incoming_slip(2,:), ca_ratio(1))';
        %     n2_cart = millerbravaisplane2cart(outgoing_slip(1,:), ca_ratio(1))';
        %     d2_cart = millerbravaisdir2cart(outgoing_slip(2,:), ca_ratio(1))';
        % end
        
        GBN_cart = GBN;
        n1_cart = incoming_slip(1,:);
        d1_cart = incoming_slip(2,:);
        n2_cart = outgoing_slip(1,:);
        d2_cart = outgoing_slip(2,:);
        
        % Apply rotation matrix to slip systems
        GBN_cart_rot = GBN_cart;
        n1_cart_rot = rot_mat1.'*n1_cart';
        d1_cart_rot = rot_mat1.'*d1_cart';
        n2_cart_rot = rot_mat2.'*n2_cart';
        d2_cart_rot = rot_mat2.'*d2_cart';
        
        % Normalization of input vectors
        n1 = n1_cart_rot ./ norm(n1_cart_rot);
        n2 = n2_cart_rot ./ norm(n2_cart_rot);
        d1 = d1_cart_rot ./ norm(d1_cart_rot);   % b_in
        b_in = d1';
        d2 = d2_cart_rot ./ norm(d2_cart_rot);   % b_out
        b_out = d2';
        
        % Check of orthogonality
        test_vectors_orthogonality(n1, d1);
        test_vectors_orthogonality(n2, d2);
        
        % Get cosine from the dot product
        cosine_n = cosFromVectors(n1, n2);
        cosine_d = cosFromVectors(d1, d2);
        cosine_n_inv = cosFromVectors(n2, n1);
        cosine_d_inv = cosFromVectors(d2, d1);
        
        % Get angle from cosine
        phi = ang_from_vectors(n1, n2);
        kappa = ang_from_vectors(d1, d2);
        
        % m' calculation
        mp = abs(cosine_n * cosine_d);
        mp_inv = abs(cosine_n_inv * cosine_d_inv);
        mp_function = mprime(n1, d1, n2, d2);
        mp_all(j,:) = [mp_function    outgoing_slip(1,:)    outgoing_slip(2,:)];
        
        % RBV calculation
        R_in  = eulers2g(euler_1);
        R_out = eulers2g(euler_2);
        rotated_b_in = R_in'*b_in';
        rotated_b_out = R_out'*b_out';
        rbv_function = residual_Burgers_vector(rotated_b_in, rotated_b_out);
        rbv_all(j,:) = [rbv_function    outgoing_slip(1,:)    outgoing_slip(2,:)];
        
        % N factor calculation
        n_function = N_factor(n1, d1, n2, d2);
        n_all(j,:) = [n_function    outgoing_slip(1,:)    outgoing_slip(2,:)];
        
        % LRB calculation
        GB_trace_angle = 90;
        if length(tilt_angle) > 12
            GB_inclination = 0;   % Twist conf.
        else
            GB_inclination = 90;   % Tilt conf.
        end
        
        rotataxis = rotation('axis',zvector,'angle',GB_trace_angle*degree);
        gbVec = rotate(xvector, rotataxis);
        gbVec = gbVec./norm(gbVec);
        perp_gb = cross(zvector, gbVec);
        rotataxis = rotation('axis',gbVec,'angle',GB_inclination*degree);
        gbVecInc = rotate(perp_gb, rotataxis) + gbVec;
        
        dGB = gbVec - gbVecInc;
        dGB = dGB./norm(dGB);
        GB_angle = dGB.xyz;
        
        l1 = cross(GB_angle, d1');
        l2 = cross(GB_angle, d2');
        l1 = l1./norm(l1);
        l2 = l2./norm(l2);
        
        LRB_function = LRB_parameter(l1, d1, l2, d2);
        LRB_all(j,:) = [LRB_function    outgoing_slip(1,:)    outgoing_slip(2,:)];
        
        % Lamba calculation
        lamda_function = lambda(n1, d1, n2, d2, 15, 45);
        lamda_all(j,:) = [lamda_function    outgoing_slip(1,:)    outgoing_slip(2,:)];
        
        % Resolved shear stress calculation
        % Index 11 = -110 slip plane 111 slip direction
        % Index 24 = 11-2 slip plabe 111 slip direction
        if length(tilt_angle) < 12    % Tilt conf.
            if ind_slip_in == 11
                sigma_a = [0 1 0; 1 0 0; 0 0 0]; % Shear XY
            elseif ind_slip_in == 24
                sigma_a = [0 0 1; 0 0 0; 1 0 0]; % Shear XZ
            else
                disp("Applied stress tensor is not appropriate on slip system considered")
            end
        else              % Twist conf.
            if ind_slip_in == 11
                sigma_a = [0 0 1; 0 0 0; 1 0 0]; % Shear XZ
            elseif ind_slip_in == 24
                sigma_a = [0 0 0; 0 0 1; 0 1 0]; % Shear YZ
            else
                disp("Applied stress tensor is not appropriate on slip system considered")
            end
        end
            
        tau_function  = resolved_shear_stress(euler_2, n2_cart, d2_cart, sigma_a);
        tau_all(j,:) = [tau_function    outgoing_slip(1,:)    outgoing_slip(2,:)];
        
        % S factor calculation
        s_function = s_factor(n1, d1, l1, n2, d2, l2);
        s_all(j,:) = [s_function    outgoing_slip(1,:)    outgoing_slip(2,:)];
        
        
    end
    % m' calculation - MAXIMIZE
    mp_sort = sortrows(mp_all,1,'descend');
    mp_max = mp_sort(1,:);
    mp_store(i,:) = [tilt_angle(i)  mp_max];
    
    % RBV calculation - MINIMIZE
    rbv_sort = sortrows(rbv_all,1,'ascend');
    rbv_min = rbv_sort(1,:);
    rbv_store(i,:) = [tilt_angle(i)  rbv_min];
    
    % N factor calculation - MAXIMIZE
    n_sort = sortrows(n_all,1,'descend');
    n_max = n_sort(1,:);
    n_store(i,:) = [tilt_angle(i)  n_max];
    
    % LRB calculation - MAXIMIZE
    LRB_sort = sortrows(LRB_all,1,'descend');
    LRB_max = LRB_sort(1,:);
    LRB_store(i,:) = [tilt_angle(i)  LRB_max];
    
    % Lamba calculation  -  critical value dependent
    lamda_sort = sortrows(lamda_all,1,'descend');
    lamda_max = lamda_sort(1,:);
    lamda_store(i,:) = [tilt_angle(i)  lamda_max];
    
    % Resolved shear stress calculation - MAXIMIZE
    tau_sort = sortrows(tau_all,1,'descend');
    tau_max = tau_sort(1,:);
    tau_store(i,:) = [tilt_angle(i)  tau_max];    
    
    % S factor calculation - MAXIMIZE
    s_sort = sortrows(s_all,1,'descend');
    s_max = s_sort(1,:);
    s_store(i,:) = [tilt_angle(i)  s_max];    
    
    
    geo_store(i,:) = [tilt(i)  tilt_angle(i)  GBE(i)  mp_max     rbv_min    n_max     LRB_max    lamda_max   s_max ];
    stress_store(i,:) = [tilt(i)  tilt_angle(i)   GBE(i)   tau_max];
    origin_tilt(i,:) = [tilt(i)   mp_max(1)     rbv_min(1)    n_max(1)     LRB_max(1)    lamda_max(1)   s_max(1)   tau_max(1) ];
    origin_gbe_unsort(i,:) = [GBE(i)    mp_max(1)     rbv_min(1)    n_max(1)     LRB_max(1)    lamda_max(1)   s_max(1)   tau_max(1) ];
    
end

% Variation of values of parameters as a function of GBE 
geo_gbe = sortrows(geo_store, 3);
stress_gbe = sortrows(stress_store, 3);
origin_gbe = sortrows(origin_gbe_unsort,1)
origin_tilt

% GBE plots
figure
plot(geo_gbe(:,3), mp_store(:,2),'*-')
title('m')
figure
plot(geo_gbe(:,3), rbv_store(:,2),'*-')
title('RBV (minimized)')
figure
plot(geo_gbe(:,3), n_store(:,2),'*-')
title('N-factor')
figure
plot(geo_gbe(:,3), LRB_store(:,2),'*-')
title('LRB')
figure
plot(geo_gbe(:,3), lamda_store(:,2),'*-')
title('Lambda')
figure
plot(geo_gbe(:,3), tau_store(:,2),'*-')
title('Resolved shear stress')
figure
plot(geo_gbe(:,3), s_store(:,2),'*-')
title('s-factor')


% Tilt plots
figure
plot(tilt, mp_store(:,2),'*-')
title('m')
figure
plot(tilt, rbv_store(:,2),'*-')
title('RBV (minimized)')
figure
plot(tilt, n_store(:,2),'*-')
title('N-factor')
figure
plot(tilt, LRB_store(:,2),'*-')
title('LRB')
figure
plot(tilt, lamda_store(:,2),'*-')
title('Lambda')
figure
plot(tilt, tau_store(:,2),'*-')
title('Resolved shear stress')
figure
plot(tilt, s_store(:,2),'*-')
title('s-factor')


%% Display
% commandwindow;
% figure;
% vis_lattice(phase, euler_1, ind_slip_in);
% figure;
% vis_lattice(phase, euler_2, ind_slip_out);
% disp(euler_1');
% disp(euler_2');
% disp(incoming_slip);
% disp(outgoing_slip);
% disp(phi);
% disp(kappa);
% disp(mp);
% disp(mp_inv);

