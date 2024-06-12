function [x_est, b_est, Ax_f_instance, Ay_f_instance, Az_f_instance, p_f_instance, q_f_instance, r_f_instance, AoA_f_instance] = SID230255137(c_k, d_k, t, dt)

%% Note: Rename the function in the format of SID + your student ID as on Blackboard (e.g. if your ID is 21010000, name the function as SID21010000 and submit it as SID21010000.m)
%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   c_k: input measurements, a N x 6 matrix with rows the recordings of different time instances and the columns being [A_x_m, A_y_m, A_z_m, p_m, q_m, r_m]
%   d_k: output measurements, a N x 12 matrix with rows the recordings of different time instances and the columns being [x_GPS_m, y_GPS_m, z_GPS_m, u_GPS_m, v_GPS_m, w_GPS_m, phi_GPS_m, theta_GPS_m, psi_GPS_m, V_TAS_m, alpha_m, beta_m]
%   t: time vector
%   dt: uniform time step size
%%%%%%%%%%%%%%%%%%%%%%%
% Output Variables
%   x_est: estimated state trajectory, a N x 12 matrix with rows the recordings of different time instances and the columns being [x_E, y_E, z_E, u, v, w, \phi, \theta, \psi, V_{wxE}, V_{wyE}, V_{wzE}]
%   b_est: estimated state trajectory, a N x 6 matrix with rows the recordings of different time instances and the columns being [b_A_x, b_A_y, b_A_z, b_p, b_q, b_r]
%   Ax_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of A_x. In case of no faults detected, please return a value of 0,
%   Ay_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of A_y. In case of no faults detected, please return a value of 0,
%   Az_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of A_z. In case of no faults detected, please return a value of 0,
%   p_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of p. In case of no faults detected, please return a value of 0,
%   q_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of q. In case of no faults detected, please return a value of 0,
%   r_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of r. In case of no faults detected, please return a value of 0,
%   AoA_f_instance: a number equal to the corresponding time instances t for the first detection of faults in measurements of angle of attack sensor. In case of no faults detected, please return a value of 0. 
%%%%%%%%%%%%%%%%%%%%%%%

%% Your code here
    % Initialization
    z_k = d_k;
    u_k = c_k;
    Ts=dt;     %time step (already provided by the data)
    N=length(t);      %total number of steps
    
    % omega_A_x omega_A_y omega_A_z omega_p omega_q omega_r
    stdw = [0.01, 0.01, 0.01, 0.01*pi()/180, 0.01*pi()/180, 0.01*pi()/180, 5, 5, 5, 1.5*pi()/180, 1.5*pi()/180, 1.5*pi()/180];  %standard deviation of w
    
    % nu_x nu_y nu_z nu_u nu_v nu_w nu_phi nu_theta nu_psi nu_V nu_alpha nu_beta
    stdv = [5, 5, 10, 0.1, 0.1, 0.1, 0.1*pi()/180, 0.1*pi()/180, 0.1*pi()/180, 0.1, 0.1*pi()/180, 0.1*pi()/180];      %standard deviation of v
    
    % x y z u v w phi theta psi v_wxE v_wyE v_wzE
    Ex_0 = [d_k(1,1) d_k(1,2) d_k(1,3) d_k(1,10) 0 0 d_k(1,7) d_k(1,8) d_k(1,9) 0 0 0 0 0 0 0 0 0];   %expected value of x_0
    stdx_0 = [1 1 1 1 100 100 0.1*pi()/180 0.1*pi()/180 0.1*pi()/180 100 100 100 200 200 200 200 200 200];  %standard deviation of x_0

    xhat_km1_km1 = Ex_0; % x(0|0) = E{x_0}
    P_km1_km1 = diag(stdx_0.^2);  % P(0|0) = P(0)
    Q=diag(stdw.^2);
    R=diag(stdv.^2);
    
    n = length(xhat_km1_km1); % n: state dimension
    m = size(u_k, 2);     % m: observation dimension
    p = size(z_k, 2);     % m: observation dimension
    u_km1 = [zeros(1,m); u_k]; % shifted to have the right indices
    
    % Preallocate storage
    stdx_cor  = zeros(N, n);  % \sigma(k-1|k-1), standard deviation of state estimation error (hint: diagonal elements of P(k-1|k-1))
    x_cor     = zeros(N, n);  % \hat{x}(k-1|k-1), previous estimation
    K       = cell(N, 1);   % K(k) Kalman Gain
    innov     = zeros(N, p);  % y(k)-y(k|k-1), innovation, with y(k|k-1)=h(\hat{x}(k|k-1),u(k|k-1),k);
    
    for k=1:N
        % Step 1: Prediction
        [t_nonlin, x_nonlin] = ode45(@(t,x) funcf(x, u_km1(k,:), t), [0 Ts], xhat_km1_km1);
        xhat_k_km1 = x_nonlin(end,:); % x(k|k-1) (prediction)
    
        % Step 2: Covariance matrix of state prediction error / Minimum
        % prediction MSE
        [Phi_km1, Gamma_km1] = funcLinDisDyn(xhat_km1_km1, u_km1(k,:), Ts); % Phi(k,k-1), Gamma(k,k-1)
        P_k_km1 = Phi_km1 * P_km1_km1 * Phi_km1' + Gamma_km1 * Q * Gamma_km1'; % P(k|k-1) (prediction)
    
        % Step 3: Kalman Gain
        H_k = funcLinDisObs(xhat_k_km1, u_km1(k,:), []);
        Ve = (H_k * P_k_km1 * H_k' + R); % Pz(k|k-1) (prediction)
        K_k = P_k_km1 * H_k' / Ve; % K(k) (gain)
        
        % Step 4: Measurement Update (Correction)
        z_k_km1 = funch(xhat_k_km1,u_km1(k,:),[]); % z(k|k-1) (prediction of output)
        xhat_k_k = xhat_k_km1 + (z_k(k,:) - z_k_km1)*K_k'; % x(k|k) (correction)
    
        % Step 5: Correction for Covariance matrix of state Estimate error /
        % Minimum MSE
        I_KH = eye(n) - K_k * H_k;
        P_k_k = I_KH * P_k_km1 * I_KH' + K_k * R * K_k'; % P(k|k) (correction)
    
        % Save data: State estimate and std dev
        stdx_cor(k,:) = sqrt(diag(P_km1_km1)); % \sigma(k-1|k-1) Standard deviation of state estimation error
        x_cor(k,:) = xhat_km1_km1; % \hat{x}(k-1|k-1), estimated state
        K{k,1} = K_k; % K(k) (gain)
        innov(k,:)= z_k(k,:) - z_k_km1;
        
        % Recursive step
        xhat_km1_km1 = xhat_k_k; 
        P_km1_km1 = P_k_k;
    end
    
    %% Assign the variables
    
    x_est = x_cor(1:N,1:12);
    b_est = x_cor(1:N,13:18);
    
    %% dynamics in continous time
    function x_dot = funcf(xin, uin, t)
        % Extract values
        x = xin(1);
        y = xin(2);
        z = xin(3);
        u = xin(4);
        v = xin(5);
        w = xin(6);
        phi = xin(7);
        theta = xin(8);
        psi = xin(9);
        v_wxE = xin(10);
        v_wyE = xin(11);
        v_wzE = xin(12);
    
        b_A_x = xin(13);
        b_A_y = xin(14);
        b_A_z = xin(15);
        b_p = xin(16);
        b_q = xin(17);
        b_r = xin(18);
    
        A_x = uin(1);
        A_y = uin(2);
        A_z = uin(3);
        p = uin(4);
        q = uin(5);
        r = uin(6);
        
        g = 9.81;
        Ts = 0.01;
    
        omega_b_A_x = 0;
        omega_b_A_y = 0;
        omega_b_A_z = 0;
        omega_b_p = 0;
        omega_b_q = 0;
        omega_b_r = 0;
    
        % original continuous-time nonlinear system equations
        x_dot =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi) - (v*cos(phi)-w*sin(phi))*sin(psi)+v_wxE;
        y_dot =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi) + (v*cos(phi)-w*sin(phi))*cos(psi)+v_wyE;
        z_dot =-u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+v_wzE;
        
        u_dot =(A_x-b_A_x)-g*sin(theta)+(r-b_r)*v-(q-b_q)*w;
        v_dot =(A_y-b_A_y)+g*cos(theta)*sin(phi)+(p-b_p)*w-(r-b_r)*u;
        w_dot =(A_z-b_A_z)+g*cos(theta)*cos(phi)+(q-b_q)*u-(p-b_p)*v;
        
        phi_dot =(p-b_p)+(q-b_q)*sin(phi)*tan(theta)+(r-b_r)*cos(phi)*tan(theta);
        theta_dot =(q-b_q)*cos(phi)-(r-b_r)*sin(phi);
        psi_dot =(q-b_q)*sin(phi)/cos(theta)+(r-b_r)*cos(phi)/cos(theta);
         
        v_wxE_dot =0;
        v_wyE_dot =0;
        v_wzE_dot =0;  
        
        b_A_x_dot = omega_b_A_x/Ts;
        b_A_y_dot = omega_b_A_y/Ts;
        b_A_z_dot = omega_b_A_z/Ts;
        b_p_dot = omega_b_p/Ts;
        b_q_dot = omega_b_q/Ts;
        b_r_dot = omega_b_r/Ts;
    
        % original, continuous-time, noiseless nonlinear dynamics
        x_dot = [x_dot; y_dot; z_dot; u_dot; v_dot; w_dot; phi_dot; theta_dot; psi_dot; v_wxE_dot; v_wyE_dot; v_wzE_dot; b_A_x_dot; b_A_y_dot; b_A_z_dot; b_p_dot; b_q_dot; b_r_dot];
    end
    
    function d = funch(xin,uin,t)
        % Extract values
        x = xin(1);
        y = xin(2);
        z = xin(3);
        u = xin(4);
        v = xin(5);
        w = xin(6);
        phi = xin(7);
        theta = xin(8);
        psi = xin(9);
        v_wxE = xin(10);
        v_wyE = xin(11);
        v_wzE = xin(12);
    
        b_A_x = xin(13);
        b_A_y = xin(14);
        b_A_z = xin(15);
        b_p = xin(16);
        b_q = xin(17);
        b_r = xin(18);
    
        A_x = uin(1);
        A_y = uin(2);
        A_z = uin(3);
        p = uin(4);
        q = uin(5);
        r = uin(6);
    
        g = 9.81;
    
        % original continuous-time nonlinear system equations
        x_GPS = x;
        y_GPS = y;
        z_GPS = z;
            
        u_GPS = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi)-(v*cos(phi)-w*sin(phi))*sin(psi)+v_wxE;
        v_GPS = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi)+(v*cos(phi)-w*sin(phi))*cos(psi)+v_wyE;
        w_GPS = -u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+v_wzE;
            
        phi_GPS = phi;
        theta_GPS = theta;
        psi_GPS = psi;
        
        V = sqrt(u^2+v^2+w^2);
        alpha = atan(w/u);
        beta = atan(v/sqrt(u^2+w^2));
        
        % Final observation dynamics
        d = [x_GPS, y_GPS, z_GPS, u_GPS, v_GPS, w_GPS, phi_GPS, theta_GPS, psi_GPS, V, alpha, beta];
    end
    
    %% lineaized dynamics in discrete time
    function [Phi,Gamma] = funcLinDisDyn(x_linpt,u_linpt,Ts)
        % Extract values
        x = x_linpt(1);
        y = x_linpt(2);
        z = x_linpt(3);
        u = x_linpt(4);
        v = x_linpt(5);
        w = x_linpt(6);
        phi = x_linpt(7);
        theta = x_linpt(8);
        psi = x_linpt(9);
        v_wxE = x_linpt(10);
        v_wyE = x_linpt(11);
        v_wzE = x_linpt(12);
    
        b_A_x = x_linpt(13);
        b_A_y = x_linpt(14);
        b_A_z = x_linpt(15);
        b_p = x_linpt(16);
        b_q = x_linpt(17);
        b_r = x_linpt(18);
    
        A_x_m = u_linpt(1);
        A_y_m = u_linpt(2);
        A_z_m = u_linpt(3);
        p_m = u_linpt(4);
        q_m = u_linpt(5);
        r_m = u_linpt(6);
    
        omega_A_x = 0;
        omega_A_y = 0;
        omega_A_z = 0;
        omega_p = 0;
        omega_q = 0;
        omega_r = 0;
    
        % F/A Matrix
        F = [[0, 0, 0, cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), sin(psi)*(w*cos(phi) + v*sin(phi)) + cos(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)),                                   -cos(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))), - sin(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - cos(psi)*(v*cos(phi) - w*sin(phi)), 1, 0, 0,  0,  0,  0,  0,                    0,                    0]
            [0, 0, 0, cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), sin(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)) - cos(psi)*(w*cos(phi) + v*sin(phi)),                                   -sin(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))),   cos(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - sin(psi)*(v*cos(phi) - w*sin(phi)), 0, 1, 0,  0,  0,  0,  0,                    0,                    0]
            [0, 0, 0,         -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta),                                               cos(theta)*(v*cos(phi) - w*sin(phi)),                                             - sin(theta)*(w*cos(phi) + v*sin(phi)) - u*cos(theta),                                                                                                     0, 0, 0, 1,  0,  0,  0,  0,                    0,                    0]
            [0, 0, 0,                   0,                                        r_m - b_r,                                        b_q - q_m,                                                                                  0,                                                                             -(981*cos(theta))/100,                                                                                                     0, 0, 0, 0, -1,  0,  0,  0,                    w,                   -v]
            [0, 0, 0,           b_r - r_m,                                                0,                                        p_m - b_p,                                                      (981*cos(phi)*cos(theta))/100,                                                                    -(981*sin(phi)*sin(theta))/100,                                                                                                     0, 0, 0, 0,  0, -1,  0, -w,                    0,                    u]
            [0, 0, 0,           q_m - b_q,                                        b_p - p_m,                                                0,                                                     -(981*cos(theta)*sin(phi))/100,                                                                    -(981*cos(phi)*sin(theta))/100,                                                                                                     0, 0, 0, 0,  0,  0, -1,  v,                   -u,                    0]
            [0, 0, 0,                   0,                                                0,                                                0,                  sin(phi)*tan(theta)*(b_r - r_m) - cos(phi)*tan(theta)*(b_q - q_m),               - cos(phi)*(b_r - r_m)*(tan(theta)^2 + 1) - sin(phi)*(b_q - q_m)*(tan(theta)^2 + 1),                                                                                                     0, 0, 0, 0,  0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta)]
            [0, 0, 0,                   0,                                                0,                                                0,                                        cos(phi)*(b_r - r_m) + sin(phi)*(b_q - q_m),                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,            -cos(phi),             sin(phi)]
            [0, 0, 0,                   0,                                                0,                                                0,              (sin(phi)*(b_r - r_m))/cos(theta) - (cos(phi)*(b_q - q_m))/cos(theta), - (cos(phi)*sin(theta)*(b_r - r_m))/cos(theta)^2 - (sin(phi)*sin(theta)*(b_q - q_m))/cos(theta)^2,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta)]
            [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
            [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
            [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
            [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
            [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
            [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
            [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
            [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
            [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]];
         
    
        % G Matrix
        G = [[ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0,   0,   0]
            [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0,   0,   0]
            [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0,   0,   0]
            [-1,  0,  0,  0,                    w,                   -v,   0,   0,   0,   0,   0,   0]
            [ 0, -1,  0, -w,                    0,                    u,   0,   0,   0,   0,   0,   0]
            [ 0,  0, -1,  v,                   -u,                    0,   0,   0,   0,   0,   0,   0]
            [ 0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta),   0,   0,   0,   0,   0,   0]
            [ 0,  0,  0,  0,            -cos(phi),             sin(phi),   0,   0,   0,   0,   0,   0]
            [ 0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta),   0,   0,   0,   0,   0,   0]
            [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0,   0,   0]
            [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0,   0,   0]
            [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0,   0,   0]
            [ 0,  0,  0,  0,                    0,                    0, 100,   0,   0,   0,   0,   0]
            [ 0,  0,  0,  0,                    0,                    0,   0, 100,   0,   0,   0,   0]
            [ 0,  0,  0,  0,                    0,                    0,   0,   0, 100,   0,   0,   0]
            [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0, 100,   0,   0]
            [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0, 100,   0]
            [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0,   0, 100]];
             
        [Phi,Gamma]=c2d(F,G,Ts);
    end
    
    function H = funcLinDisObs(x_linpt,u_linpt,t)
        % Extract values
        x = x_linpt(1);
        y = x_linpt(2);
        z = x_linpt(3);
        u = x_linpt(4);
        v = x_linpt(5);
        w = x_linpt(6);
        phi = x_linpt(7);
        theta = x_linpt(8);
        psi = x_linpt(9);
        v_wxE = x_linpt(10);
        v_wyE = x_linpt(11);
        v_wzE = x_linpt(12);
    
        b_A_x = x_linpt(13);
        b_A_y = x_linpt(14);
        b_A_z = x_linpt(15);
        b_p = x_linpt(16);
        b_q = x_linpt(17);
        b_r = x_linpt(18);
    
        A_x = u_linpt(1);
        A_y = u_linpt(2);
        A_z = u_linpt(3);
        p = u_linpt(4);
        q = u_linpt(5);
        r = u_linpt(6);
    
        H = [[1, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            [0, 1, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            [0, 0, 1,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            [0, 0, 0,                              cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), sin(psi)*(w*cos(phi) + v*sin(phi)) + cos(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)), -cos(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))), - sin(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - cos(psi)*(v*cos(phi) - w*sin(phi)), 1, 0, 0, 0, 0, 0, 0, 0, 0]
            [0, 0, 0,                              cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), sin(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)) - cos(psi)*(w*cos(phi) + v*sin(phi)), -sin(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))),   cos(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - sin(psi)*(v*cos(phi) - w*sin(phi)), 0, 1, 0, 0, 0, 0, 0, 0, 0]
            [0, 0, 0,                                      -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta),                                               cos(theta)*(v*cos(phi) - w*sin(phi)),           - sin(theta)*(w*cos(phi) + v*sin(phi)) - u*cos(theta),                                                                                                     0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
            [0, 0, 0,                                                0,                                                0,                                                0,                                                                                  1,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            [0, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               1,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            [0, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            [0, 0, 0,                        u/(u^2 + v^2 + w^2)^(1/2),                        v/(u^2 + v^2 + w^2)^(1/2),                        w/(u^2 + v^2 + w^2)^(1/2),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            [0, 0, 0,                           -w/(u^2*(w^2/u^2 + 1)),                                                0,                              1/(u*(w^2/u^2 + 1)),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            [0, 0, 0, -(u*v)/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)),      1/((u^2 + w^2)^(1/2)*(v^2/(u^2 + w^2) + 1)), -(v*w)/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]];
    end


    %% CUSUM Test
    y = x_cor(1:N,13:18);
    AoA_k_km1_fault = atan(x_cor(1:N,6)./x_cor(1:N,4));
    Innov_AoA_fault = z_k(1:N,11) - AoA_k_km1_fault;
    y = [y, Innov_AoA_fault];
    y(1:10, 6) = 0; %% remove initial spike
    y(1:10, 7) = 0; %% remove initial spike

    theta_0 = [0.5007    0.3004   -0.0947    0.0052   -0.0070    0.0105   1.8171e-06];
    sigma_0 = [ 3.2143    2.9422    2.9578    0.0027    0.0074    0.0073   4.5582e-04];

    % Thresholds for positive and negative tests
    threshold_pos = [ 20,   50,  100,  20000, 1000,   1000,  500]; % Example thresholds for positive tests for each measurement
    threshold_neg = [-100, -50, -100, -1000, -1000, -1000, -200]; % Example thresholds for negative tests for each measurement
    
    % Initialize variables to count alarms
    n_alarm_pos = zeros(1, size(y, 2)); % number of alarms for positive test
    n_alarm_neg = zeros(1, size(y, 2)); % number of alarms for negative test
    
    % Initialize arrays to store time indices of alarms
    k_alarm_pos = cell(1, size(y, 2)); % time indices k for alarms of positive test
    k_alarm_neg = cell(1, size(y, 2)); % time indices k for alarms of negative test
    
    % Initialize CUSUM variables
    g_pos = zeros(size(y));
    g_neg = zeros(size(y));
    s = (y - theta_0) ./ sigma_0; 

    % Define leakage terms for each signal
    leakage_terms = [0.5, 0.5, 0.5, 0.5, 2.5, 9, 8]; % Example leakage terms for each signal

    for m = 1:size(y, 2)
        for k = 1:size(y, 1)-1
            g_pos(k+1, m) = max(0, g_pos(k, m) + s(k, m))-leakage_terms(m); % positive test
            g_neg(k+1, m) = min(0, g_neg(k, m) + s(k, m))+leakage_terms(m); % negative test
    
            % Check for positive alarms
            if g_pos(k+1, m) > threshold_pos(m)
                n_alarm_pos(m) = n_alarm_pos(m) + 1;
                k_alarm_pos{m} = [k_alarm_pos{m} k+1];
                g_pos(k+1, m) = 0; % reset
            end
    
            % Check for negative alarms
            if g_neg(k+1, m) < threshold_neg(m)
                n_alarm_neg(m) = n_alarm_neg(m) + 1;
                k_alarm_neg{m} = [k_alarm_neg{m} k+1];
                g_neg(k+1, m) = 0; % reset
            end
        end
    end
    
    % Check if k_alarm_pos and k_alarm_neg have elements, if not, assign zero
    if isempty(k_alarm_pos{1}) && isempty(k_alarm_neg{1})
        Ax_f_instance = 0;
    else
        Ax_f_instance = [cell2mat(k_alarm_pos(1)), cell2mat(k_alarm_neg(1))];
    end
    
    if isempty(k_alarm_pos{2}) && isempty(k_alarm_neg{2})
        Ay_f_instance = 0;
    else
        Ay_f_instance = [cell2mat(k_alarm_pos(2)), cell2mat(k_alarm_neg(2))];
    end
    
    if isempty(k_alarm_pos{3}) && isempty(k_alarm_neg{3})
        Az_f_instance = 0;
    else
        Az_f_instance = [cell2mat(k_alarm_pos(3)), cell2mat(k_alarm_neg(3))];
    end
    
    if isempty(k_alarm_pos{4}) && isempty(k_alarm_neg{4})
        p_f_instance = 0;
    else
        p_f_instance = [cell2mat(k_alarm_pos(4)), cell2mat(k_alarm_neg(4))];
    end
    
    if isempty(k_alarm_pos{5}) && isempty(k_alarm_neg{5})
        q_f_instance = 0;
    else
        q_f_instance = [cell2mat(k_alarm_pos(5)), cell2mat(k_alarm_neg(5))];
    end
    
    if isempty(k_alarm_pos{6}) && isempty(k_alarm_neg{6})
        r_f_instance = 0;
    else
        r_f_instance = [cell2mat(k_alarm_pos(6)), cell2mat(k_alarm_neg(6))];
    end
    
    if isempty(k_alarm_pos{7}) && isempty(k_alarm_neg{7})
        AoA_f_instance = 0;
    else
        AoA_f_instance = [cell2mat(k_alarm_pos(7)), cell2mat(k_alarm_neg(7))];
    end

end