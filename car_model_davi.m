% Universidade Federal de Santa Maria
% Dinâmica de Estruturas e Aeroelasticidade
% Nome: Davi 
% Sobrenome: Lima Mendes dos Santos
% Matrícula: 202012282

%%
close all;
m = 840;
mf = 53;
m1 = mf;
m2 = mf;
mr = 76;
m3 = mr;
m4 = mr;
Iy = 820;
Iz = 1100;
sig = 0.5;
kf = 10000;
kr = 13000;
ktf = 200000;
ktr = 200000;
cf = 2*sig*sqrt(kf*m);
cr = 2*sig*sqrt(kr*m);
a1 = 1.4;
b1 = 0.7;
a2 = 1.47;
b2 = 0.75;

 y1 = 0;
 y2 = 0;
 y3 = 0;
 y4 = 0;
%%
M = [m, 0, 0, 0, 0, 0, 0;
     0, Iy, 0, 0, 0, 0, 0;
     0, 0, Iz, 0, 0, 0, 0;
     0, 0, 0, m1, 0, 0, 0;
     0, 0, 0, 0, m2, 0, 0;
     0, 0, 0, 0, 0, m3, 0;
     0, 0, 0, 0, 0, 0, m4]
 
 C = [2*cf+2*cr, b1*cf-b2*cf-b2*cr+b1*cr, -a1*cf-a1*cf+a2*cr+a2*cr, -cf, -cf, -cr, -cr;
     cf*b1-cf*b2-cr*b2+cr*b1, b1^2*cf+b2^2*cf+b2^2*cr+b1^2*cr, -a1*b1*cf+a1*b2*cf-a2*b2*cr+a2*b1*cr, -b1*cf, b2*cf, b2*cr, -b1*cr;
     -cf*a1-cf*a1+cr*a2+cr*a2, -a1*b1*cf+a1*b2*cf-a2*b2*cr+a2*b1*cr, 2*a1^2*cf+2*a2^2*cr, a1*cf, a1*cf, -a2*cr, -a2*cr;
     -cf, -b1*cf, a1*cf, cf, 0, 0, 0;
     -cf, b2*cf, a1*cf, 0, cf, 0, 0;
     -cr, b2*cr, -a2*cr, 0, 0, cr, 0;
     -cr, -b1*cr, -a2*cr, 0, 0, 0, cr]
 
  K = [2*kf+2*kr, b1*kf-b2*kf-b2*kr+b1*kr, -2*a1*kf+2*a2*kr, -kf, -kf, -kr, -kr;
     b1*kf-b2*kf-b2*kr+b1*kr, b1^2*(kf+kr)+b2^2*(kf+kr), kf*(-a1*b1+a1*b2)+kr*(-a2*b2+a2*b1), -b1*kf, b2*kf, b2*kr, -b1*kr;
     -a1*kf-a1*kf+a2*kr+a2*kr, kf*(-a1*b1+a1*b2)+kr*(-a2*b2+a2*b1), 2*a1^2*kf+2*a2^2*kr, a1*kf, a1*kf, -a2*kr, -a2*kr;
     -kf, -b1*kf, a1*kf, ktf+kf, 0, 0, 0;
     -kf, b2*kf, a1*kf, 0, ktf+kf, 0, 0;
     -kr, b2*kr, -a2*kr, 0, 0, ktr+kr, 0;
     -kr, -b1*kr, -a2*kr, 0, 0, 0, ktr+kr]
 
 F = [0;
     0;
     0;
     y1*ktf;
     y2*ktf;
     y3*ktr;
     y4*ktr];

 %%
 
 close all;

t0 = 0;         %[s]
v_kmh = 30;
v = v_kmh/3.6; %[m/s]
D = 10; %[m]
A1 = 1; %[m]
A2 = 0.5; %[m]
B1 = 6e-02; %[m]
B2 = -8e-02; %[m]

% -------------------------------------------------------------------------
% Simulação 1: tff = (D+A1)/v; tfr = (D+a1+a2+A1)/v;
tff = (D+A1)/v;
tfr = (D+a1+a2+A1)/v;

t0f = D/v;
t0r = (D+a1+a2)/v;
tf = tfr + 4;
tf = round(tf, 3);

deltat = 0.001;  %[s]
t = t0:deltat:tf;
points = tf/deltat;
%% 
% Simulação 1: Y = B1; lambda = A1*2;
Y = B1;
lambda = A1*2;

f = v/lambda;
omega = 2*pi*f;

X0 = zeros(7,1);
X0_dot = zeros(7,1);

X0_2dot = M\(F(:,1) - C*X0_dot - K*X0);
Xneg1 = (deltat^2)*X0_2dot/2 - deltat*X0_dot + X0;
X = zeros(7, points);
X(:,1) = X0;

Y_vec = zeros(4, points);
i = 1;

% o 'while' a seguir só é válido para a 
% 1ª simulação 

while i < (points+1)
    
    t_now = t(1,i);
    
    if i == 1
        Y_vec(:,i) = 0;
        X(:,i+1) = ((M/deltat^2) + (C/(2*deltat)))\( F(:,1) - (K - 2*M/(deltat^2))*X0 -  (M/(deltat^2) - C/(2*deltat))*Xneg1 );
    else
        if t0f<t_now && t_now < tff
            y1 = Y*sin(omega*(t_now - t0f));
            y2 = y1;
            y3 = 0;
            y4 = y3;
            F = [0; 0; 0; y1*ktf; y2*ktf; y3*ktr; y4*ktr];
            Y_vec(1,i) = y1;
            Y_vec(2,i) = y2;
            Y_vec(3,i) = y3;
            Y_vec(4,i) = y4;
            X(:, i+1) = inv((M/deltat^2) + (C/(2*deltat)))*( F(:,1) - (K - 2*M/deltat^2)*X(:,i) -  (M/deltat^2 - C/(2*deltat))*X(:,i-1) );
        elseif t0r<t_now && t_now<tfr        
            y1 = 0;
            y2 = y1;
            y3 = Y*sin(omega*(t_now - t0r));
            y4 = y3;
            F = [0; 0; 0; y1*ktf; y2*ktf; y3*ktr; y4*ktr];
            Y_vec(1,i) = y1;
            Y_vec(2,i) = y2;
            Y_vec(3,i) = y3;
            Y_vec(4,i) = y4;
            X(:, i+1) = inv((M/deltat^2) + (C/(2*deltat)))*( F(:,1) - (K - 2*M/deltat^2)*X(:,i) -  (M/deltat^2 - C/(2*deltat))*X(:,i-1) );
        else
            y1 = 0;
            y2 = y1;
            y3 = 0;
            y4 = y3;
            F = [0; 0; 0; y1*ktf; y2*ktf; y3*ktr; y4*ktr];
            Y_vec(:,i) = 0;
            X(:, i+1) = inv((M/deltat^2) + (C/(2*deltat)))*( F(:,1) - (K - 2*M/deltat^2)*X(:,i) -  (M/deltat^2 - C/(2*deltat))*X(:,i-1) );
        end
        
    end
    
    i=i+1;
end

Y_vec(:,i) = 0;
%% 
% SIMULATION 1 PLOTS

% plotting x
figure;
subplot(2,2,1);
plot(t, X(1,:), 'b-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('x(t) [m]');
grid on;
% plotting phi
subplot(2,2,2);
plot(t, X(2,:), 'r-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('\phi(t) [rad]');
grid on;
% plotting theta
subplot(2,2,3);
plot(t, X(3,:), 'g-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('\theta(t) [rad]');
grid on;
sgtitle('Simulation 1');

% plotting x1 
figure; 
subplot(2,2,1);
plot(t, X(4,:),'Color', [0.5, 0.2, 0.8]); % purple
xlabel('t [s]');
ylabel('x_1(t) [m]');
grid on;
% plotting x2
subplot(2,2,2);
plot(t, X(5,:), 'r-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('x_2(t) [m]');
grid on;
% plotting x3
subplot(2,2,4);
plot(t, X(6,:), 'g-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('x_3(t) [m]');
grid on;
% plotting x4
subplot(2,2,3);
plot(t, X(7,:), 'b-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('x_4(t) [m]');
grid on;
sgtitle('Simulation 1');

% plotting y1, y2, y3 and y4
figure;
subplot(2,2,1);
plot(t, Y_vec(1,:), 'b-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('y_1(t) [m]');
grid on;
subplot(2,2,2);
plot(t, Y_vec(2,:), 'b-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('y_2(t) [m]');
grid on;
subplot(2,2,4);
plot(t, Y_vec(3,:), 'r-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('y_3(t) [m]');
grid on;
subplot(2,2,3);
plot(t, Y_vec(4,:), 'r-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('y_4(t) [m]');
grid on;
sgtitle('Simulation 1');

%% ------------------------------------------------------------------------
% Simulação 2: tff = (D+A2)/v; tfr = (D+a1+a2+A2)/v;
tff = (D+A2)/v;
tfr = (D+a1+a2+A2)/v;

t0f = D/v;
t0r = (D+a1+a2)/v;
tf = tfr + 4;
tf = round(tf, 3);

deltat = 0.001;  %[s]
t = t0:deltat:tf;
points = tf/deltat;
%% 
% Simulação 2: Y = B2; lambda = A2*2;
Y = B2;
lambda = A2*2;

f = v/lambda;
omega = 2*pi*f;

X0 = zeros(7,1);
X0_dot = zeros(7,1);

X0_2dot = M\(F(:,1) - C*X0_dot - K*X0);
Xneg1 = (deltat^2)*X0_2dot/2 - deltat*X0_dot + X0;
X = zeros(7, points);
X(:,1) = X0;

Y_vec = zeros(4, points);
i = 1;

% o 'while' a seguir só é válido para a 
% 2ª simulação

while i < (points+1)
    
    t_now = t(1,i);
    
    if i == 1
        Y_vec(:,i) = 0;
        X(:,i+1) = ((M/deltat^2) + (C/(2*deltat)))\( F(:,1) - (K - 2*M/(deltat^2))*X0 -  (M/(deltat^2) - C/(2*deltat))*Xneg1 );
    else
        if t0f<t_now && t_now < tff
            y1 = Y*sin(omega*(t_now - t0f));
            y2 = 0;
            y3 = 0;
            y4 = 0;
            F = [0; 0; 0; y1*ktf; y2*ktf; y3*ktr; y4*ktr];
            Y_vec(1,i) = y1;
            Y_vec(2,i) = y2;
            Y_vec(3,i) = y3;
            Y_vec(4,i) = y4;
            X(:, i+1) = inv((M/deltat^2) + (C/(2*deltat)))*( F(:,1) - (K - 2*M/deltat^2)*X(:,i) -  (M/deltat^2 - C/(2*deltat))*X(:,i-1) );
        elseif t0r<t_now && t_now<tfr        
            y1 = 0;
            y2 = 0;
            y3 = 0;
            y4 = Y*sin(omega*(t_now - t0r));
            F = [0; 0; 0; y1*ktf; y2*ktf; y3*ktr; y4*ktr];
            Y_vec(1,i) = y1;
            Y_vec(2,i) = y2;
            Y_vec(3,i) = y3;
            Y_vec(4,i) = y4;
            X(:, i+1) = inv((M/deltat^2) + (C/(2*deltat)))*( F(:,1) - (K - 2*M/deltat^2)*X(:,i) -  (M/deltat^2 - C/(2*deltat))*X(:,i-1) );
        else
            y1 = 0;
            y2 = 0;
            y3 = 0;
            y4 = 0;
            F = [0; 0; 0; y1*ktf; y2*ktf; y3*ktr; y4*ktr];
            Y_vec(:,i) = 0;
            X(:, i+1) = inv((M/deltat^2) + (C/(2*deltat)))*( F(:,1) - (K - 2*M/deltat^2)*X(:,i) -  (M/deltat^2 - C/(2*deltat))*X(:,i-1) );
        end
        
    end
    
    i=i+1;
end

Y_vec(:,i) = 0;

%%
% SIMULATION 2 PLOTS

% plotting x
figure;
subplot(2,2,1);
plot(t, X(1,:), 'b-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('x(t) [m]');
grid on;
% plotting phi
subplot(2,2,2);
plot(t, X(2,:), 'r-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('\phi(t) [rad]');
grid on;
% plotting theta
subplot(2,2,3);
plot(t, X(3,:), 'g-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('\theta(t) [rad]');
grid on;
sgtitle('Simulation 2');

% plotting x1 
figure; 
subplot(2,2,1);
plot(t, X(4,:),'Color', [0.5, 0.2, 0.8]); % purple
xlabel('t [s]');
ylabel('x_1(t) [m]');
grid on;
% plotting x2
subplot(2,2,2);
plot(t, X(5,:), 'r-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('x_2(t) [m]');
grid on;
% plotting x3
subplot(2,2,4);
plot(t, X(6,:), 'g-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('x_3(t) [m]');
grid on;
% plotting x4
subplot(2,2,3);
plot(t, X(7,:), 'b-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('x_4(t) [m]');
grid on;
sgtitle('Simulation 2');

% plotting y1, y2, y3 and y4
figure;
subplot(2,2,1);
plot(t, Y_vec(1,:), 'b-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('y_1(t) [m]');
grid on;
subplot(2,2,2);
plot(t, Y_vec(2,:), 'b-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('y_2(t) [m]');
grid on;
subplot(2,2,4);
plot(t, Y_vec(3,:), 'r-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('y_3(t) [m]');
grid on;
subplot(2,2,3);
plot(t, Y_vec(4,:), 'r-', 'LineWidth', 1);
xlabel('t [s]');
ylabel('y_4(t) [m]');
grid on;
sgtitle('Simulation 2');

