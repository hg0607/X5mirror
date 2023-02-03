clear
close all
format short
file = 'MirrorData2023_2_1_17_51_40_664.csv';

time = readmatrix(file,'Range','E2:E100000');
n = length(time);
msec = zeros(n,1);
for i = 1:1:n-1
    if(time(i)>time(i+1))
        msec(i) = time(i+1)+1000-time(i);
    else
        msec(i) = time(i+1)-time(i);
    end    
end
mean(msec)
std(msec)

X1 = readmatrix(file,'Range','F2:F100000');
Y1 = readmatrix(file,'Range','G2:G100000');
Z1 = readmatrix(file,'Range','H2:H100000');
X2 = readmatrix(file,'Range','I2:I100000');
Y2 = readmatrix(file,'Range','J2:J100000');
Z2 = readmatrix(file,'Range','K2:K100000');
X3 = readmatrix(file,'Range','L2:L100000');
Y3 = readmatrix(file,'Range','M2:M100000');
Z3 = readmatrix(file,'Range','N2:N100000');
q1 = readmatrix(file,'Range','O2:O100000');
q2 = readmatrix(file,'Range','P2:P100000');
q3 = readmatrix(file,'Range','Q2:Q100000');
pos1 = readmatrix(file,'Range','R2:R100000');
pos2 = readmatrix(file,'Range','S2:S100000');
pos3 = readmatrix(file,'Range','T2:T100000');
vel1 = readmatrix(file,'Range','U2:U100000');
vel2 = readmatrix(file,'Range','V2:V100000');
vel3 = readmatrix(file,'Range','W2:W100000');
n = length(X1);
t = 0.033*linspace(1,n,n);
lambda = 0.1;
[b,a] = butter(2,lambda);
L1 = zeros(n,1);
L2 = zeros(n,1);
for i=1:1:n
L1(i) = norm([X1(i)-X2(i), Y1(i)-Y2(i), Z1(i)-Z2(i)]);    
L2(i) = norm([X3(i)-X2(i), Y3(i)-Y2(i), Z3(i)-Z2(i)]);            
end
figure(1)
plot(t,L1,'r--')
hold on
plot(t,L2,'b--')
L1_mean = ones(1,n)*mean(L1);
L1_std = ones(1,n)*std(L1);
L2_mean = ones(1,n)*mean(L2);
L2_std = ones(1,n)*std(L2);
plot(t,L1_mean,'r')
fill([t,fliplr(t)],[L1_mean-L1_std, fliplr(L1_mean+L1_std)], [.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5)
text(t(floor(n/2)),L1_mean(1)+3*L1_std(1),strcat(num2str(L1_mean(1)) , '\pm' , num2str(L1_std(1))))
plot(t,L2_mean,'b')
fill([t,fliplr(t)],[L2_mean-L2_std, fliplr(L2_mean+L2_std)], [.6797 .7031 .9805], 'linestyle', 'none', 'FaceAlpha',0.5)
text(t(floor(n/2)),L2_mean(1)+3*L2_std(1),strcat(num2str(L2_mean(1)) , '\pm' , num2str(L2_std(1))))
title('Arm length')

X12 = X2 - X1;
Y12 = Y2 - Y1;
Z12 = Z2 - Z1;
X13 = X3 - X1;
Y13 = Y3 - Y1;
Z13 = Z3 - Z1;
X23 = X3 - X2;
Y23 = Y3 - Y2;
Z23 = Z3 - Z2;
myq1 = zeros(n,1);
myq2 = zeros(n,1);
myq3 = zeros(n,1);
for i=1:1:n
   v12 = [X12(i); Y12(i); Z12(i)];
   v23 = [X23(i); Y23(i); Z23(i)];
   vtmp1 = cross(v12, [0;-1;0]);
   myq1(i) = atan2d(norm(cross(vtmp1,[-1;0;0])),dot(vtmp1,[-1;0;0]));
   myq2(i) = -atan2d(norm(cross(v12,[0;-1;0])),dot(v12,[0;-1;0])); 
   myq3(i) = -atan2d(norm(cross(v12,v23)),dot(v12,v23)); 
end
figure
subplot(3,1,1)
plot(t,q1,'k');
hold on;
plot(t,myq1,'r')
title('q1')
subplot(3,1,2)
plot(t,q2,'k');
hold on;
plot(t,myq2,'r')
title('q2')
subplot(3,1,3)
plot(t,q3,'k');
hold on;
plot(t,myq3,'r')
title('q3')

X12T = zeros(n,1);Y12T = zeros(n,1);Z12T = zeros(n,1);
X13T = zeros(n,1);Y13T = zeros(n,1);Z13T = zeros(n,1);
for i=1:n
    T2 = MDH(-pi/2,0,0,deg2rad(myq1(i))+pi/2)*MDH(pi/2,0,0,deg2rad(myq2(i)+90))*MDH(-pi/2,mean(L1),0,deg2rad(myq3(i)));
    T3 = T2*MDH(0,mean(L2),0,0);
   
    X12T(i) = T2(1,4);Y12T(i) = -T2(2,4);Z12T(i) = T2(3,4);
    X13T(i) = T3(1,4);Y13T(i) = -T3(2,4);Z13T(i) = T3(3,4);
end
figure
subplot(3,1,1)
plot(X12,'k')
hold on
plot(X12T,'r')
title('X12')
subplot(3,1,2)
plot(Y12,'k')
hold on
plot(Y12T,'r')
title('Y12')
subplot(3,1,3)
plot(Z12,'k')
hold on
plot(Z12T,'r')
title('Z12')

figure
subplot(3,1,1)
plot(X13,'k')
hold on
plot(X13T,'r')
title('X13')
subplot(3,1,2)
plot(Y13,'k')
hold on
plot(Y13T,'r')
title('Y13')
subplot(3,1,3)
plot(Z13,'k')
hold on
plot(Z13T,'r')
title('Z13')

%% 
% Your initial state guess at time k, utilizing measurements up to time k-1: xhat[k|k-1]
initialStateGuess = [0;0;0; 0;0;0; -0.2;-0.27;1.4; 0;0;0; 0.25;0.2]; % xhat[k|k-1]
% Construct the filter
ukf = unscentedKalmanFilter(...
    @vdpNonAdditiveNoiseStateFcn_X5,... % State transition function
    @vdpMeasurementFcn_X5,... % Measurement function
    initialStateGuess,...
    'HasAdditiveProcessNoise',false);

ukf.ProcessNoise = diag([1 1 1 1 1 1 1e-6 1e-6]);
ukf.MeasurementNoise  = 1e-3;
Nsteps = n;
timeVector = 0.033*1:Nsteps;
yMeas = [X1 -Y1 Z1 X2 -Y2 Z2 X3 -Y3 Z3];
Nstates = 14;
e = zeros(Nsteps,9); % Residuals (or innovations)
xCorrectedUKF = zeros(Nsteps,Nstates); % Corrected state estimates
PCorrected = zeros(Nsteps,Nstates,Nstates); % Corrected state estimation error covariances
for k=1:Nsteps
    % Let k denote the current time.
    %
    % Residuals (or innovations): Measured output - Predicted output
    e(k,:) =  yMeas(k,:) - vdpMeasurementFcn_X5(ukf.State); % ukf.State is x[k|k-1] at this point
    % Incorporate the measurements at time k into the state estimates by
    % using the "correct" command. This updates the State and StateCovariance
    % properties of the filter to contain x[k|k] and P[k|k]. These values
    % are also produced as the output of the "correct" command.
    [xCorrectedUKF(k,:), PCorrected(k,:,:)] = correct(ukf,yMeas(k,:));
    % Predict the states at next time step, k+1. This updates the State and
    % StateCovariance properties of the filter to contain x[k+1|k] and
    % P[k+1|k]. These will be utilized by the filter at the next time step.
    predict(ukf);
end

figure();
subplot(3,1,1);
plot(timeVector,rad2deg(xCorrectedUKF(:,1))+180,'r',timeVector,myq1,'k');
ylabel('x_1');
subplot(3,1,2);
plot(timeVector,-rad2deg(xCorrectedUKF(:,2)),'r',timeVector,myq2,'k');
ylabel('x_2');
subplot(3,1,3);
plot(timeVector,-rad2deg(xCorrectedUKF(:,3)),'r',timeVector,myq3,'k');
ylabel('x_3');

figure();
subplot(9,1,1);
plot(timeVector, e(:,1));
subplot(9,1,2);
plot(timeVector, e(:,2));
subplot(9,1,3);
plot(timeVector, e(:,3));
subplot(9,1,4);
plot(timeVector, e(:,4));
subplot(9,1,5);
plot(timeVector, e(:,5));
subplot(9,1,6);
plot(timeVector, e(:,6));
subplot(9,1,7);
plot(timeVector, e(:,7));
subplot(9,1,8);
plot(timeVector, e(:,8));
subplot(9,1,9);
plot(timeVector, e(:,9));
xlabel('Time [s]');
ylabel('Residual (or innovation)');

eStates = deg2rad([myq1-180 -myq2 -myq3])-xCorrectedUKF(:,1:3);
figure();
subplot(3,1,1);
plot(timeVector,eStates(:,1),...               % Error for the first state
    timeVector, sqrt(PCorrected(:,1,1)),'r', ... % 1-sigma upper-bound
    timeVector, -sqrt(PCorrected(:,1,1)),'r');   % 1-sigma lower-bound
ylabel('Error for state 1');
title('State estimation errors');
subplot(3,1,2);
plot(timeVector,eStates(:,2),...               % Error for the second state
    timeVector,sqrt(PCorrected(:,2,2)),'r', ...  % 1-sigma upper-bound
    timeVector,-sqrt(PCorrected(:,2,2)),'r');    % 1-sigma lower-bound
ylabel('Error for state 2');
subplot(3,1,3);
plot(timeVector,eStates(:,3),...               % Error for the second state
    timeVector,sqrt(PCorrected(:,3,3)),'r', ...  % 1-sigma upper-bound
    timeVector,-sqrt(PCorrected(:,3,3)),'r');    % 1-sigma lower-bound
xlabel('Time [s]');
ylabel('Error for state 3');
legend('State estimate','1-sigma uncertainty bound',...
    'Location','Best');
%% 
function x = vdpNonAdditiveNoiseStateFcn_X5(x,u) 
% x = [q1 q2 q3 dq1 dq2 dq3 x1 y1 z1 dx1 dy1 dz1 L1 L2];
% u = [ddq1 ddq2 ddq3 ddx1 ddy1 ddz1 dL1 dL2];
% Euler integration of continuous-time dynamics x'=f(x) with sample time dt
dt = 0.033; % [s] Sample time
x = x + [x(4); u(1); x(5); u(2); x(6); u(3); x(10); u(4); x(11); u(5); x(12); u(6);0;0]*dt ...
+ [u(1); 0; u(2); 0; u(3); 0; u(4); 0; u(5); 0; u(6); 0;0;0]*dt^2/2 ...
+ [0;0;0;0;0;0; 0;0;0;0;0;0; u(7); u(8)];
end

function y = vdpMeasurementFcn_X5(x) 
y = zeros(1,9);
L1 = x(13);
L2 = x(14);
q1 = x(1);
q2 = x(2);
q3 = x(3);
y(1) = x(7);
y(2) = x(8);
y(3) = x(9);
y(4) = x(7) + L1*sin(q1)*sin(q2);
y(5) = (x(8) + L1*cos(q2));
y(6) = x(9) + L1*cos(q1)*sin(q2);
y(7) = x(7) + L1*sin(q1)*sin(q2) - L2*(cos(q1)*sin(q3) - cos(q3)*sin(q1)*sin(q2));
y(8) = (x(8) + L1*cos(q2) + L2*cos(q2)*cos(q3));
y(9) = x(9) + L2*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + L1*cos(q1)*sin(q2);
end

function T = MDH(alpha,a,d,theta)

T = (trotx(alpha)*transl([a,0,0])*transl([0,0,d])*trotz(theta));
end
