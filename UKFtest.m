format short


% Your initial state guess at time k, utilizing measurements up to time k-1: xhat[k|k-1]
initialStateGuess = [0;0;0]; % xhat[k|k-1]
% Construct the filter
ukf = unscentedKalmanFilter(...
    @vdpNonAdditiveNoiseStateFcn_,... % State transition function
    @vdpMeasurementFcn_,... % Measurement function
    initialStateGuess,...
    'HasAdditiveProcessNoise',false);
ukf.ProcessNoise = 1*diag([1 1 1]);
ukf.MeasurementNoise  = 1e-4;

Nsteps = 10000;
timeVector = 0.033*1:Nsteps;
yMeas = zeros(Nsteps,6);
v = 0.01;
w = 0.01;
L1 = 1;
L2 = 1;
for i=1:Nsteps
    yMeas(i,1) = v*i+1e-4*randn;
    yMeas(i,3) = v*i + L1*cos(w*i)+1e-4*randn;
    yMeas(i,4) = L1*sin(w*i)+1e-4*randn;
    yMeas(i,5) = v*i + L1*cos(w*i) + L2*cos(2*w*i)+1e-4*randn;
    yMeas(i,6) = L1*sin(w*i) + L2*sin(2*w*i)+1e-4*randn;
end
e = zeros(Nsteps,6); % Residuals (or innovations)
xCorrectedUKF = zeros(Nsteps,3); % Corrected state estimates
PCorrected = zeros(Nsteps,3,3); % Corrected state estimation error covariances
for k=1:Nsteps
    % Let k denote the current time.
    %
    % Residuals (or innovations): Measured output - Predicted output
    e(k,:) =  yMeas(k,:) - vdpMeasurementFcn_(ukf.State); % ukf.State is x[k|k-1] at this point
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
plot(timeVector,xCorrectedUKF(:,1));
ylabel('x_1');
subplot(3,1,2);
plot(timeVector,xCorrectedUKF(:,2));
ylabel('x_2');
subplot(3,1,3);
plot(timeVector,xCorrectedUKF(:,3));
ylabel('x_3');

syms alpha1 alpha2 alpha3 alpha4 a1 a2 a3 a4 d1 d2 d3 d4 q1 q2 q3 q4 L1 L2 real
T2 = simplify(MDH(alpha1,a1,d1,q1+pi/2)*MDH(alpha2,a2,d2,q2+pi/2)*MDH(alpha3,a3,d3,q3));
xyz2 = T2(1:3,4);
subs(xyz2,[alpha1,alpha2,alpha3,d1,d2,d3,a1,a2,a3],...
          [-pi/2,pi/2,-pi/2,    0,0,0,   0,0,L1])
T3 = simplify(T2*MDH(alpha4,a4,d4,q4));
xyz3 = T3(1:3,4);
subs(xyz3,[alpha1,alpha2,alpha3,alpha4, d1,d2,d3,d4, a1,a2,a3,a4],...
          [-pi/2,pi/2,-pi/2,0           0,0,0,0,     0,0,L1,L2])
function T = MDH(alpha,a,d,theta)

T = simplify(trotx(alpha)*transl([a,0,0])*transl([0,0,d])*trotz(theta));
end

%% 
function x = vdpNonAdditiveNoiseStateFcn_(x,u)
% Euler integration of continuous-time dynamics x'=f(x) with sample time dt
dt = 0.033; % [s] Sample time
x = x + [u(1); u(2); u(3)]*dt ;
end

function y = vdpMeasurementFcn_(x)
y = zeros(1,6);
L1 = 1;
L2 = 1;
y(1) = x(1);
y(3) = x(1) + L1*cos(x(2));
y(4) = 0 + L1*sin(x(2));
y(5) = x(1) + L1*cos(x(2)) + L2*cos(x(2)+x(3));
y(6) = 0 + L1*sin(x(2))+ L2*sin(x(2)+x(3));
end