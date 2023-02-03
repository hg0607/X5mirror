
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
% 500Hz
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
fc = 3;
fs = 60;
[b,a] = butter(2,fc/(fs/2));
% X1 = filter(b,a,X1);
% X2 = filter(b,a,X2);
% X3 = filter(b,a,X3);
% Y1 = filter(b,a,Y1);
% Y2 = filter(b,a,Y2);
% Y3 = filter(b,a,Y3);
% Z1 = filter(b,a,Z1);
% Z2 = filter(b,a,Z2);
% Z3 = filter(b,a,Z3);
% X1 = medfilt1(X1,9);
% Y1 = medfilt1(Y1,9);
% Z1 = medfilt1(Z1,9);
% X2 = medfilt1(X2,9);
% Y2 = medfilt1(Y2,9);
% Z2 = medfilt1(Z2,9);
% X3 = medfilt1(X3,9);
% Y3 = medfilt1(Y3,9);
% Z3 = medfilt1(Z3,9);

% %3d trajectory
% figure('NumberTitle','off','Name','3d trajectory');   
% for i = 1:1:n    
%     plot3([X1(i),X2(i)],[Y1(i),Y2(i)],[Z1(i),Z2(i)],'R','markersize',10);
%     axis([-1 0 -0.2 0.5 1.0 1.8]);
%     xlabel('x');ylabel('y');zlabel('z');
%     hold on; 
%     grid on;
%     plot3([X2(i),X3(i)],[Y2(i),Y3(i)],[Z2(i),Z3(i)],'B','markersize',10);
%     if i >= 30
%         plot3(X1(i-29:i), Y1(i-29:i), Z1(i-29:i),'R--')
%         plot3(X2(i-29:i), Y2(i-29:i), Z2(i-29:i),'G--')
%         plot3(X3(i-29:i), Y3(i-29:i), Z3(i-29:i),'B--')
%     end
%     hold off; 
% pause(1/30);
% end

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
q1_ = rad2deg(atan2(sqrt(X12.^2+Z12.^2),-Y12)); 
% q1_ = rad2deg(acos(-Y12./sqrt(X12.^2+Y12.^2+Z12.^2))); 
% q2_ = rad2deg(atan2(-Z12,-Y12./cosd(q1_))); 
q2_ = rad2deg(acos(Z12./sqrt(X12.^2+Z12.^2))); 

myq1 = zeros(n,1);
for i=1:1:n
   v12 = [X12(i); Y12(i); Z12(i)];
   vtmp = cross(v12, [0;-1;0]);
   myq1(i) = atan2d(norm(cross(vtmp,[-1;0;0])),dot(vtmp,[-1;0;0]));
    
end




tmpY3 = -(X13.*cosd(q1_) + Y13.*sind(q1_));
tmpX3 = Y13.*cosd(q1_).*sind(q2_) - Z13.*cosd(q2_) - X13.*sind(q1_).*sind(q2_);
q3_ = rad2deg(atan2(tmpY3,tmpX3)); 

tmpY4 = tmpY3./sind(q3_);
tmpX4 = X13.*sind(q1_).*cosd(q2_) - Y13.*cosd(q1_).*cosd(q2_) - Z13.*sind(q2_) - L1_mean(1);
q4_ = rad2deg(atan2(tmpY4,tmpX4)); 



figure
plot(t,q1,'k');
hold on;
plot(t,pos1,'b')
plot(t,-q2_+180,'r')
plot(t,myq1,'m--')
title('q1')

figure
plot(t,q2,'k');
hold on;
plot(t,-pos2,'b')
plot(t,-q1_,'r')
title('q2')

figure
plot(t,q3,'k');
hold on;
plot(t,-pos3,'b')
plot(t,-q4_,'r')
title('q3')

figure
plot([0;diff(q1)*30],'r')
hold on
plot(vel1,'k')
plot((q1-pos1),'b')
vel1d = [0;diff(q1)*30];
figure
plot([0;diff(vel1)*30],'r')
amax = 0.25*360*1000000/524288; %degree/s^2

q1d = zeros(n,1);q2d = zeros(n,1);q3d = zeros(n,1);
for i=1:1:n
q1d(i) = CalcAngleFromKinectOriginalData([X2(i),Y2(i),Z2(i)], [X1(i),Y1(i),Z1(i)], [ 0,0,-1 ], [ 0,0,0 ]);
q2d(i) = CalcAngleFromKinectOriginalData([X2(i),Y2(i),Z2(i)], [X1(i),Y1(i),Z1(i)], [ 0,-1,0 ], [ 0,0,0 ]);
q3d(i) = CalcAngleFromKinectOriginalData([X2(i),Y2(i),Z2(i)], [X1(i),Y1(i),Z1(i)], [X3(i),Y3(i),Z3(i)], [X2(i),Y2(i),Z2(i)]);
end

figure()
plot(myq1,'k')
hold on
plot(q1d,'r')
plot(q1,'b')
[c,lags] = xcorr(pos1, q1);
figure()
stem(lags,c)

close all
Plant = tf([1],[1 0]);
Ts = 0.033;
MV = struct('Min',-60,'Max',60);
p = 10;
m = 2;
acc1 = 171.66; % 0.25 
MPCobj = mpc(Plant,Ts, p, m, [], MV);
% MPCobj.W.OutputVariables = 10;
% open_system('velsim.slx')
% out = sim('velsim.slx');
% 
% [c,lags] = xcorr(out.MPC, q1);
% figure()
% stem(lags,c)
% title('xcorr of MPC')
% 
% [c,lags] = xcorr(out.ff, q1);
% figure()
% stem(lags,c)
% title('xcorr of feedforward')

function ang = CalcAngleFromKinectOriginalData(start1, end1, start2, end2)
        v1 = MakeVector3d(start1, end1);
        v2 = MakeVector3d(start2, end2);
        ang = RadianOfTwoVector(v1, v2);
end

function y = MakeVector3d(start1, end1)
    y = [start1(1) - end1(1),start1(2) - end1(2),start1(3) - end1(3)];   
end

function ang = RadianOfTwoVector(v1, v2)
        n1 = norm(v1);
        n2 = norm(v2);
        ang = acosd(dot(v1,v2) / (n1 * n2));
end

