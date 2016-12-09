%Kinematics of a 3R planar serial manipulator
%Compact formulation Method
clc;clear all;close all;
syms a1 a2 a3; % Joint Angles
L1=input('Enter the length of Link 1 '); % Length of Link 1
L2=input('Enter the length of Link 2 '); % Length of Link 2
L3=input('Enter the length of Link 3 '); % Length of Link 3
xd=input('Enter the desired position of the end effector '); %Desired position of the end effector
R=L1+L2+L3; % length of the manipulator in fully extended position
f=sqrt(xd(1,1)^2+xd(2,1)^2); % Location of the given coordinates 
if (f>R)
    disp('The Desired Point lies outside the workspace\n'); % To check if the desired coordinate of the end effector lles outside the workspace 
else
x_1=L1*cos(a1)+L2*cos(a1+a2)+L3*cos(a1+a2+a3); % Equation of x in terms of the joint angles
y_1=L1*sin(a1)+L2*sin(a1+a2)+L3*sin(a1+a2+a3); % Equation of y in terms of the joint angles
J_1=[diff(x_1,a1) diff(x_1,a2) diff(x_1,a3); diff(y_1,a1) diff(y_1,a2) diff(y_1,a3)] %jacbian matrix
a=input('Enter the initial jont angles of the manipulator '); % Initial Joint angles
x_2=subs(x_1,[a1,a2,a3],[a(1,1),a(1,2),a(1,3)]);  % Substituting the initial joint angles to x equation
y_2=subs(y_1,[a1,a2,a3],[a(1,1),a(1,2),a(1,3)]);  % Substituting the initial joint angles to y equation
fprintf('The initial position of the end effector is '); 
x_2=eval(x_2) % printing the initial x coordinate of the end effector 
y_2=eval(y_2) % printing the initial y coordinate of the end effector 
x=[x_2;y_2];
e=abs(xd-x); %error
I=[1]; %Identity Matrix of size (1x1)
while(e(1,1)>0.0000001 || e(2,1)>0.0000001)   % while error is very less, the loop runs
    
    J=subs(J_1,[a1,a2,a3],[a(1,1),a(1,2),a(1,3)]); %passing the values of joint angles to the Jacobian Matrix
    J=eval(J); % Evaluating Jacobian Matrix
    J1=J(:,1:2); % extracting the square matric from the Jacobian
    J2=J(:,3:3); % Remaining Jacobian Matrix
    N=[((-inv(J1))*J2);I]; % matrix that maps into Null Space
    O=[0]; %Zero Coloumn vector of size(n-m) 
    a_p=[((inv(J1))*e);O]; %particular solution
    ad=a_p+((N*(inv(N'*N)*N'))*(-a_p)); % Change in Joint angles
    x_2=subs(x_1,[a1,a2,a3],[ad(1,1)+a(1,1),ad(2,1)+a(1,2),ad(3,1)+a(1,3)]); % Substituting the new value of Joint angles to find the new x coordinate of the end effector
    y_2=subs(y_1,[a1,a2,a3],[ad(1,1)+a(1,1),ad(2,1)+a(1,2),ad(3,1)+a(1,3)]); % Substituting the new value of Joint angles to find the new y coordinate of the end effector
    x_2=eval(x_2); % evaluating x coordinate
    y_2=eval(y_2); % evaluating y coordinate
    x=[x_2;y_2];
    e=abs(xd-x); % Finding the new error
    a(1,1)=ad(1,1)+a(1,1); % Updating the Joint angles
    a(1,2)=ad(2,1)+a(1,2); % Updating the Joint angles
    a(1,3)=ad(3,1)+a(1,3); % Updating the Joint angles
end
    
% To Plot the Manipulator Orientation 

a_1=a(1,1);
a_2=a(1,2);
a_3=a(1,3);
x1=0:L1*cos(a_1)/100:L1*cos(a_1);
y1=x1*tan(a_1);
w1=size(x1,2);
v1=size(y1,2);
x2=L1*cos(a_1):(L2*cos(a_1+a_2))/100:(L2*cos(a_1+a_2)+L1*cos(a_1));
y2=x2*tan(a_1+a_2)+(y1(v1)-(x1(w1)*tan(a_1+a_2)));
w2=size(x2,2);
v2=size(y2,2);
x3=(L2*cos(a_1+a_2)+L1*cos(a_1)):(L3*cos(a_1+a_2+a_3))/100:(L1*cos(a_1)+L2*cos(a_1+a_2)+L3*cos(a_1+a_2+a_3));
y3=(x3*tan(a_1+a_2+a_3))+(y2(v2)-(x2(w2)*tan(a_1+a_2+a_3)));
plot(x1,y1);
hold on;
plot(x1(w1),y1(v1),'*');
hold on;
plot(x2,y2);
hold on;
w2=size(x2,2);
v2=size(y2,2);
plot(x2(w2),y2(v2),'*');
hold on;
plot(x3,y3);
hold on;
w3=size(x3,2);
v3=size(y3,2);
plot(x3(w3),y3(v3),'*');
hold on;

% To covert the solution in radians to Degree

for i=1:3
     a(1,i)=(a(1,i)*180)/pi;
     number=a(1,i)/360;
     integ=fix(number);
     fract=abs(number-integ);
     a_d(1,i)=fract*360;
end
for q=1:3
if a_d(1,q)<0
   a_d(1,q)=abs(a_d(1,q)); 
    
end
if a_d(1,q)>180
  a_d(1,q)=360-a_d(1,q)
elseif a_d<180 & a_d>90
  a_d(1,q)=180-a_d(1,q)
end
end
a_d
end
 
 