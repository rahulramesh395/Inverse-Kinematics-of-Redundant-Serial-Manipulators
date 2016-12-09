%Kinematics of a 3R planar serial manipulator
%Pseudo Inverse Method with Null Space Included and Joint Limits
clc;clear all;close all;
syms a1 a2 a3 ad; % Joint angles
syms xd yd x y; %position of the end effector;
L1=input('Enter the length of Link 1 '); % Length of Link 1
L2=input('Enter the length of Link 2 '); % Length of Link 2
L3=input('Enter the length of Link 3 '); % Length of Link 3
xd=input('Enter the desired position of the end effector '); %Desired position of the end effector
R=L1+L2+L3;  % length of the manipulator in fully extended position
f=sqrt(xd(1,1)^2+xd(2,1)^2); % Location of the given coordinates 
if (f>R)
    disp('The Desired Point lies outside the workspace\n'); % To check if the desired coordinate of the end effector lles outside the workspace 
else
x_1=L1*cos(a1)+L2*cos(a1+a2)+L3*cos(a1+a2+a3);  % Equation of x interms of the joint angles
y_1=L1*sin(a1)+L2*sin(a1+a2)+L3*sin(a1+a2+a3);  % Equation of x interms of the joint angles
J_1=[diff(x_1,a1) diff(x_1,a2) diff(x_1,a3); diff(y_1,a1) diff(y_1,a2) diff(y_1,a3)] %jacbian matrix
a_min=input('Enter min '); %Minimum Joint Limit of the manipulator
a_max=input('Enter max '); %Maximum Joint Limit of the manipulator
a=[a_min a_min a_min]; %initializing the minimum angle to be the initial joint angle
a_d=[0 0 0];
x_2=subs(x_1,[a1,a2,a3],[a(1,1),a(1,2),a(1,3)]); % Substituting the initial joint angles to x equation
y_2=subs(y_1,[a1,a2,a3],[a(1,1),a(1,2),a(1,3)]); % Substituting the initial joint angles to y equation
fprintf('The initial position of the end effector is ');
x_2=eval(x_2)  % printing the initial x coordinate of the end effector 
y_2=eval(y_2)  % printing the initial y coordinate of the end effector
x=[x_2;y_2];
e=abs(xd-x); %error
sum=0; 
syms a_t1 a_t2 a_t3; %temporary variable
a_t=[a_t1 a_t2 a_t3]; %converting the temporary variables to matrix form 
K=0.01; % positive gradient 
qd=a_max-a_min; %rangle of the joint angles
qc2=(a_max+a_min)/2; %midpoint of the joint angle limits
for i=1:3
    sum=sqrt(sum+(((K*(a_t(1,i)-qc2)/qd))^2));  %cost function for JOINT LIMIT AVOIDANCE
end

for y=1:3
    Q(1,y)=diff(sum,a_t(1,y)); % differentiating the cost function wrt to joint angles  
end
Q_2=matlabFunction(Q,'vars',[a_t1 a_t2 a_t3]); %converting into matlabFunction
I=[1 0 0;0 1 0;0 0 1]; % Indentity matrix (3x3)
Q_1=[0 0 0];

% To plot the initial orientation of the end effector

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
      
while(e(1,1)>0.00001 || e(2,1)>0.00001)   % while error is very less the loop runs
  
    Q_1=subs(Q,[a_t1,a_t2,a_t3],[a(1,1),a(1,2),a(1,3)]); % Passing the values of joint angles into the function after updating
    
    for w=1:3
        Q_1=Q_2(a(1,1),a(1,2),a(1,3));
        J=subs(J_1,[a1,a2,a3],[a(1,1),a(1,2),a(1,3)]);  %passing the values of joint angles to the Jacobian Matrix
        J=eval(J); % Evaluating the Jacobian Matrix
        J_i=pinv(J); %Finding the Pseudoinverse of the jacobian
        ad = (J_i)*e-(I-(J_i*J))*Q_1'; % Change in Joint angles
        x_2=subs(x_1,[a1,a2,a3],[ad(1,1)+a(1,1),ad(2,1)+a(1,2),ad(3,1)+a(1,3)]); % Substituting the new value of Joint angles to find the new x coordinate of the end effector
        y_2=subs(y_1,[a1,a2,a3],[ad(1,1)+a(1,1),ad(2,1)+a(1,2),ad(3,1)+a(1,3)]); % Substituting the new value of Joint angles to find the new y coordinate of the end effector
        x_2=eval(x_2); % evaluating x coordinate
        y_2=eval(y_2); % evaluating y coordinate
        x=[x_2;y_2];
        e=abs(xd-x);   % Finding the new error
        a(1,1)=ad(1,1)+a(1,1); % Updating the Joint angles
        a(1,2)=ad(2,1)+a(1,2); % Updating the Joint angles
        a(1,3)=ad(3,1)+a(1,3); % Updating the Joint angles
    end
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
pause(0.5);

% To covert the solution in radians to Degree

for i=1:3
     a(1,i)=(a(1,i)*180)/pi;
     number=a(1,i)/360;
     integ=fix(number);
     fract=abs(number-integ);
     a_d(1,i)=fract*360
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