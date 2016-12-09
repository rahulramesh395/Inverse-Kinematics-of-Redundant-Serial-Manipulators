% Inverse Kinematics of a 3R planar serial manipulator
% Minimization of Joint Angles 
clc;clear all;close all;
j=1;
syms a1 a2 a3; %Joint angles
syms x y; %position of the end effector;
L1=input('Enter the length of Link 1 ');   % Length of Link 1
L2=input('Enter the length of Link 2 ');   % Length of Link 2
L3=input('Enter the length of Link 3 ');   % Length of Link 3
x=input('Enter the x coordinate of the end effector ');  %Desired x coordinate of the end effector
y=input('Enter the y coordinate of the end effetctor '); %Desired Y coordinate of the end effector
R=L1+L2+L3; % length of the manipulator in fully extended position
flag=0;
f=sqrt(x^2+y^2)  % Location of the given coordinates 
if (f>R)
    disp('The Desired Point lies outside the workspace\n'); % To check if the desired coordinate of the end effector lles outside the workspace 
else
    x_1=L1*cos(a1)+L2*cos(a1+a2)+L3*cos(a1+a2+a3);  % Equation of x interms of the joint angles
    y_1=L1*sin(a1)+L2*sin(a1+a2)+L3*sin(a1+a2+a3);  % Equation of y interms of the joint angles
    g1=L1*cos(a1)+L2*cos(a1+a2)+L3*cos(a1+a2+a3)-x  % function relating x and joint angles
    g2=L1*sin(a1)+L2*sin(a1+a2)+L3*sin(a1+a2+a3)-y  % function relating y and joint angles
    f=(a1)^2+(a2)^2+(a3)^2;  %Candidate function for optimization
    T=[diff(f,'a1') diff(g1,'a1') diff(g2,'a1'); diff(f,'a2') diff(g1,'a2') diff(g2,'a2'); diff(f,'a3') diff(g1,'a3') diff(g2,'a3')] %diffrentiating f g1 and g2 wrt joint angles
    h=det(T) % the 3rd equation relating the joint angles obtained using Legrange Multipliers  
    [sola1,sola2,sola3]=solve(g1,g2,h) %%Solving the 3 eqations to find the joint angles 
 
% To Plot the Manipulator Orientation    
 
    a1=sola1;
    a2=sola2;
    a3=sola3;
    a=[a1 a2 a3];
    x1=0:L1*cos(a1)/100:L1*cos(a1);
    y1=x1*tan(a1);
    w1=size(x1,2);
    v1=size(y1,2);
    x2=L1*cos(a1):(L2*cos(a1+a2))/100:(L2*cos(a1+a2)+L1*cos(a1));
    y2=x2*tan(a1+a2)+(y1(v1)-(x1(w1)*tan(a1+a2)));
    w2=size(x2,2);
    v2=size(y2,2);
    x3=(L2*cos(a1+a2)+L1*cos(a1)):(L3*cos(a1+a2+a3))/100:(L1*cos(a1)+L2*cos(a1+a2)+L3*cos(a1+a2+a3));
    y3=(x3*tan(a1+a2+a3))+(y2(v2)-(x2(w2)*tan(a1+a2+a3)));
    plot(x1,y1,'b');
    hold on;
    plot(x1(w1),y1(v1),'*');
    hold on;
    plot(x2,y2,'g');
    hold on;
    w2=size(x2,2);
    v2=size(y2,2);
    plot(x2(w2),y2(v2),'*');
    hold on;
    plot(x3,y3,'r');
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
eval(a_d)
end