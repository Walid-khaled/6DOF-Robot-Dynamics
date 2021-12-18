clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%% Symbolic Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 6;
syms d1 a1 a2 b1 b2 b4z b5x b5y b5z b6z gs...
    m1 m2 m3 m4 m5 m6...
    q1 q2 q3 q4 q5 q6...
    dq1 dq2 dq3 dq4 dq5 dq6...
    Ixx1 Ixy1 Ixz1 Iyx1 Iyy1 Iyz1 Izx1 Izy1 Izz1...
    Ixx2 Ixy2 Ixz2 Iyx2 Iyy2 Iyz2 Izx2 Izy2 Izz2...
    Ixx3 Ixy3 Ixz3 Iyx3 Iyy3 Iyz3 Izx3 Izy3 Izz3...
    Ixx4 Ixy4 Ixz4 Iyx4 Iyy4 Iyz4 Izx4 Izy4 Izz4...
    Ixx5 Ixy5 Ixz5 Iyx5 Iyy5 Iyz5 Izx5 Izy5 Izz5...
    Ixx6 Ixy6 Ixz6 Iyx6 Iyy6 Iyz6 Izx6 Izy6 Izz6

x1 = 0;
y1 = 0;
z1 = b1;

x2 = a1*cos(q1)+b2*cos(q1+q2);
y2 = a1*sin(q1)+b2*sin(q1+q2);
z2 = d1;

x3 = a1*cos(q1)+a2*cos(q1+q2);
y3 = a1*sin(q1)+a2*sin(q1+q2);
z3 = d1-q3;

x4 = x3;
y4 = y3;
z4 = z3-b4z;

x5 = x4+b5x;
y5 = y3+b5y;
z5 = z4-b5z;

x6 = x4;
y6 = y4;
z6 = z5-b6z;

%%%%%%%%%%%%%%%%%%%%%%%% Inertia or Mass Materix %%%%%%%%%%%%%%%%%%%%%%%%%%

%Linear Velocity Jacobians COM
for i=1:n
    eval(strcat('J',string(i),'v = sym(zeros(3,n));'))
    for j=1:n
    eval(strcat('J',string(i),'v(1,',string(j),') = diff(x',string(i),',q',string(j),');'))%J1v(1,j) = diff(x1,qj)
    eval(strcat('J',string(i),'v(2,',string(j),') = diff(y',string(i),',q',string(j),');'))%J1v(2,j) = diff(y1,qj)
    eval(strcat('J',string(i),'v(3,',string(j),') = diff(z',string(i),',q',string(j),');'))%J1v(3,j) = diff(z1,qj)
    end
    eval(strcat('J',string(i),'v;'))
end

%Jacobian of angular velocity of COM
Jw = sym(zeros(3,n));
for i=1:n
    if i == 3
        Jw((1:3),i)=[0;0;0];%prismatic joint
    else
        Jw((1:3),i)=[0;0;1];%revolute joint
    end
    eval(strcat('J',string(i),'w = Jw;'))%Jiw = Jw
    
end

%Evaluating M matrix
M = sym(zeros(n,n));
for i=1:n
    eval(strcat('R',string(i),' = [cos(q',string(i),') -sin(q',string(i),') 0; sin(q',string(i),') cos(q',string(i),') 0; 0 0 1];'))
    eval(strcat('I',string(i),' = [Ixx',string(i),' Ixy',string(i),' Ixz',string(i),'; Iyx',string(i),' Iyy',string(i),' Iyz',string(i),'; Izx',string(i),' Izy',string(i),' Izz',string(i),'];'))
    eval(strcat('M',string(i),' = m',string(i),'*transpose(J',string(i),'v)*J',string(i),'v+transpose(J',string(i),'w)*R',string(i),'*I',string(i),'*transpose(R',string(i),')*J',string(i),'w;'))
    eval(strcat('M = M+M',string(i),';'))
end
M = simplify(M);
fprintf('M matrix symbolic:\n')
disp(M)

%%%%%%%%%%%%%%%%%%%%%%%%%%% Coriolis forces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Evaluating C matrix
C = sym(zeros(n,n));
for i=1:n
    for j=1:n
        eval(strcat('m',string(i),string(j),' = M(',string(i),',',string(j),');')) %mij = M(i,j)
        eval(strcat('c',string(i),string(j),' = sym(0);'))%cij = sym(0)
        for k=1:n
            eval(strcat('m',string(i),string(k),' = M(',string(i),',',string(k),');'))%mik = M(i,k)
            eval(strcat('m',string(j),string(k),' = M(',string(j),',',string(k),');'))%mjk = M(j,k)
           
            eval(strcat('c',string(i),string(j),string(k),' = 0.5*(diff(m',string(i),string(j),',q',string(k),')+diff(m',string(i),string(k),',q',string(j),')-diff(m',string(j),string(k),',q',string(i),'));')) 
            eval(strcat('c',string(i),string(j),' = c',string(i),string(j),' + c',string(i),string(j),string(k),'*dq',string(k),';'))%cij = cij + cijk*dqk;
        end
        eval(strcat('C(',string(i),',',string(j),') = c',string(i),string(j),';'))%C(i,j) = cij
    end
end
C = simplify(C);
fprintf('C matrix symbolic:\n')
disp(C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gravity forces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Evaluating g vector
g0 = sym([0;-gs;0]);
g = sym(zeros(n,1));
for i=1:n
    eval(strcat('g',string(i),' = sym(0);'))%gi = sym(0)
    for k=1:n
        eval(strcat('g',string(i),' = g',string(i),'-transpose(J',string(k),'v((1:3),',string(i),'))*m',string(k),'*g0;'))%gi = gi-transpose(Jkvi((1:3),k))*mk*g0
    end
    eval(strcat('g(',string(i),') = g',string(i),';'))%g(i) = gi
end
g = simplify(g);
fprintf('g vector symbolic:\n')
disp(g)

%%%%%%%%%%%%%%%%%%%%%%%%% Substitution Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%

d1 = 0.6;
a1 = 0.1;
a2 = 0.4;
b1 = 0.35;
b2 = 0.2;
b4z = 0.1;
b5x = 0.1;
b5y = 0.1;
b5z = 0.1;
b6z = 0.1;
gs = 9.81;
m1 = 3; 
m2 = 2;
m3 = 1.5;
m4 = 1;
m5 = 0.9;
m6 = 0.8;
q1 = deg2rad(50);
q2 = deg2rad(30);
q3 = 0.15;
q4 = deg2rad(90);
q5 = deg2rad(-30);
q6 = deg2rad(10);
dq = [0.2; 0.3; 0.15; 0.1; 0.22; 0.18];
dq1 = dq(1);
dq2 = dq(2);
dq3 = dq(3);
dq4 = dq(4);
dq5 = dq(5);
dq6 = dq(6);
Izz1 = 0.01;
Izz2 = 0.015;
Izz3 = 0.001;
Izz4 = 0.0015;
Izz5 = 0.02;
Izz6 = 0.002;

M = subs(M);%substitue in symbolic form 
C = subs(C);
g = subs(g);
M = double(M);%convert from sym to double
C= double(C);
g = double(g);

fprintf('M matrix:\n')
disp(M)
fprintf('C matrix:\n')
disp(C)
fprintf('g matrix:\n')
disp(g)

%%%%%%%%%%%%%%%%%%%%%%%%%%% Torque Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Dynamics Equation: M(q)ddq + C(q,dq)dq + g(q) = T\n')
ddq = [0.1; 0.15; 0.1; 0.06; 0.07; 0.12];
T = M*ddq + C*dq + g; 
fprintf('Desired Torques:\n')
disp(T)
