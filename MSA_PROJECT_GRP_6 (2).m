clear all; clc; close all;

JointCoord = xlsread("MSA_PROJECT_PROBLEM_DATA_GRP_6.xlsx","Joint Coordinates");
MemConn = xlsread("MSA_PROJECT_PROBLEM_DATA_GRP_6.xlsx","Member Connectivity");

A=3600;
E=2e5;
I=54110085;

ResDoF = [20*3-2,20*3-1,20*3,124*3-1]; % Restrained Degree of Freedom

NumNodes=numel(JointCoord(:,1));
NumMems=numel(MemConn(:,1));

Kunrestrained=zeros(NumNodes*3,NumNodes*3); %Full stiffness matrix including restrained DOFS

%FreeDof = NumNodes*3 - numel(ResDoF); %Total Dof minus restrained DoF's

Forces=zeros(NumNodes*3,1);

Kloc=zeros(6);

for nmem= 1:NumMems
    
    i=MemConn(nmem, 2);
    j=MemConn(nmem, 3);
    w=(1)*MemConn(nmem, 4); % Scaling of loading w

    ix=JointCoord(i, 2); 
    iy=JointCoord(i, 3);
    jx=JointCoord(j,2);
    jy=JointCoord(j,3);

    L=sqrt((jx-ix)^2+(jy-iy)^2);%length
    c=(jx-ix)/L; % cos theta
    s=(jy-iy)/L; % sin theta

    T=[c,s,0,0,0,0; -s,c,0,0,0,0; 0,0,1,0,0,0; 0,0,0,c,s,0; 0,0,0,-s,c,0; 0,0,0,0,0,1];%Transformation Matrix

    %Kloc of member
    Kloc(1,1)= E*A/L; Kloc(1,2)= 0; Kloc(1,3)= 0; Kloc(1,4)= -E*A/L; Kloc(1,5)= 0; Kloc(1,6)= 0;
    Kloc(2,1)= 0; Kloc(2,2)= 12*E*I/(L)^3; Kloc(2,3)= 6*E*I/(L)^2; Kloc(2,4)= 0; Kloc(2,5)= -12*E*I/(L)^3; Kloc(2,6)= 6*E*I/(L)^2;
    Kloc(3,1)= 0; Kloc(3,2)= 6*E*I/(L)^2; Kloc(3,3)= 4*E*I/L; Kloc(3,4)= 0; Kloc(3,5)= -6*E*I/(L)^2; Kloc(3,6)= 2*E*I/L;
    Kloc(4,1)= -E*A/L; Kloc(4,2)= 0; Kloc(4,3)= 0; Kloc(4,4)= E*A/L; Kloc(4,5)= 0; Kloc(4,6)= 0;
    Kloc(5,1)= 0; Kloc(5,2)= -12*E*I/(L)^3; Kloc(5,3)= -6*E*I/(L)^2; Kloc(5,4)= 0; Kloc(5,5)= 12*E*I/(L)^3; Kloc(5,6)= -6*E*I/(L)^2;
    Kloc(6,1)= 0; Kloc(6,2)= 6*E*I/(L)^2; Kloc(6,3)= 2*E*I/L; Kloc(6,4)= 0; Kloc(6,5)= -6*E*I/(L)^2; Kloc(6,6)= 4*E*I/L;

    Kg = T'*Kloc*T;  %Kg of member
    
    GlobDOF(1) = (i-1)*3 + 1; % 3*i-2
    GlobDOF(2) = (i-1)*3 + 2; %3*i-1
    GlobDOF(3) = (i-1)*3 + 3; %3*i
    GlobDOF(4) = (j-1)*3 + 1; %3*j-2
    GlobDOF(5) = (j-1)*3 + 2; %3*j-1
    GlobDOF(6) = (j-1)*3 + 3; %3*j

    %Calculate contribution of member to global (full) stiffness matrix
    for m= 1:6
        for n= 1:6
            
            Kunrestrained(GlobDOF(m),GlobDOF(n))= Kunrestrained(GlobDOF(m),GlobDOF(n)) + Kg(m,n);

        end
    end
    
    %Calculate fixed end forces member wise
    
    Floc= zeros(6,1);

    Floc(2)= w*L/2;
    Floc(3)= w*(L)^2/12;
    Floc(5)= w*L/2;
    Floc(6)= -w*(L)^2/12;
    
    Fglob= T'*Floc;
    
    for p=1:6
    
        Forces(GlobDOF(p))= Forces(GlobDOF(p)) - Fglob(p);

    end

end

%wind loading on member 1,2,11,21,31 of 10 in +x dirn
wind_members = [1,2,11,21,31];
%wind_members = (1)*MemConn(:,1);%Scaling of Wind loading
num_wind_mems = numel(wind_members);

for n_wind_mem = 1:num_wind_mems 

    i= MemConn(n_wind_mem,2);
    j= MemConn(n_wind_mem,3);
    wind_load= MemConn(n_wind_mem,5); 

    ix=JointCoord(i, 2); 
    iy=JointCoord(i, 3);
    jx=JointCoord(j,2);
    jy=JointCoord(j,3);

    L=sqrt((jx-ix)^2+(jy-iy)^2);%length
    c=(jx-ix)/L; % cos theta
    s=(jy-iy)/L; % sin theta

    wind_shear= wind_load*c;
    wind_perpendicular= wind_load*s;

    T=[c,s,0,0,0,0; -s,c,0,0,0,0; 0,0,1,0,0,0; 0,0,0,c,s,0; 0,0,0,-s,c,0; 0,0,0,0,0,1];%Transformation Matrix

    GlobDOF(1) = (i-1)*3 + 1; % 3*i-2
    GlobDOF(2) = (i-1)*3 + 2; %3*i-1
    GlobDOF(3) = (i-1)*3 + 3; %3*i
    GlobDOF(4) = (j-1)*3 + 1; %3*j-2
    GlobDOF(5) = (j-1)*3 + 2; %3*j-1
    GlobDOF(6) = (j-1)*3 + 3; %3*j

    F_wind_loc= zeros(6,1);

    F_wind_loc(1)= -wind_shear*L/2;
    F_wind_loc(2)= wind_perpendicular*L/2;
    F_wind_loc(3)= wind_perpendicular*(L)^2/12;
    F_wind_loc(4)= -wind_shear*L/2;
    F_wind_loc(5)= wind_perpendicular*L/2;
    F_wind_loc(6)= -wind_perpendicular*(L)^2/12;

    F_wind_glob= T'*F_wind_loc;

    for p= 1:6

        Forces(GlobDOF(p))= Forces(GlobDOF(p)) - F_wind_glob(p);

    end

end

Kres= Kunrestrained;
deleted= 0;
for q= 1:numel(ResDoF)
    Kres(ResDoF(q)-deleted,:)= [];
    Kres(:,ResDoF(q)-deleted)= [];
    deleted= deleted + 1;
end

Fres= Forces;
deleted= 0;
for r= 1:numel(ResDoF)
    Fres(ResDoF(r)-deleted,:)= [];
    deleted= deleted + 1;
end

% EXTERNAL CONCENTRATED FORCE & MOMENT 
%numel(ResDoF) is 4
%ResDoF Node No.s = [58 59 60 371]
%          1  2  3   4
% Only for all Node No.s except these ResDof Node No.
% If Node no.(let x) calculated is in ranges:         x<58 58<x<59 59<x<60 60<x<371 x>371
% then deducting no.(let y) from the Node No.(x), y:   0      1       2       3       4
% New Node No. or Node No. in Restrained Matrix Fres is    (x - y) instead of x

%External F(+x)   F(+y)   M(ACW)  are considered +ve
%         3*i-2   3*i-1   3*i

% 500 magnitude of forces and moments with alternating directions 
 Fres(144*3 -4) = Fres(144*3 -4) + 50; % 428 :  M(ACW) on 144

 Fres(135*3-2 -4) = Fres(135*3-2 -4) + (-50); % 399 : F(+x) on 135
 Fres(135*3-1 -4) = Fres(135*3-1 -4) + 50; % 400 : F(+y) on 135
 Fres(135*3 -4) = Fres(135*3 -4) + (-50);     % 401 : M(ACW) on 135

 Fres(124*3 -4) = Fres(124*3 -4) + 50; % 368 :  M(ACW) on 124

 % 371

 Fres(124*3-2 -3) = Fres(124*3-2 -3) + (-50); % 367 : F(+x) on 124

 Fres(100*3-2 -3) = Fres(100*3-2 -3) + 50; % 295 : F(+x) on 100
 Fres(100*3-1 -3) = Fres(100*3-1 -3) + (-50); % 296 : F(+y) on 100
 Fres(100*3 -3) = Fres(100*3 -3) + 50;     % 297 : M(ACW) on 100

 Fres(21*3-2 -3) = Fres(21*3-2 -3) + (-50); % 58 : F(+x) on 21

 % 60

 % 59

 % 58

 Fres(19*3 -0) = Fres(19*3 -0) + 50; % 57 : M(ACW) on 19

 Fres(15*3-2 -0) = Fres(15*3-2 -0) + (-50); % 43 : F(+x) on 15
 Fres(15*3-1 -0) = Fres(15*3-1 -0) + 50; % 44 : F(+y) on 15
 Fres(15*3 -0) = Fres(15*3 -0) + (-50);     % 45 : M(ACW) on 15

 Fres(1*3-2 -0) = Fres(1*3-2 -0) + 50;  % 1 : F(+x) on 1

Dres = Kres\Fres;
Rescount = 1;
Count = 1;
Displacement = zeros(3*NumNodes,1);
for DOF=1:3*NumNodes
      if Rescount <= numel(ResDoF) && DOF == ResDoF(Rescount)
            Displacement(DOF) = 0;
            Rescount= Rescount + 1;
      else
            Displacement(DOF) = Dres(Count);
            Count = Count + 1;
      end
end

% Calculate reactions at restrained DOFs
ReactionForces = Kunrestrained * Displacement - Forces;

% Extract only the reaction forces corresponding to the restrained DOFs
ReactionsAtRestraints = ReactionForces(ResDoF);

% Display the reaction forces and moments
disp('Reaction Forces and Moments at Restrained DOFs:');
for i = 1:numel(ResDoF)
    fprintf('DOF %d: Reaction = %.2f\n', ResDoF(i), ReactionsAtRestraints(i));
end

% Display the displacement of important periphery joints of Howrah Bridge
disp('Displacement x & y && Rotation θ of important periphery joints of Howrah Bridge in 10^(-6) units\n\n');
perphery_Node_Nos = [1, 2, 20, 21, 52, 53, 72, 73, 92, 93, 124, 125, 143, 144];
% Initialize data as an empty matrix to store results
data = [];
for i = 1:numel(perphery_Node_Nos)
    node = perphery_Node_Nos(i);
    fprintf('Node %d: Displacement x: %d Displacement y: %d Rotation θ in degrees: %d\n', ...
        node, Displacement(3*node - 2), Displacement(3*node - 1),(180/pi)*Displacement(3*node));
     % Extract the displacements and rotation for the current node
    disp_x = Displacement(3*node - 2) * 1e6; % Convert to 10^(-6) units
    disp_y = Displacement(3*node - 1) * 1e6;
    rotation = (180/pi)*Displacement(3*node) * 1e6;
    
    % Append row to data array with formatted values
    data = [data; node, disp_x, disp_y, rotation];
end

% Create a table with proper column names
resultTable = array2table(data, ...
    'VariableNames', {'Node', 'Displacement_x_10^-6', 'Displacement_y_10^-6', 'Rotation_theta_10^-6'});

% Display the table
disp(resultTable);

scale = 1000000; % Overal Scaling of Deflection
    DisplacedCoords = zeros(NumNodes,3);
    for nodenumber = 1:NumNodes
        DisplacedCoords(nodenumber,1)= nodenumber;
        DisplacedCoords(nodenumber,2)= JointCoord(nodenumber,2) + scale*Displacement(3*nodenumber-2); %x
        DisplacedCoords(nodenumber,3)= JointCoord(nodenumber,3) + scale*Displacement(3*nodenumber-1); %y
    end
    for memnumber=1:NumMems
        i=MemConn(memnumber,2);%Start Node No.
        j=MemConn(memnumber,3);%End Node No.
        orgxi = JointCoord(i,2);% x-coordinate of original start node
        orgyi = JointCoord(i,3);% y-coordinate of original start node
        orgxj = JointCoord(j,2);% x-coordinate of original end node
        orgyj = JointCoord(j,3);% y-coordinate of original end node
        x = [orgxi,orgxj];%Range/variation of x coordinate from original start node to end node
        y = [orgyi,orgyj];%Range/variation of y coordinate from original start node to end node
        plot(x,y, 'BLACK-', 'LineWidth', 0.1)
        hold on
        finxi = DisplacedCoords(i,2);% x-coordinate of final start node
        finyi = DisplacedCoords(i,3);% y-coordinate of final start node
        finxj = DisplacedCoords(j,2);% x-coordinate of final end node
        finyj = DisplacedCoords(j,3);% y-coordinate of final end node
        x = [finxi,finxj];%Range/variation of x coordinate from final start node to end node
        y = [finyi,finyj];%Range/variation of y coordinate from final start node to end node
        plot(x,y,'r-', 'LineWidth', 0.1);
        hold on
    end

    xlim([-50 450])
    ylim([-10 50])
    Fall = Kunrestrained*Displacement;

    % Label the plot
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('Original and Displaced Structure');
legend('Original Structure', 'Displaced Structure');
grid on;
axis equal;
hold off;