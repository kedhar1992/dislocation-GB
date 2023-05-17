% Convert Miller indices to orientation matrix and Euler angles
% Euler angle of tilted/twisted crystal  (?1 ? ?2) ZXZ = (t1 t t2)

%  1 1 1 orient y -1 1 0  orient z  -1 -1 2
% 11: 20 12 16 orient y -40 56 8 orient z  -1 -1 2 
% 16: 35 71 53 orient y -195 123 -36 orient z  -1 -1 2         
% 23: 11 35 23 orient y -93 45 -24 orient z  -1 -1 2
% 33: 5 53 29 orient y -135 39 -48 orient z  -1 -1 2     
% 39: 0 16 8 orient y -40 8 -16 orient z  -1 -1 2
% 44: -2 22 10 orient y -54 6 -24 orient z  -1 -1 2
% 57:  -17 55 19 orient y -129 -15 -72 orient z  -1 -1 2
% 63: -14 34 10 orient y -78 -18 -48 orient z  -1 -1 2  
% 78: -10 14 2 orient y -30 -18 -24 orient z  -1 -1 2     
% 88: -59 61 1 orient y -123 -117 -120 orient z  -1 -1 2
% 
% 
%  1 1 1 orient y 1 1 -2  orient z  -1 1 0
% 8: 4 4 3 orient y 3 3 -8 orient z  -1 1 0
% 11: 8 8 5 orient y 5 5 -16 orient z  -1 1 0  
% 20: -26 -26 -10 orient y -10 -10 52 orient z  -1 1 0  
% 27: -23 -23 -5 orient y -5 -5 46 orient z  -1 1 0 
% 39: -1 -1 -5 orient y -5 -5 2 orient z  -1 1 0
% 50: -26 -26 10 orient y 10 10 52 orient z  -1 1 0 
% 58: -1 -1 19 orient y 19 19 2 orient z  -1 1 0    
% 70: 1 1 -1 orient y -1 -1 -2 orient z  -1 1 0
% 83: -11 -11 29 orient y 29 29 22 orient z  -1 1 0       
% 86: -13 -13 23 orient y 23 23 26 orient z  -1 1 0   
% 
% %  Z               X                 Y
%  1 1 1 orient y -1 1 0  orient z  -1 -1 2
% 8: 9.43: 1 1 1      y  -33 40 -7      z -47 -26 73
% 16: 1 1 1    y -55 39 16      z  -23 -71 94
% 21: 1 1 1   y  -8 5 3	array	-2	-11	13
% 28: 1 1 1   y  -15 8 7			-1 -22 23
% 37: 1 1 1   y  -8 3 5	array	2	-13	11
% 43: 1 1 1   y  -16 55 -39	array	-94	23	71
% 50: 1 1 1   y  -40 7 33	array	26	-73	47
% 60: 1 1 1   y  60 -24	0	24	array	24	-48	24
% 69: 1 1 1   y  -33 -7 40	array	47	-73	26	
% 73: 1 1 1   y  -16 -5 21	array	26	-37	11
% 81: 1 1 1   y  -5 -3 8	array	11	-13	2
% 87: 1 1 1   y  	-8 -7 15	array	22	-23	1


% Idea from TSL : 
% Tilt GBs : [tilt axis] = (hkl) = Z; [GBN] = [uvw] = X
% Twist GBs : [twist axis] = (hkl) = X (but here its Z in the IPF plot);  [uvw] = can be either X or Y  
clc
clear all
k = 1;
d = 1;

% MTEX 
cs = crystalSymmetry('cubic',[3.3,3.3,3.3],'mineral','Tantalum')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TILT 1
% Slip plane = {1 1 0}
rotation_axis = [-1 -1 2];
GBN = [1 1 1    20 12 16    35 71 53   11 35 23  5 53 29   0 16 8   -2 22 10     -17 55 19  -14 34 10   -10 14 2   -59 61 1];

for i = 1:3:length(rotation_axis)
    hkl = rotation_axis(i:(i+2));
    for j = 1:3:length(GBN)
        uvw = GBN(j:(j+2));   
        m(d) = [Miller(uvw(1),uvw(2),uvw(3),cs)];
        % Check orthogonality 
        CosTheta = max(min(dot(hkl,uvw)/(norm(hkl)*norm(uvw)),1),-1);
        ThetaInDegrees = real(acosd(CosTheta));
        if ThetaInDegrees == 90
            t1 = asind((uvw(3)*sqrt(hkl(1)^2 + hkl(2)^2 + hkl(3)^2))/(sqrt(uvw(1)^2 + uvw(2)^2 + uvw(3)^2) * sqrt(hkl(1)^2 + hkl(2)^2)));
            t = acosd(hkl(3)/ sqrt(hkl(1)^2 + hkl(2)^2 + hkl(3)^2));
            t2 = acosd(hkl(2)/ sqrt(hkl(1)^2 + hkl(2)^2)); 
            euler = [t1 t t2];
            euler_angle1(k,1:3) = euler;
            % MTEX 
            ori1(k) = orientation.byEuler([euler(1),euler(2),euler(3)]*degree,cs);
            
            ori = orientation.byEuler([euler(1),euler(2),euler(3)]*degree,cs);
            plotIPDF(discreteSample(ori,100),xvector,'antipodal','points','all','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k'); 
            if d > 0
                str = num2str(d);
                annotate(ori,'label',str,'MarkerSize',1,'MarkerFaceColor','r','MarkerEdgeColor','k')
            end
            d = d + 1;
            k = k + 1;
            hold on
        end
    end
end

% MTEX 
% figure
% plotIPDF(discreteSample(ori1,100),xvector,'antipodal','points','all','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');  
% figure
% plotPDF(ori1,Miller({1,0,0},ori1.CS));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  1 1 1 orient y 1 1 -2  orient z  -1 1 0
% 8: 4 4 3 orient y 3 3 -8 orient z  -1 1 0
% 11: 8 8 5 orient y 5 5 -16 orient z  -1 1 0  
% 20: -26 -26 -10 orient y -10 -10 52 orient z  -1 1 0  
% 27: -23 -23 -5 orient y -5 -5 46 orient z  -1 1 0 
% 39: -1 -1 -5 orient y -5 -5 2 orient z  -1 1 0
% 50: -26 -26 10 orient y 10 10 52 orient z  -1 1 0 
% 58: -1 -1 19 orient y 19 19 2 orient z  -1 1 0    
% 70: 1 1 -1 orient y -1 -1 -2 orient z  -1 1 0
% 83: -11 -11 29 orient y 29 29 22 orient z  -1 1 0       
% 86: -13 -13 23 orient y 23 23 26 orient z  -1 1 0   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
k = 1;
% TILT 2
% Slip plane = {1 1 2}
rotation_axis = [-1 1 0];
%rotation_axis = [1 1 -2   3 3 -8   5 5 -16   -10 -10 52    -5 -5 46   -5 -5 2   10 10 52   19 19 2   -1 -1 -2    29 29 22   23 23 26];   
GBN = [1 1 1    4 4 3    8 8 5    -26 -26 -10   -23 -23 -5     -1 -1 -5     -26 -26 10    -1 -1 19    1 1 -1   -11 -11 29     -13 -13 23];

for i = 1:3:length(rotation_axis)
    hkl = rotation_axis(i:(i+2));
    for j = 1:3:length(GBN)
        uvw = GBN(j:(j+2));  
        m(d) = [Miller(uvw(1),uvw(2),uvw(3),cs)];
        % Check orthogonality 
        CosTheta = max(min(dot(hkl,uvw)/(norm(hkl)*norm(uvw)),1),-1);
        ThetaInDegrees = real(acosd(CosTheta));
        if ThetaInDegrees == 90
            t1 = asind((uvw(3)*sqrt(hkl(1)^2 + hkl(2)^2 + hkl(3)^2))/(sqrt(uvw(1)^2 + uvw(2)^2 + uvw(3)^2) * sqrt(hkl(1)^2 + hkl(2)^2)));
            t = acosd(hkl(3)/ sqrt(hkl(1)^2 + hkl(2)^2 + hkl(3)^2));
            t2 = acosd(hkl(2)/ sqrt(hkl(1)^2 + hkl(2)^2)); 
            euler = [t1 t t2];
            euler_angle2(k,1:3) = euler;
            % MTEX 
            ori2(k) = orientation.byEuler([euler(1),euler(2),euler(3)]*degree,cs);
            
            ori = orientation.byEuler([euler(1),euler(2),euler(3)]*degree,cs);
            plotIPDF(discreteSample(ori,100),xvector,'antipodal','points','all','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k'); 
            if d > 0
                str = num2str(d);
                annotate(ori,'label',str,'MarkerSize',1,'MarkerFaceColor','r','MarkerEdgeColor','k')
            end
            d = d + 1;
            k = k + 1;
            hold on
        end
    end
end

% MTEX 
% plotIPDF(discreteSample(ori2,100),xvector,'antipodal','points','all','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');  
% figure
% plotPDF(ori2,Miller({1,0,0},ori2.CS));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  1 1 1 orient y -1 1 0  orient z  -1 -1 2
% 8:  9.43: 1 1 1      y  -33 40 -7      z -47 -26 73
% 16: 16.43: 1 1 1    y -55 39 16      z  -23 -71 94
% 21: 1 1 1   y  -8	5	3	array	-2	-11	13
% 28: 27.8: 1 1 1   y  -15 8 7			-1 -22 23
% 37: 1 1 1   y  -8	3	5	array	2	-13	11
% 43: 1 1 1   y  -16	55	-39	array	-94	23	71
% 50: 1 1 1   y  -40	7	33	array	26	-73	47
% 60: 1 1 1   y  -24	0	24	array	24	-48	24
% 69: 1 1 1   y  -33	-7	40	array	47	-73	26	
% 73: 1 1 1   y  -16	-5	21	array	26	-37	11
% 81: 1 1 1   y  -5	-3	8	array	11	-13	2
% 87: 1 1 1   y  	-8	-7	15	array	22	-23	1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Z      X      Y
figure
k = 1;
% TWIST 
% Slip plane = {1 1 1} & {1 1 2}
%GBN = [-1 -1 2   -47 -26 73     -23 -71 94    -2 -11 13     -1 -22 23    2 -13 11      -94	23	71     26 -73 47     24 -48 24     47 -73 26	 26	-37	11    11 -13 2    22 -23 1];
rotation_axis = [1 1 1];  
%%GBN = [-1 -1 2   -9 -5 14   -1 -3 4     -2 -11 13    -1 -21 22     2 -13 11    -4 1 3     5 -14 9   1 -2 1     9 -14  5     7 -10 3   11 -13 2    21 -22 1];
%ori = [90 35.26 135; 77 36 241; 112 38 198; 118 41 190;  123 44 183; 131 50 171; 46 54 284; 138 59 160; 141 66 153; 143 73 147; 144 76 145; 144 83 140;  145 88 136];
GBN = [-1 1 0   -33 40 -7     -55 39 16     -8	5	3     -15 8 7     -8	3	5     -16	55	-39     -40	7	33      -24	0	24     -33	-7	40      -16	-5	21     -5	-3	8      -8	-7	15];  


for i = 1:3:length(rotation_axis)
    hkl = rotation_axis(i:(i+2));
    for j = 1:3:length(GBN)
        uvw = GBN(j:(j+2));  
        m(d) = [Miller(uvw(1),uvw(2),uvw(3),cs)];
        % Check orthogonality 
        CosTheta = max(min(dot(hkl,uvw)/(norm(hkl)*norm(uvw)),1),-1);
        ThetaInDegrees = real(acosd(CosTheta));
        if ThetaInDegrees == 90
            t1 = asind((uvw(3)*sqrt(hkl(1)^2 + hkl(2)^2 + hkl(3)^2))/(sqrt(uvw(1)^2 + uvw(2)^2 + uvw(3)^2) * sqrt(hkl(1)^2 + hkl(2)^2)));
            t = acosd(hkl(3)/ sqrt(hkl(1)^2 + hkl(2)^2 + hkl(3)^2));
            t2 = acosd(hkl(2)/ sqrt(hkl(1)^2 + hkl(2)^2)); 
            euler = [t1 t t2];
            euler_angle3(k,1:3) = euler;
            % MTEX 
            ori3(k) = orientation.byEuler([euler(1),euler(2),euler(3)]*degree,cs);
            
            ori = orientation.byEuler([euler(1),euler(2),euler(3)]*degree,cs);
            plotIPDF(discreteSample(ori,100),xvector,'antipodal','points','all','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k'); 
            if d > 33 && d < 38
                str = num2str(d);
                annotate(ori,'label',str,'MarkerSize',1,'MarkerFaceColor','r','MarkerEdgeColor','k')
            end
            d = d + 1;
            k = k + 1;
            hold on
        end
    end
end

% MTEX 
% plotIPDF(discreteSample(ori3,100),yvector,'antipodal','points','all','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');  
% figure
% plotPDF(ori3,Miller({1,0,0},ori3.CS));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Various plots 

% IPF 
figure
plotIPDF(discreteSample(ori1,100),xvector,'antipodal','points','all','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k');  
hold on
plotIPDF(discreteSample(ori2,100),xvector,'antipodal','points','all','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','k');  
hold on
plotIPDF(discreteSample(ori3,100),xvector,'antipodal','points','all','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','k');  
hold off

% PF
figure
plotPDF(ori1,Miller({1,0,0},ori1.CS));
hold on
plotPDF(ori2,Miller({1,0,0},ori2.CS));
hold on
plotPDF(ori3,Miller({1,0,0},ori3.CS));
hold off

% ODF 
figure
annotate(ori1,'MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k')
hold on
annotate(ori2,'MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','k')
hold on
annotate(ori3,'MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



