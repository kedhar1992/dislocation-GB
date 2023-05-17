


clc
clear all

function euler = hkl2g2euler(zdir, xdir, varargin)
% hkl/uvw to rotation matrix and Euler angle

% TILT 112
%zdir = [-1 -1 2];  % Z
%xdir = [1 1 1    20 12 16    35 71 53   11 35 23  5 53 29   0 16 8   -2 22 10     -17 55 19  -14 34 10   -10 14 2   -59 61 1];   % X

% TILT 110
%zdir = [-1 1 0];  % Z
%rotation_axis = [1 1 -2   3 3 -8   5 5 -16   -10 -10 52    -5 -5 46   -5 -5 2   10 10 52   19 19 2   -1 -1 -2    29 29 22   23 23 26];     % Y
%xdir = [1 1 1    4 4 3    8 8 5    -26 -26 -10   -23 -23 -5     -1 -1 -5     -26 -26 10    -1 -1 19    1 1 -1   -11 -11 29     -13 -13 23];  % X

% TWIST
%zdir = [1 1 1];    % Z
%%GBN = [-1 -1 2   -9 -5 14   -1 -3 4     -2 -11 13    -1 -21 22     2 -13 11    -4 1 3     5 -14 9   1 -2 1     9 -14  5     7 -10 3   11 -13 2  21 -22 1];     % Y
%xdir = [-1 1 0   -33 40 -7     -55 39 16     -8	5	3     -15 8 7     -8	3	5     -16	55	-39     -40	7	33      -24	0	24     -33	-7	40      -16	-5	21     -5	-3	8      -8	-7	15];    % X

k = 1;
for i = 1:3:length(zdir)
    for j = 1:3:length(xdir)
        z = zdir(i:(i+2));
        x = xdir(j:(j+2));
        
        y = cross(z,x);
        r_col_1 = x/norm(x);
        r_col_2 = y/norm(y);
        r_col_3 = z/norm(z);
        
        r = [r_col_1' r_col_2' r_col_3'];
        
        % Check the properties of Rotation matrix           (1). det(r) = 1   (2). inv(r) = tranpose(r)
        if det(r) == 1
            if inv(r) == transpose(r)
                disp("Rotation matrix achieved")
            end
        else
            disp("Check/Alter R values")
        end
        
        euler(k,:) = g2eulers(r);
        k = k + 1;
        
    end
end

angle = real(euler(:,1)');

end



