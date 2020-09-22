function [Vec_RF_M, Vec_RF_A] = Polar_nav_sub(Vec)

Angle = [0 pi/4 pi/2 3*pi/4 pi 5*pi/4 6*pi/4 7*pi/4 2*pi];

Vec_Y=0; Vec_X=0; 

for p=1:8
    Vec_Y= Vec_Y + Vec(p)*sin(Angle(p));
    Vec_X= Vec_X + Vec(p)*cos(Angle(p));
end



if (Vec_Y>0 && Vec_X>0) || (Vec_Y<0 && Vec_X>0)
    Vec_RF_A = atan(Vec_Y/Vec_X);
elseif (Vec_Y<0 && Vec_X<0)  || (Vec_Y>0 && Vec_X<0)
    Vec_RF_A = atan(Vec_Y/Vec_X)+pi;
end

Vec_RF_M = sqrt((Vec_Y*Vec_Y) + (Vec_X*Vec_X));


end