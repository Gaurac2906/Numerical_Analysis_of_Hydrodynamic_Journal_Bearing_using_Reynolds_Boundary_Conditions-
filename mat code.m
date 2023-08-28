clear all;
clc;

m = 88;
n = 15;

del_theta = 2*pi/m;
del_z = 1/n;
del = del_theta/del_z;

% initialization of variables
E=0.1;   %Varied the value of E to get different curves%
LBD = 1;
theta = zeros( m+2, n+2 );
h_bar = zeros( m+2, n+2);
p_bar = zeros( m+2, n+2 );
orf = 1.7;
p_new = zeros( m+2,n+2 );
itteration = 0;
conv = 1;
old_sum = 0;
err = zeros( m+2, n+2 );
a = zeros( m+2, n+2 );
b = zeros( m+2, n+2 );
e = zeros( m+2, n+2 );
f = zeros( m+2, n+2 );


% plotting film thickness 
for i = 2: m+1 
    for j = 2: n+1
        theta(i,j) = del_theta * ( i - 2 );
        h_bar(i,j) = 1 +E*cos(theta(i,j));
    end
end
%surf(h_bar)
%ylabel('Angle')
%zlabel('Film Thickness')
%title('FILM THICKNESS VARIATION')

% pressure distribution 
while conv > 1e-4
 for i = 2: m+1
     for j = 2: n+1
         p_bar ( m+2,j ) = p_bar ( 2,j ); % closing of mesh showing boundary condition
         p_bar ( i,1 ) = p_bar ( i,3 ); % initialising the p (i,j-1) value in the starting of mesh

         a(i,j) = ( 1 - 3*E*del_theta*sin(theta(i,j))/(2*h_bar(i,j)));
         b(i,j) = ( 1 + 3*E*del_theta*sin(theta(i,j))/(2*h_bar(i,j)));
         e(i,j) = E*sin(theta(i,j))*power(del_theta,2)/( h_bar(i,j)*h_bar(i,j)*h_bar(i,j) );
         f(i,j) = 2*(1+power(LBD*del,2));


         p_new( i,j ) = ( b(i,j)*p_bar(i-1,j) + a(i,j)*p_bar(i+1,j) + power(LBD*del,2)*(p_bar(i,j+1) + p_bar(i,j-1)) + e(i,j) ) / f(i,j);

         err( i,j ) = p_new( i,j ) - p_bar( i,j );

         p_bar( i,j ) = p_bar( i,j ) + err( i,j )*orf;

         if ( p_bar( i,j ) < 0 )
             p_bar( i,j ) = 0;
         end
     end
 end
 new_sum = sum(sum(p_bar));
 conv = ( new_sum - old_sum ) / new_sum ;
 old_sum = new_sum;
 itteration = itteration + 1;
end

surf(p_bar)
ylabel('Angle')
zlabel('Pressure')
title('PRESSURE VARIATION')


 
