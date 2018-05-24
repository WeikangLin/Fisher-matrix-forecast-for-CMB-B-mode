%6th after running fisher fortran code, you can copy and paste the fisher
%matrix to the following.

fisher=[ 1.189703E+06     -1.898943E+04     -1.664760E+05     -2.129516E+06      4.463605E+05     -2.329325E+02      5.349484E+13     -1.893503E+05
-1.898943E+04      1.466861E+06     -3.104180E+06      1.680948E+07     -4.373875E+06     -9.680872E+03      7.138396E+14      3.602760E+03
-1.664760E+05     -3.104180E+06      1.042583E+07     -1.851921E+07      1.375469E+07      2.270726E+04     -2.411279E+15      4.379349E+04
-2.129516E+06      1.680948E+07     -1.851921E+07      1.478089E+09     -3.931714E+08     -1.071596E+06      3.560581E+15      4.784659E+05
 4.463605E+05     -4.373875E+06      1.375469E+07     -3.931714E+08      2.079552E+08      5.111470E+05     -2.961896E+15     -1.022139E+05
-2.329325E+02     -9.680872E+03      2.270726E+04     -1.071596E+06      5.111470E+05      1.364885E+03     -5.457212E+12      6.257674E+01
 5.349484E+13      7.138396E+14     -2.411279E+15      3.560581E+15     -2.961896E+15     -5.457212E+12      5.798293E+23     -1.265456E+13
-1.893503E+05      3.602760E+03      4.379349E+04      4.784659E+05     -1.022139E+05      6.257674E+01     -1.265456E+13      4.111772E+04
];

r=0.01;
nu_0=0.035;

convariance=fisher^-1;
d_r=convariance(1,1)^0.5
d_ns=convariance(2,2)^0.5;
d_tau=convariance(3,3)^0.5;
d_ombh2=convariance(4,4)^0.5;
d_omch2=convariance(5,5)^0.5;
d_H=convariance(6,6)^0.5;
d_As=convariance(7,7)^0.5;
d_nu_0=convariance(8,8)^0.5
d_r_ns=convariance(1,2);
d_r_nu_0=convariance(1,8);


%for r and nu_0

a=sqrt((d_r^2+d_nu_0^2)/2+sqrt((d_r^2-d_nu_0^2)^2/4+d_r_nu_0^2));
b=sqrt((d_r^2+d_nu_0^2)/2-sqrt((d_r^2-d_nu_0^2)^2/4+d_r_nu_0^2));
theta=atan(2*d_r_nu_0/(d_nu_0^2-d_r^2))/2;

t=linspace(0,2*pi); 
plot(nu_0+1.52*a*cos(t)*cos(theta)-1.52*b*sin(t)*sin(theta),r+1.52*b*sin(t)*cos(theta)+1.52*a*cos(t)*sin(theta))
hold on
plot(nu_0+2.48*a*cos(t)*cos(theta)-2.48*b*sin(t)*sin(theta),r+2.48*b*sin(t)*cos(theta)+2.48*a*cos(t)*sin(theta))
plot(nu_0+3.44*a*cos(t)*cos(theta)-3.44*b*sin(t)*sin(theta),r+3.44*b*sin(t)*cos(theta)+3.44*a*cos(t)*sin(theta))
axis([0,0.09,0,0.02])
plot(nu_0,r,'*')
x_y_axis=axis;
text(x_y_axis(2)/10,x_y_axis(4)*0.9,'COrE','fontsize',18)
str=['\nu_{0,fid}',num2str(nu_0)];
text(nu_0/1.5,r/5,str,'fontsize',18)
xlabel('\nu_0','fontsize',22)
ylabel('r','fontsize',22)
set(gca,'FontSize',14)