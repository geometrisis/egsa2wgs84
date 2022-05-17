format long g
x=444444;
y=4444444;
h=137.944;
h=zeros(size(x));
% Conversion from EGSA Plane coordinates (x,y,h) to EGSA 87 Ellipsoidal
% coordinates (Lat, Lon, Height)
x=x-500000;
k0=0.9996;a = 6378137.000000 ; b = 6356752.314245;
long0=24;
e=sqrt(1-(b/a)^2);
e2=e^2;
e4=e^4;
e6=e^6;
M=y/k0;
mu = M/(a*(1-e2/4-3*e4/64 - 5*e6/256));
e1=(1-(1-e2)^(1/2))/(1+(1-e2)^(1/2));
J1=3*e1/2-27*e1^3/32;
J2=21*e1^2/16-55*e1^4/32;
J3=151*e1^3/96;
J4=1097*e1^4/512;
fp=mu+J1*sin(2*mu)+J2*sin(4*mu)+J3*sin(6*mu)+J4*sin(8*mu);
e2a=e^2/(1-e^2);
C1=e2a.*(cos(fp)).^2;
T1=(tan(fp)).^2;
R1=a*(1-e2)./(1-e2.*(sin(fp).^2)).^(3/2);
N1=a./(1-e2*(sin(fp).^2)).^(1/2);
D=x./(N1.*k0);
Q1=N1.*tan(fp)./R1;
Q2=(D.^2./2);
Q3=(5+3.*T1+10.*C1-4.*C1.^2-9.*e2a).*D.^4/24;
Q4=(61+90.*T1+298.*C1+45.*T1.^2-3.*C1.^2-252.*e2a).*D.^6/720;
Q5=D;
Q6=(1+2.*T1+C1).*D.^3/6;
Q7=(5-2.*C1+28.*T1-3.*C1.^2+8.*e2a+24.*T1.^2).*D.^5/120;
lat_egsa=fp-Q1.*(Q2-Q3+Q4);
long_egsa=(deg2rad(long0)+(Q5-Q6+Q7)./cos(fp));
% Conversion from EGSA 87 Ellipsoidal Coordinates(Lat, Lon, Height) to
% EGSA 87 cartesian coordinates (XYZ)
v=a./sqrt(1-e2*sin(lat_egsa).*sin(lat_egsa));
X_egsa=(v+h).*cos(lat_egsa).*cos(long_egsa);
Y_egsa=(v+h).*cos(lat_egsa).*sin(long_egsa);
Z_egsa=(v.*(1-e2)+h).*sin(lat_egsa);
% Datum shift from EGSA 87 datum to WGS84 datum
X_wgs84=X_egsa-199.723;
Y_wgs84=Y_egsa+74.03;
Z_wgs84=Z_egsa+246.018;
% Conversion from WGS84 cartesian coordinates (XYZ) to WGS84 Ellipsoidal
% Coordinates(Lat, Lon, Height)
elat=1.e-12;
eht=1.e-5;
p=sqrt(X_wgs84.*X_wgs84+Y_wgs84.*Y_wgs84);
lat_wgs84=atan2(Z_wgs84,p./(1-e2));
h_wgs84=0;
dh=1;
dlat=1;
while sum(dlat>elat) || sum(dh>eht)
  lat0=lat_wgs84;
  h0=h_wgs84;
  v=a./sqrt(1-e2.*sin(lat_wgs84).*sin(lat_wgs84));
  h_wgs84=p./cos(lat_wgs84)-v;
  lat_wgs84=atan2(Z_wgs84, p.*(1-e2.*v./(v+h_wgs84)));
  dlat=abs(lat_wgs84-lat0);
  dh=abs(h_wgs84-h0);
end
long_wgs84=atan2(Y_wgs84,X_wgs84);
% Final WGS84 coordinates (Lat, Lon, H). They can be used in Google Earth
LAT_WGS84=180/pi()*lat_wgs84;
LON_WGS84=180/pi()*(long_wgs84);
LAT_WGS84
LON_WGS84
