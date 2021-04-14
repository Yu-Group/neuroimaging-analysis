% calculation for the FWHM of the beam that was used to image AOLLSM and
% ExLLSM data
lambda = .488;
ri = 1.33;
NAmin = 0.517;
NAmax = 0.55;
mind = asind(NAmin/ri);
maxd = asind(NAmax/ri);
lambda = 0.488;
alpha = [mind, maxd];
x = -05:.075:05;
z = -05:.075:05;
[EField_xpol, EField_zpol] = Calc_Bessel_EField(alpha, x, z);

%%
EField = EField_xpol;
sv_xy = squeeze(EField(1,:,:) .* conj(EField(1,:,:)) + EField(2,:,:) .* conj(EField(2,:,:)) + EField(3,:,:) .* conj(EField(3,:,:)));
figure, imagesc(sv_xy)

figure
for i = 67
    plot(sv_xy(i,:)*10^4)
    hold on
%     pause
end

x_r = [(1:134)]*.01*lambda;
y_r = sv_xy(i, 1:134)*10^4;

FWHM_manual = 89*.0075*.488
FWHM_fitGa = 2.35* 0.185

Y = z;
sv = sv_xy;
sv = max(sv,[],1)/max(max(sv,[],1));
m = max(sv);
r = Y(sv>=m/2);
vq = interp1(sv(sv>=m/2-0.1 & Y > 0),Y(sv>=m/2-0.1 & Y > 0),0.5);
FWHM = vq*2*lambda;