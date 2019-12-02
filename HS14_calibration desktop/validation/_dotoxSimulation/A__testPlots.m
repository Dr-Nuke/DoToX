


% plot cross-section values
figure
plot(macroXSenergies,macroXSall,'LineWidth',3); grid on; xlim([30 inf]); ylim([0 3.5])
xlabel('Energy [keV]');ylabel('attenuation [1/cm]');legend('CHCl3', 'H2O', 'Al');

% plot source data
figure;
E = sourceDataAll(:,1);output = sourceDataAll(:,2:end);
plot(E,output);grid on;xlabel('Energy [keV]');ylabel('#/keV/cm^2/mAs @ 1 meter');
legend('0','0.2','0.5','1.0','1.5','2.0')