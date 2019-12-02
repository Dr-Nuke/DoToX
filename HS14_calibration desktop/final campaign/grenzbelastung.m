% Matlab-Skript zur Darstellung des Problems der hohen Grenzbelastungen bei
% Geringverdienern und die Abmilderung durch "weichere"
% Transferleistungsreduktionen bei steigendem Einkommen

Emin=0;         % minimales betrachtetes Einkommen
Emax=50000;     % maximales betrachtetes Einkommen

KG = 340*12;    % jährliches Kindergeld für zwei Kinder
KGEmax =28092;  % maximales Einkommen, bis wohin Kindergeld gezahlt wird

nE=101;         %Anzahl der betrachteten Einkommensschritte
dE=(Emax-Emin)/(nE-1); % Einkommensänderung pro Schritt

E=linspace(Emin,Emax,nE); % Basis-Eenkommen

% Wohngeld in Abhängigkeit vom Eenkommen
EKind=KG*(E<KGEmax); %klassisch "heavyside"

Breite=7000; % Breite der Fehlerfunktion
EKindNeu=0.5*KG*(-erf((E-KGEmax)/Breite)+1); % neu, kontinuierlich



% plot
pf=1000; %Plot-Faktor zur Normierung auf 1000€
fig=figure(1);clf;
fig.Position=[100 162 453 775];
subplot(3,1,1)
plot(E/pf,EKind,'DIsplayname','bisher')
hold on
plot(E/pf,EKindNeu,'DIsplayname','neu')
grid on

xlabel('Brutto in [1000 €]')
ylabel('Netto')
legend()
title('Kindergeld in Abhängigkeit vom Bruttoeinkommen')



subplot(3,1,2)
plot(E/pf,(E+EKind)/pf,'DIsplayname','bisher')
hold on
plot(E/pf,(E+EKindNeu)/pf,'DIsplayname','neu')
grid on

xlabel('Brutto in [1000 €]')
ylabel('Netto in [1000 €]')
legend()
title('Einkommen incl. Kindergeld')


subplot(3,1,3)
GrenzBelastung=(-gradient(E+EKind)/dE+1)*100;
plot(E/pf,GrenzBelastung,'DIsplayname','bisher')
hold on
GrenzBelastungNeu=(-gradient(E+EKindNeu)/dE+1)*100;
plot(E/pf,GrenzBelastungNeu,'DIsplayname','neu')

grid on
xlabel('Brutto in [1000 €]')
ylabel('Grenzbelastung [%]')
legend()
str1='Beitrag des Wohngeldes zur Grenzbelastung';
title(sprintf('%s\n neu maximal %2.0f%s',str1,max(GrenzBelastungNeu),'%'))

 print('Grenzbelastung','-dpdf')
  print('Grenzbelastung','-dpng')