
%plot output temperatures (return water and exhaust air)
figure
tiledlayout(2,2); hAx = nexttile;
for n = 1:4 %looping through each subsystem
	if n > 1, nexttile; end
	
	plot( t, squeeze( x2(:,n,:) ) - 273.15, 'color', matRed  ), hold on
	ylim([0 30])
	
	set(gca, 'XLimSpec', 'Tight');
	grid on, grid minor
	xlabel('time [s]')
	labelY = sprintf('Temp. %s', '[$^\circ$C]');
	ylabel(labelY)
	leg = { sprintf('T$_%i$', n) };
	legend( leg, 'Location','southeast' );
end

%plot output temperatures (return water and exhaust air)
figure
tiledlayout(2,2); hAx = nexttile;
for n = 1:4 %looping through each subsystem
	if n > 1, nexttile; end
	
	plot( t, squeeze( x1(:,n,:) ) - 273.15, 'color', matBlue  ), hold on
	ylim([0 30])
	
	set(gca, 'XLimSpec', 'Tight');
	grid on, grid minor
	xlabel('time [s]')
	labelY = sprintf('Temp. %s', '[$^\circ$C]');
	ylabel(labelY)
	leg = { sprintf('$\\theta_%i$', n) };
	legend( leg, 'Location','southeast' );
end

sys1 = squeeze(x2(:,2,:) - 273.15)';
min1 = min( sys1, [], 2 );
max1 = max( sys1, [], 2 );


%figure, plot(t, sys1), hold on

%plot(t, min1, 'r', 'linewidth', 2)
%plot(t, max1, 'r', 'linewidth', 2)

%h_fill = fill(X,Y_sys1,[.8 0 0],'edgecolor','none', 'facealpha', '.2');


%%
X = [ t'          ;
			flipud(t') ];

Y_sys1 = [ min1          ;
				   flipud(max1) ];

figure
h_fill = fill(X,Y_sys1,[0 .7 0],'edgecolor','none', 'facealpha', '.2');
hold on

min_idx     = find( ismember(Q_var,Q_min, 'rows'), 1 );
nominal_idx = find( ismember(Q_var,Q_midt,'rows'), 1 ) + 81;
max_idx     = find( ismember(Q_var,Q_max, 'rows'), 1 ) + 81*2;

plot( t, sys1(:,min_idx),     'color', matRed, 'linewidth', 1.4 )
plot( t, sys1(:,nominal_idx), 'color', matRed, 'linewidth', 1.4 )
plot( t, sys1(:,max_idx),     'color', matRed, 'linewidth', 1.4 )

set(gca, 'XLimSpec', 'Tight');
grid on, grid minor
xlabel('time [s]')
labelY = sprintf('Temp. %s', '[$^\circ$C]');
ylabel(labelY)
%leg = { sprintf('$\\theta_%i$', n) };
%legend( leg, 'Location','northeast' );



figure
h_fill = fill(X,Y_sys1,[0 .7 0],'edgecolor','none', 'facealpha', '.2');
hold on

plot( t, sys1(:,:), 'color', matRed)

set(gca, 'XLimSpec', 'Tight');
grid on, grid minor
xlabel('time [s]')
labelY = sprintf('Temp. %s', '[$^\circ$C]');
ylabel(labelY)

















