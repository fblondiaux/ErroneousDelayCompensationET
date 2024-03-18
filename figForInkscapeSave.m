%Change font, legend size, lines thickness, ... before saving the figure 

function figForInkscapeSave(F,fullFileName)

constantsPlots;

figure(F);
set(gca,'linewidth',thinLine)
set(gca, 'Color', 'none'); % Sets axes background
savefig(fullFileName)
set(gca, 'box', 'off')
set(findall(gcf,'-property','FontSize'),'FontWeight','Normal','FontSize',6,'fontname','Arial')

% Update the font size of the legend
legendFontSize = 4;
hLegend = findobj(gcf, 'Type', 'legend');
set(hLegend, 'FontSize', legendFontSize);

set(F, 'PaperPositionMode', 'auto');
%xlim([-0.005  0.7550])
saveas(gcf,fullFileName,'svg')
end