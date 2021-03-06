%% Crop and save MatLAB figure as PDF
function savePDF(filename)
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,'-dpdf','-painters','-r600','-bestfit',filename);
end