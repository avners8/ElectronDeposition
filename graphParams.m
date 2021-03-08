function graphParams(ptitle, pxlabel, pylabel, pzlabel, twoD)
    grid on; 
    title(ptitle); 
    xlabel(pxlabel); 
    ylabel(pylabel); 
    set(gca, 'FontSize', 14); 
    set(gcf,'color','w'); 
    set(gca,'linewidth',2); 
    legend('show');
    if twoD
        h = colorbar;
        colormap pink;
        set(get(h,'title'),'string',pzlabel,'Rotation',0.0);
    end
end