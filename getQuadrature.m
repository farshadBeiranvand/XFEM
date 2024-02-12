function [xgp,ygp,wgp,xshear,yshear,wshear]=getQuadrature(Q,p,normEl,heavyEl,heavytipEl,neighbortipEl,tipEl,geom,indice,x1c,y1c,tc,xc,yc)

t90=[-360,-270,-180,-90,0,90,180,270,360];
if ismember(p,normEl)==1
    %type='normEl';
    [xgp, ygp, wgp] = IntPoints2DRefElemQuad(2,2);
    [xshear, yshear, wshear] = IntPoints2DRefElemQuad(2,2);
    if Q==4
        [xshear, yshear, wshear] = IntPoints2DRefElemQuad(1,1);
    end
    
elseif ismember(p,heavyEl)==1
    %type='heavyEl';
    [xcrack,ycrack,xend,yend,xxElem,yyElem]=firstendpoint(Q,geom,indice,x1c,y1c,tc);
    
    if ismember(tc,t90)==1
        [xgp, ygp, wgp] = IntPoints2DLevelSet(xxElem,yyElem,xcrack,ycrack,xend,yend,4,2);
        [xshear, yshear, wshear] = IntPoints2DLevelSet(xxElem,yyElem,xcrack,ycrack,xend,yend,4,2);
        if Q==4
            [xshear, yshear, wshear] = IntPoints2DLevelSet(xxElem,yyElem,xcrack,ycrack,xend,yend,1,1);
        end
    else
        [xgp, ygp, wgp] = IntPoints2DRefElemQuad(51,51);
        [xshear, yshear, wshear] = IntPoints2DRefElemQuad(51,51);
        if Q==4
            [xshear, yshear, wshear] = IntPoints2DLevelSet(xxElem,yyElem,xcrack,ycrack,xend,yend,1,1);
        end
    end
    
    
elseif ismember(p,neighbortipEl)==1
    %type='neighbortipEl';
    [xgp, ygp, wgp] = IntPoints2DRefElemQuad(5,5);
    [xshear, yshear, wshear] = IntPoints2DRefElemQuad(5,5);
    if Q==4
        [xshear, yshear, wshear] = IntPoints2DRefElemQuad(1,1);
    end
    
elseif ismember(p,heavytipEl)==1
    %type='heavytipEl';
    [xcrack,ycrack,xend,yend,xxElem,yyElem]=firstendpoint(Q,geom,indice,x1c,y1c,tc);
    
    if ismember(tc,t90)==1
        [xgp, ygp, wgp] = IntPoints2DLevelSet(xxElem,yyElem,xcrack,ycrack,xend,yend,4,2);
        [xshear, yshear, wshear] = IntPoints2DLevelSet(xxElem,yyElem,xcrack,ycrack,xend,yend,4,2);
        if Q==4
            [xshear, yshear, wshear] = IntPoints2DLevelSet(xxElem,yyElem,xcrack,ycrack,xend,yend,1,1);
        end
    else
        [xgp, ygp, wgp] = IntPoints2DRefElemQuad(51,51);
        [xshear, yshear, wshear] = IntPoints2DRefElemQuad(51,51);
        if Q==4
            [xshear, yshear, wshear] = IntPoints2DLevelSet(xxElem,yyElem,xcrack,ycrack,xend,yend,1,1);
        end
    end
    
elseif ismember(p,tipEl)==1
    %type='tipEl';
    [xcrack,ycrack,xend,yend,xxElem,yyElem]=firstendpoint(Q,geom,indice,x1c,y1c,tc);
    [xgp, ygp, wgp] = IntPoints2DCrackTipElem(xxElem,yyElem, xcrack,ycrack,xend,yend, xc, yc,3);
    [xshear, yshear, wshear] = IntPoints2DCrackTipElem(xxElem,yyElem, xcrack,ycrack,xend,yend, xc, yc,3);
    if Q==4
        [xshear, yshear, wshear] = IntPoints2DCrackTipElem(xxElem,yyElem, xcrack,ycrack,xend,yend, xc, yc,1);
    end
    
end


end




