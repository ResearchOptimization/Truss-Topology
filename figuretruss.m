figure;
hold on;
gb = 10;

PG = 3; PA = 1;
N1 = 4; N2 = 5;
XX = 1; YY = 2; ZZ = 3;
Aux = PropGeo;

if (size(Mnodo,2)==2)
    % 2D figure
    for ele = 1:size(Melem,1)
        Nod1 = Melem(ele,N1);
        Nod2 = Melem(ele,N2);
        
        Ve = PropGeo(Melem(ele,PG),PA);
        if (Ve < vtol), Ve = 0; end
        
        x1 = Mnodo(Nod1,XX); y1 = Mnodo(Nod1,YY);
        x2 = Mnodo(Nod2,XX); y2 = Mnodo(Nod2,YY);
        
        Le = sqrt((x2-x1)^2 + (y2-y1)^2);
        Ae = Ve/Le;
        
        Aux(Melem(ele,PG),PA) = Ae;
    end
    
    Amax = max(Aux);
    
    for ele = 1:size(Melem,1)
        Nod1 = Melem(ele,N1);
        Nod2 = Melem(ele,N2);
        
        Ae = Aux(Melem(ele,PG),PA);
        
        x1 = Mnodo(Nod1,XX); y1 = Mnodo(Nod1,YY);
        x2 = Mnodo(Nod2,XX); y2 = Mnodo(Nod2,YY);
        
        if (Ae > 0)
            plot( [x1 x2], [y1 y2], 'b-', 'linewidth', gb * sqrt(Ae/Amax) );
        end
    end
    
else
    % 3D figure
    for ele = 1:size(Melem,1)
        Nod1 = Melem(ele,N1);
        Nod2 = Melem(ele,N2);
        
        Ve = PropGeo(Melem(ele,PG),PA);
        if (Ve < vtol), Ve = 0; end
        
        x1 = Mnodo(Nod1,XX); y1 = Mnodo(Nod1,YY); z1 = Mnodo(Nod1,ZZ);
        x2 = Mnodo(Nod2,XX); y2 = Mnodo(Nod2,YY); z2 = Mnodo(Nod2,ZZ);
        
        Le = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2);
        Ae = Ve/Le;
        
        Aux(Melem(ele,PG),PA) = Ae;
    end
    
    Amax = max(Aux);
    
    for ele = 1:size(Melem,1)
        Nod1 = Melem(ele,N1);
        Nod2 = Melem(ele,N2);
        
        Ae = Aux(Melem(ele,PG),PA);
        
        x1 = Mnodo(Nod1,XX); y1 = Mnodo(Nod1,YY); z1 = Mnodo(Nod1,ZZ);
        x2 = Mnodo(Nod2,XX); y2 = Mnodo(Nod2,YY); z2 = Mnodo(Nod2,ZZ);
        
        if (Ae > 0)
            plot3( [x1 x2], [y1 y2], [z1 z2], 'b-', 'linewidth', gb * sqrt(Ae/Amax) );
        end
    end
end

axis equal;
clear Ae Amax Aux Le N1 N2 Nod1 Nod2 PA PG Ve XX YY ZZ ele gb x1 x2 y1 y2 z1 z2;
