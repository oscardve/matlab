function a00638139_ProyectoFinal9
%{
--------------------------------------------------------------------------
DESCRIPCION
Esta rutina calcula los ángulos a los cuales se debe de apuntar un rifle
para que la bala llegue al blanco bajo condiciones adversas de viento,
empleando el método de Monte C

CONSIDERACIONES
1. Se consideró una fuerza de fricción de aire proporcional a la velocidad.
2. El blanco se encuentra a 900m en dirección z. 
3. La gravedad actúa en dirección -y. 
4. La fricción del aire actúa en las direcciones -z y -y. 
5. La fuerza del viento se calculó mediante la ecuación de fricción del 
   aire, como si fuera la bala la que se mueve en lugar del aire. 
   La función de viento representa a la variación de la constante de 
   fricción. El viento actúa en dirección +-x, y varía a lo largo de z 
   y en el tiempo.

AUTOR
Oscar David Villarreal <oscardvillarreal@gmail.com>
Versión 1.00 - Noviembre, 2011
--------------------------------------------------------------------------
%}

%--------------------------------------------------------------------------
%CONSTANTES Y VALORES INICIALES DE LAS VARIABLES Y GRÁFICAS:

clear,clc,clf %borra los valores antiguos

%Constantes modificables:
k = 0.001; %cte. de arrastre
ciclos = 10;
mod = pi/1000;
T = 0.8; %Temperatura para el factor de Boltzmann
tolerancia = 0.11;

%Constantes no modificables:
distBlanco = 900;
m = 0.15;
g = 9.81;
v = 800;
zvector = 0:0.01:distBlanco;

%Valores iniciales de las variables a optimizar:
theta = 27*pi/14000 + (rand - 0.5) * 13*pi/9500; %ángulo entre los ejes z->y
phi = -pi/2000 + (rand - 0.5) * 3*pi/1000; %ángulo entre los ejes z->x
angmin = -pi/2000 - 3*pi/1000;
angmax = 27*pi/14000 + 13*pi/9500;

%Valores iniciales de las demás variables:
t = 0;
z0 = 0;
y0 = 0;
x0 = 0;

%Gráfica inicial de la función y(x):
subplot(2,3,1)
hold on
grafy = plot(x0,y0, 'b*');
grafy2 = plot(6,6, 'r*'); %para comparar con el valor previo
grid on
set(ezplot('x^2+y^2=1'), 'color', 'b')
ezplot('x^2+y^2=4')
ezplot('x^2+y^2=9')
ezplot('x^2+y^2=16')
ezplot('x^2+y^2=25')
plot(0,0,'bo')
axis([-5,5,-5,5])
title('Blanco') 
xlabel('x [m]')
ylabel('y [m]')
hold off

%Gráfica inicial de la función y(z):
subplot(2,3,2)
hold on
grafz = plot(z0,y0, 'b*');
grafz2 = plot(-1,-1, 'r*'); %para comparar con el valor previo
grafx = plot(z0,x0, 'y*');
grafx2 = plot(-1,-1, 'r*');
grid on
axis([0,distBlanco,-5,5])
title('Trayecto: en y (azul) y en x (amarillo)') 
xlabel('z [m]')
ylabel('y,x [m]')
hold off

%Gráfica inicial de los valores de los ángulos:
subplot(2,3,3)
hold on
plot(0,theta,'b*');
plot(0,phi,'g*');
axis([0,ciclos,angmin,angmax])
grid on
xlabel('# de iteración')
title('Theta (azul) y Phi (verde)')
hold off

%Gráfica inicial de la función de la constante k del viento:
vientox = 1000;
ax0 = 5;
bx = 1000000;
cx = 450;
kx = ax0 .* exp(-1./bx .* (zvector - vientox.*t - cx).^2);
subplot(2,3,4)
hold on
kvientox = plot(zvector,kx);
grid on
axis([0,distBlanco,-1,5])
xlabel('z [m]')
title('Viento')
hold off

%Gráfica inicial de la magnitud del error:
subplot(2,3,5)
grid on
xlabel('# de iteración')
title('Magnitud del error')

%Anota los ángulos finales:
subplot(2,3,6)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
title('Ángulos óptimos de disparo')

%Agiliza la graficación:
set(kvientox,'EraseMode','xor')
set(grafy,'EraseMode','xor')
set(grafz,'EraseMode','xor')
set(grafx,'EraseMode','xor')

%--------------------------------------------------------------------------
%MÉTODO DE MONTE CARLO:

%Evalúa la efectividad de los ángulos iniciales y registra sus valores:
vz0 = v*cos(theta);
vy0 = v*sin(theta);
vx0 = v*sin(phi);
y = (vy0 + g/k) * distBlanco/vz0 + g/k^2 * log(1 - distBlanco*k/vz0);
tf = -1/k * log(1- distBlanco*k/vz0);
ax = ax0*sin(tf);
kx0 = ax .* exp(-1./bx .* (zvector - vientox.*tf - cx).^2);
x = kx0(114)*tf^2/(2*m) + vx0*tf;
errori = abs(x)+abs(y);
errorii = 0;        

j = 0;
while errori > tolerancia
    j = j+1;

    %Modifica los ángulos, asegurándose de que queden dentro de sus rangos:
    thetaf = theta + (rand - 0.5).*mod;
    thetaf(thetaf<149*pi/266000) = 877*pi/266000;
    thetaf(thetaf>877*pi/266000) = 149*pi/266000;    
    phif = phi + (rand - 0.5).*mod;
    phif(phif<-7*pi/2000) = pi/400;
    phif(phif>pi/400) = -7*pi/2000;
    
    %Evalúa la efectividad de los nuevos ángulos:
    vz0f = v*cos(thetaf);
    vy0f = v*sin(thetaf);
    vx0f = v*sin(phif);
    yf = (vy0f + g/k) * distBlanco/vz0f + g/k^2 * log(1 - distBlanco*k/vz0f);
    tf = -1/k * log(1- distBlanco*k/vz0f);
    ax = ax0*sin(tf);
    kx0 = ax .* exp(-1./bx .* (zvector - vientox.*tf - cx).^2);
    xf = kx0(114)*tf^2/(2*m) + vx0f*tf;
    errorf = abs(xf)+abs(yf);
    set(subplot(2,3,3),'NextPlot','add')
    axis([0,j,angmin,angmax])
    set(subplot(2,3,5),'NextPlot','add')
    axis([1,j+1,0,7])
    drawnow

    %Decide si adoptar los nuevos ángulos empleando el factor de Boltzmann:
    if errorf < errori
        theta = thetaf;
        phi = phif;
        set(subplot(2,3,5),'NextPlot','add')
        plot([j,j+1],[errori,errorf])
        errori = errorf;
        y = yf;
        vy0 = vy0f;
        vz0 = vz0f;
        x = xf;
        vx0 = vx0f;
        set(subplot(2,3,3),'NextPlot','add')
        plot(j,theta,'b*')
        plot(j,phi,'g*')      
    else       
        boltzmann = exp(-abs(errorf-errori)/T);
        ran = rand;
        if ran < boltzmann
            theta = thetaf;
            phi = phif;
            set(subplot(2,3,5),'NextPlot','add')
            plot([j,j+1],[errori,errorf])
            errori = errorf;
            y = yf;
            vy0 = vy0f;
            vz0 = vz0f;
            x = xf;
            vx0 = vx0f;
            set(subplot(2,3,3),'NextPlot','add')
            plot(j,theta,'b*')
            plot(j,phi,'g*')     
        else
            set(subplot(2,3,5),'NextPlot','add')
            plot([j,j+1],[errori,errori])
        end        
    end
    
    %----------------------------------------------------------------------
    %EVOLUCIÓN TEMPORAL: Graficar cada "ciclos" iteraciones
    if rem(j, ciclos) == 0 && errori ~= errorii || j==1 || errori < tolerancia

        errorii = errori;        
        t = 0;
        i = 1;
        tf = -1/k * log(1- distBlanco*k/vz0);
        
        while t<tf
            t = t + 0.01;

            %La constante k de fuerza del viento como función del tiempo:
            ax = ax0*sin(t);
            kxvector = ax .* exp(-1./bx .* (zvector - vientox.*t - cx).^2);
            set(kvientox,'XData',zvector,'YData',kxvector)
            kx = kxvector(i);

            %Coordenadas de la bala:      
            x = kx*t^2/(2*m) + vx0*t;
            y = -g*t/k + (k*vy0 + g) / (k^2) * (1- exp(-k*t));
            set(grafy,'XData',x,'YData',y)
            z = (vz0/k) * (1- exp(-k * t));
            set(grafz,'XData',z,'YData',y)
            set(grafx,'XData',z,'YData',x)
            drawnow
            i = i+1;
        end
        set(grafy2,'XData',x,'YData',y)
        set(grafx2,'XData',z,'YData',x)
        set(grafz2,'XData',z,'YData',y)
    end
    %
    if abs(errori) < 0.5
        T = 0.001;
        mod = pi/10000;
    end  
    %}
end

%Escribe los ángulos finales:
theta_letras = num2str(theta);
phi_letras = num2str(phi);
set(subplot(2,3,6),'NextPlot','add')
text(0.1,0.8,'Theta final =');
text(0.1,0.7,theta_letras);
text(0.1,0.5,'Phi final =');
text(0.1,0.4,phi_letras);


%FIN
%--------------------------------------------------------------------------
