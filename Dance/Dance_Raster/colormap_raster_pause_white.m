%% create colour map %%%%
Blue = [0 0 1];
Cyan = [0 1 1];
White = [1 1 1];
Yellow = [1 1 0];
Red = [1 0 0];
Dblue = [0 0 0.5];
Dred = [0.5 0 0];
step1 = 20;
step2 = 60;
Grad2 = [linspace(Blue(1),Cyan(1),step1+step2);linspace(Blue(2),Cyan(2),step1+step2);linspace(Blue(3),Cyan(3),step1+step2)];
Grad3 = [linspace(Yellow(1),Red(1),step2);linspace(Yellow(2),Red(2),step2);linspace(Yellow(3),Red(3),step2)];
Grad4 = [linspace(Red(1),Dred(1),step1);linspace(Red(2),Dred(2),step1);linspace(Red(3),Dred(3),step1)];
colourmap =[Grad2'; White; Grad3'; Grad4'];
%%%% create colour map end %%%%