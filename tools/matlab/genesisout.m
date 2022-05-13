%author: Biaobin Li
%date:   2022/05/13

classdef genesisout < handle
    properties
        t
        current
        power_t
        z

        s
        power_s
        pulse_energy

        lamda
        amp

    end

    methods
        function obj=genesisout(path)
            if nargin==0
                path='.';
            end

            file1 = [path '/tip.txt'];
            file2 = [path '/spe.txt'];
            file3 = [path '/lamamp.txt'];

            tip = importdata(file1);
            spe = importdata(file2);
            lamamp = importdata(file3);

            obj.t = tip(:,1);              %[fs]
            obj.current = tip(:,2);        %[A]
            obj.power_t = tip(:,3);        %[W]
            tmp = 2.998e8*tip(:,1)*1e15; 
            obj.z = tmp-mean(tmp);         %[m]

            obj.s = spe(:,1);              %[m]
            obj.power_s = spe(:,2);        %[W]
            obj.pulse_energy = spe(:,3);   %[J]

            obj.lamda = lamamp(:,1);       %[nm]
            obj.amp   = lamamp(:,2);       % 1
        end
    end
end





