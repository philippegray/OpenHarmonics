%PWRCONVERTERS  class type
%
%Only one object of "pwrconverters" type is created per simulation study.
%This "pwrconverters" object contains each 
%converter object modelled for the tested distribution network.
%
%"pwrconverters" has two properties: converterArr and numConverters
%
%converterArr: is an array of "converter" type objects. Each element in
%this array contains a converter object that the user has chosen to be
%included in the simulated network. The converter object must be specified
%in the inputs.m file - not in the OpenDSS circuit file.
%
%numConverters: is an integer variable containing the number of converter
%objects modelled in the system.
%

classdef pwrconverters
    properties
        converterArr %--- array of type "converter"
        numConverters%--- integer variable%
    end
    methods
        %constructor - initializes: the number of converters to 0 and
        %the array containing the converter objects to an empty array of 
        %100 "converter" type elements
        function conv_objs = pwrconverters()
            conv_objs.numConverters = 0;
            conv_objs.converterArr = converter.empty(100,0);
        end
        %function - adds a converter object to converterArr with all
        %parameters initialized to their correct values. OpenDSS code is
        %updated to include a monitor, spectrum, and isource object at 3
        %nodes - each corresponding to one of the three phases
        %corresponding to this particular converter object.
        %'name' = arbitrary name that the user chooses for that converter -
        %each object MUST have a unique name.
        %'bus' = the PCC bus name. must exactly match a name of a bus in
        %the OpenDSS circuit file.
        %'Srated' = the rated kVA of the converter
        %'Vdc' = For the VSC: the dc-link voltage in kV for the vsc, and
        %for line-commutated converters: the dc-side emf voltage source in
        %kV
        %'P' = the real power absorbed by the converter in kW
        %'Q' = the reactive power absorbed by the converter kvar
        %Also, the numConverters count is increased by 1 when this function
        %is called, because a new converter is being added to the system!
        function conv_objs = add(conv_objs,type,name,bus,Srated,Vdc,P,Q)
            conv_objs.numConverters = conv_objs.numConverters + 1; %increasing system converter count by 1
            conv_objs.converterArr(conv_objs.numConverters) = ...
                converter(type,name,bus,Srated,Vdc,P,Q); %adding converter object to converterArr
        end
        %function - cycles through each of the converters in the
        %converterArr array. The harmonic current spectra for each
        %of the converters is solved given the solution of the PCC voltage
        %fund + harmonics in the preceding iteration. The solved injection
        %currents are then used to update spectrum and isource objects
        %corresponding to each of the solved converters.
        function solveConverters(conv_objs)
            for i = 1:conv_objs.numConverters
                solveConverter(conv_objs.converterArr(i)); %the injection currents of the converter referenced in the ith converterArr element are solved.
            end     
        end
    end
end
        