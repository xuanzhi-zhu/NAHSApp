function [inC,inD]=fcn_decide(xi, ref, parameters)

timer=xi(13);

if timer>=0
    inC=1;
else
    inC=0;
end

if timer<=0
    inD=1;
else
    inD=0;
end

end